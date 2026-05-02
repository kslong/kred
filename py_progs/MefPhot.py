#!/usr/bin/env python
"""MefPhot - Multi-Extension FITS Forced Photometry

Space Telescope Science Institute

Command Line Usage
------------------

::

    MefPhot.py [-h] [-np N] [-r RADIUS] [-b INNER OUTER] [-cat gaia|smash] [-out ROOT] file1 [file2 ...]

**Required Arguments:**

file1, file2, ...
    One or more FITS files to process, or a text/tab file containing
    a 'Filename' column

**Optional Arguments:**

-h
    Display help message and exit

-np N
    Number of parallel processes (default: 8)

-r RADIUS
    Aperture radius in pixels for flux extraction (default: 6).
    Also sets background annulus to RADIUS+1 and RADIUS+4

-b INNER OUTER
    Inner and outer radius of background annulus in pixels (default: 8, 12)

-cat gaia|smash
    Reference catalog for source positions (default: gaia).
    Use ``smash`` for Magellanic Cloud fields to calibrate against the
    DECam-native SMASH DR2 photometric system.

-out ROOT
    Root name for output tables (default: derived from input filename)

Output
------
Results are written to ``TabPhot/`` as two-extension FITS files.  The
output filename encodes the reference catalog so that Gaia and SMASH
runs on the same input do not overwrite each other::

    TabPhot/<root>.gaia.fits    (default)
    TabPhot/<root>.smash.fits   (with -cat smash)

Extension 0 carries the primary header copied from the input MEF file.
Extension 1 is a FITS table with header keywords DATE, FILE, RADIUS,
B_IN, B_OUT, and CATALOG, and columns including:

* Source positions (pixel and sky coordinates)
* Raw and background-subtracted fluxes with uncertainties
* FWHM and eccentricity measurements
* Instrumental magnitudes (zero point = 28)
* Reference catalog photometry (RA, Dec, G, R, ...)
* EXT, CCD, Filter, Exptime, MAGZERO, SEEING, Filename
* Catalog: 'Gaia' or 'SMASH'

Examples
--------
Process with Gaia (default)::

    MefPhot.py observation.fits
    MefPhot.py -np 16 -r 5 -b 7 10 file1.fits file2.fits

Process using SMASH catalog (Magellanic Cloud fields)::

    MefPhot.py -cat smash observation.fits

Process files listed in a text file::

    MefPhot.py -out my_photometry file_list.txt

Dependencies
------------

* photutils >= 2.3.0 (earlier versions may hang)
* astropy, numpy, scipy, matplotlib
* GaiaCat module (for Gaia DR3 catalog access)
* Smash module (for SMASH DR2 catalog access)
* ImageSum module (for WCS and image utilities)

Notes
-----

* For Gaia: ``Gaia_MagClouds.fits`` must be present locally or in
  ``$KRED/xdata/``
* For SMASH: uses ``Smash_MagClouds.fits`` if present locally or in
  ``$KRED/xdata/``; falls back to a live NOAO Data Lab query (requires
  the ``dl`` package) if the file is absent.  SMASH only covers the
  Magellanic Cloud footprint.
* **Parallel SMASH memory**: ``Smash_MagClouds.fits`` is ~15 GB and each
  worker process loads it independently.  Before spawning workers,
  ``do_many`` checks whether ``-np N`` would exceed available RAM and
  prompts for a lower value if so.  Pre-populating the per-tile cache
  with a serial run (``-np 1``) first avoids the issue entirely on
  subsequent runs.
* Processing time is approximately 8 minutes per MEF file on an M1 Mac
* Background is estimated using sigma-clipped statistics in an annulus
* Sources outside image boundaries (with margin) are automatically excluded

Version History
---------------

2025-11-28 ksl
    Initial coding

2025-12-05 ksl
    Tested with photutils 2.3.0

2026-04-14 ksl
    Added SMASH DR2 as an alternative reference catalog (-cat smash).
    Output filenames now include catalog suffix (.gaia.fits / .smash.fits).
    MEF primary header copied to extension 0 of output file.
    Extraction parameters written to extension 1 header.

Author
------
Space Telescope Science Institute

.. moduleauthor:: KSL
"""


import sys
import os
import numpy as np
from astropy.io import fits, ascii
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
from photutils.aperture import (aperture_photometry, CircularAperture, 
                                 CircularAnnulus, ApertureStats)
from astropy.stats import SigmaClip
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.table import Table, join, hstack, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import timeit
import time
from astropy.time import Time
import pathlib
import os.path as path
import requests
from scipy.spatial import KDTree
from astropy.wcs import NoConvergence
from astropy.wcs._wcs import InvalidCoordinateError
from multiprocessing import Pool
from tqdm import tqdm
import multiprocessing as mp
import traceback

from GaiaCat import get_gaia
from Smash import get_smash
import ImageSum


def get_total_memory():
    """Return total installed physical RAM in bytes, cross-platform.

    Returns
    -------
    int or None
        Total RAM in bytes, or None if it cannot be determined.
    """
    try:
        import psutil
        return psutil.virtual_memory().total
    except ImportError:
        pass

    import platform
    system = platform.system()

    if system == 'Linux':
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        return int(line.split()[1]) * 1024   # kB -> bytes
        except OSError:
            pass

    elif system == 'Darwin':
        try:
            import subprocess
            result = subprocess.run(['sysctl', '-n', 'hw.memsize'],
                                    capture_output=True, text=True, check=True)
            return int(result.stdout.strip())
        except Exception:
            pass

    return None


def get_available_memory():
    """Return available (free + reclaimable) memory in bytes, cross-platform.

    Tries psutil first (works on Linux and macOS), then falls back to
    platform-specific methods so the function works without psutil.

    Returns
    -------
    int or None
        Available memory in bytes, or None if it cannot be determined.
    """
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        pass

    import platform
    system = platform.system()

    if system == 'Linux':
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemAvailable:'):
                        return int(line.split()[1]) * 1024   # kB -> bytes
        except OSError:
            pass

    elif system == 'Darwin':
        try:
            import subprocess
            result = subprocess.run(['vm_stat'], capture_output=True, text=True, check=True)
            page_size = 4096
            free = inactive = speculative = 0
            for line in result.stdout.splitlines():
                if 'page size of' in line:
                    page_size = int(line.split('page size of')[1].split()[0])
                elif line.startswith('Pages free:'):
                    free = int(line.split(':')[1].strip().rstrip('.'))
                elif line.startswith('Pages inactive:'):
                    inactive = int(line.split(':')[1].strip().rstrip('.'))
                elif line.startswith('Pages speculative:'):
                    speculative = int(line.split(':')[1].strip().rstrip('.'))
            return (free + inactive + speculative) * page_size
        except Exception:
            pass

    return None


def check_smash_memory(n_processes, catalog, n_files=None):
    """Warn and prompt when launching multiple SMASH workers would exhaust RAM.

    Each worker process independently loads the pre-assembled SMASH catalog
    into memory.  The actual number of workers spawned is
    ``min(n_processes, n_files)``, so both values are taken into account.
    Always prints a one-line memory summary when the pre-assembled file is
    present, and prompts the user to lower the process count if the total
    would exceed available RAM.

    Parameters
    ----------
    n_processes : int
        Requested number of parallel worker processes.
    catalog : str
        Catalog name ('gaia' or 'smash').  Returns immediately for 'gaia'.
    n_files : int or None
        Number of files to be processed.  If provided, the effective worker
        count is capped at ``min(n_processes, n_files)``.

    Returns
    -------
    int
        The number of processes to actually use (may be adjusted by the user).
    """
    if catalog != 'smash':
        return n_processes

    # Locate the pre-assembled file (mirrors the search in Smash.get_smash_from_file)
    smash_file = None
    kred = os.environ.get('KRED', '')
    for candidate in ['Smash_MagClouds.fits',
                       os.path.join(kred, 'xdata', 'Smash_MagClouds.fits') if kred else None]:
        if candidate and os.path.isfile(candidate):
            smash_file = candidate
            break

    if smash_file is None:
        print('check_smash_memory: pre-assembled file not found; archive path will be used.',
              flush=True)
        return n_processes

    # Effective workers = min(requested, files to process)
    effective = n_processes if n_files is None else min(n_processes, n_files)

    # Workers that hit an already-cached per-tile file skip loading the big table.
    # Count existing Smash/Smash.*.fits files as a proxy for already-cached tiles.
    import glob as _glob
    n_cached   = len(_glob.glob(os.path.join('Smash', 'Smash.*.fits')))
    n_uncached = max(0, n_files - n_cached) if n_files is not None else effective
    # At most `effective` workers run simultaneously; only loaders need the big file.
    loaders = min(effective, n_uncached)

    # Per-worker overhead: Python runtime + astropy imports + FITS data in memory.
    OVERHEAD_GB = 1.5
    file_gb     = os.path.getsize(smash_file) / 1e9
    total_gb    = loaders * file_gb + effective * OVERHEAD_GB
    total_bytes = int(total_gb * 1e9)

    total_ram   = get_total_memory()
    # Warn when estimated usage exceeds 60% of total installed RAM
    WARN_FRAC   = 0.6
    warn_bytes  = int(total_ram * WARN_FRAC) if total_ram is not None else None

    if total_ram is not None:
        total_ram_gb = total_ram / 1e9
        warn_gb      = warn_bytes / 1e9
        # Safe process count: workers that fit within the 60% budget
        safe_np  = max(1, int((warn_bytes / 1e9 - effective * OVERHEAD_GB) / file_gb))
        mem_ok   = total_bytes <= warn_bytes
    else:
        total_ram_gb = None
        warn_gb      = None
        safe_np      = 1
        mem_ok       = False    # can't determine, so warn

    # Always show the memory summary so the user can see the check ran
    cache_note = (f'; {n_cached} tile(s) already cached'
                  f', so {loaders} worker(s) will load big file')
    if total_ram_gb is not None:
        mem_str = (f'{total_gb:.1f} GB estimated vs {warn_gb:.1f} GB limit'
                   f' (60% of {total_ram_gb:.1f} GB total RAM)')
    else:
        mem_str = f'{total_gb:.1f} GB estimated; total RAM unknown'
    print(f'SMASH memory check: {mem_str}{cache_note}', flush=True)

    if mem_ok:
        return n_processes

    print(f'WARNING: insufficient RAM for {effective} parallel SMASH worker(s).',
          flush=True)
    print(f'  Suggested: {safe_np} process(es)', flush=True)
    print(flush=True)

    try:
        response = input(f'Enter number of processes to use [{safe_np}]: ').strip()
        chosen = int(response) if response else safe_np
    except (EOFError, ValueError):
        print(f'Non-interactive or invalid input — using {safe_np} process(es).', flush=True)
        chosen = safe_np

    return max(1, chosen)


def random_rows(tab, nrows, seed=None):
    """
    Randomly select rows from an Astropy Table without duplicates.

    Parameters
    ----------
    tab : Table
        Input astropy Table from which to select rows.
    nrows : int
        Number of rows to randomly select. Must be <= len(tab).
    seed : int, optional
        Random seed for reproducibility. Default is None.

    Returns
    -------
    Table
        Astropy Table containing the randomly selected rows.

    Notes
    -----
    If nrows exceeds the table length, returns the full table with a warning
    message printed to stdout.
    
    Uses numpy's default_rng for random selection, which provides better
    statistical properties than the legacy numpy.random functions.

    Examples
    --------
    >>> from astropy.table import Table
    >>> tab = Table({'a': [1, 2, 3, 4, 5], 'b': [10, 20, 30, 40, 50]})
    >>> subset = random_rows(tab, 3, seed=42)
    >>> print(len(subset))
    3
    """
    if nrows > len(tab):
        print("Requested more rows than available in table")
        return tab

    rng = np.random.default_rng(seed)
    indices = rng.choice(len(tab), size=nrows, replace=False)
    return tab[indices]


def read_table(filename):
    """
    Generic table reader supporting FITS and ASCII formats.

    Attempts to read a table using multiple format detection strategies,
    trying FITS first and falling back to ASCII.

    Parameters
    ----------
    filename : str
        Path to the table file.

    Returns
    -------
    Table
        The loaded astropy Table.

    Raises
    ------
    IOError
        If the file does not exist or cannot be read in any supported format.

    Notes
    -----
    The function attempts reading in this order:
    
    1. As a FITS table using Table.read()
    2. As an ASCII table using ascii.read()
    
    This allows for flexible input formats without requiring the user to
    specify the format explicitly.

    Examples
    --------
    >>> tab = read_table('sources.fits')
    >>> tab = read_table('sources.txt')
    >>> tab = read_table('sources.csv')
    """
    if not os.path.isfile(filename):
        raise IOError(f'read_table: {filename} does not appear to exist')

    try:
        xtable = Table.read(filename)
    except:
        try:
            xtable = ascii.read(filename)
        except:
            raise IOError(f'read_table: {filename} exists, but could not be read')
    return xtable


def do_forced_photometry(filename='LMC_c48_T08.r.t060.fits', image_ext=1,
                         object_file='objects.txt', nrows_max=-1,
                         rstar=6, b_in=8, b_out=12, add_psf_metrics=True):
    """
    Perform forced aperture photometry at specified sky positions.

    Extracts photometry for sources with known positions (typically from Gaia)
    using aperture photometry with local background subtraction.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file containing image data. Default: 'LMC_c48_T08.r.t060.fits'.
    image_ext : int, optional
        FITS extension number containing the image. Default: 1.
    object_file : str, optional
        Path to file containing source positions (RA, Dec columns required).
        Default: 'objects.txt'.
    nrows_max : int, optional
        Maximum number of sources to process. If -1, process all sources.
        Default: -1.
    rstar : float, optional
        Aperture radius in pixels for source extraction. Default: 6.
    b_in : float, optional
        Inner radius of background annulus in pixels. Default: 8.
    b_out : float, optional
        Outer radius of background annulus in pixels. Default: 12.
    add_psf_metrics : bool, optional
        If True, adds columns useful for PSF star selection:
        SNR, Concentration, BkgContam. Default: True.

    Returns
    -------
    Table or str
        Astropy Table containing photometry results with columns:
        
        * id : int - Source identifier
        * xcenter, ycenter : float - Pixel coordinates
        * Raw : float - Raw aperture sum
        * Bkg : float - Total background in aperture
        * BkgMean : float - Mean background level
        * BkgStd : float - Background standard deviation
        * Net : float - Background-subtracted flux
        * ErrNet : float - Flux uncertainty
        * FWHM : float - Full-width at half-maximum
        * Eccentricity : float - Source eccentricity (0=circular, 1=linear)
        * Max, Min : float - Maximum and minimum pixel values in aperture
        * phot_mag : float - Instrumental magnitude (negative if Net < 0)
        * phot_mag_raw : float - Magnitude from raw flux
        * SNR : float - Signal-to-noise ratio (if add_psf_metrics=True)
        * Concentration : float - Peak/mean flux density (if add_psf_metrics=True)
        * BkgContam : float - Background contamination normalized to median (if add_psf_metrics=True)
        * Original columns from object_file (RA, Dec, Gaia data, etc.)

        Returns 'Error' string if file cannot be opened.

    Notes
    -----
    **Background Estimation:**
    
    Background is estimated using sigma-clipped statistics in an annulus
    defined by b_in and b_out. The sigma clipping uses 3-sigma rejection
    with up to 10 iterations.
    
    **Source Selection:**
    
    Sources are automatically filtered to exclude those:
    
    * Outside image boundaries (with margin = 2 * b_out)
    * With non-finite pixel coordinates
    * That fail coordinate transformation
    
    **FWHM Calculation:**
    
    FWHM is calculated on locally background-subtracted cutouts around each
    source using ApertureStats.
    
    **Magnitude System:**
    
    Instrumental magnitudes use a zero point of 28:
    
    mag = 28 - 2.5 * log10(flux)
    
    Negative fluxes result in negative magnitudes.
    
    **Error Estimation:**
    
    Flux uncertainty includes:
    
    * Poisson noise from the source
    * Background variance in the aperture
    * Uncertainty in background estimation

    Examples
    --------
    >>> phot = do_forced_photometry('image.fits', image_ext=1,
    ...                             object_file='gaia_sources.fits',
    ...                             rstar=5, b_in=7, b_out=12)
    >>> print(phot['Net', 'ErrNet', 'phot_mag'])
    >>> bright = phot[phot['Net'] > 1000]
    """
    try:
        x = fits.open(filename)
    except:
        print(f'Error: get_photometry: could not open {filename}')
        return 'Error'

    image_wcs = WCS(x[image_ext].header)
    image = x[image_ext].data
    
    # Create mask for invalid pixels
    image_mask = (image == 0) | ~np.isfinite(image)
    NAXIS1 = x[image_ext].header['NAXIS1']
    NAXIS2 = x[image_ext].header['NAXIS2']

    # Load source catalog
    sources = read_table(object_file)
    if 'G' in sources.colnames:
        good = np.isfinite(sources['R'])
        sources = sources[good]

    coords = SkyCoord(ra=np.array(sources['RA'])*u.deg, dec=np.array(sources['Dec'])*u.deg)

    # Transform coordinates to pixel space
    sources['xcentroid'] = np.nan
    sources['ycentroid'] = np.nan

    try:
        x_pix, y_pix = image_wcs.world_to_pixel(coords)
        sources['xcentroid'] = x_pix
        sources['ycentroid'] = y_pix
    except (NoConvergence, InvalidCoordinateError) as e:
        if isinstance(e, NoConvergence):
            sources['xcentroid'] = e.best_solution[0]
            sources['ycentroid'] = e.best_solution[1]
            print(f"Warning: {len(e.divergent)} coordinates failed to converge")
        else:
            print(f"Warning: Severe coordinate transformation error - skipping bad coordinates")

    # Filter sources within image boundaries
    margin = 2 * b_out
    mask = (
        np.isfinite(sources['xcentroid']) &
        np.isfinite(sources['ycentroid']) &
        (sources['xcentroid'] >= margin) &
        (sources['xcentroid'] < NAXIS1 - margin) &
        (sources['ycentroid'] >= margin) &
        (sources['ycentroid'] < NAXIS2 - margin)
    )
    sources = sources[mask]

    # Optionally limit number of sources
    if nrows_max > 0 and len(sources) > nrows_max:
        sources = random_rows(sources, nrows=nrows_max, seed=None)
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    # Define apertures
    apertures = CircularAperture(positions, r=rstar)
    annulus_apertures = CircularAnnulus(positions, r_in=b_in, r_out=b_out)

    # Calculate background statistics with sigma clipping
    sigclip = SigmaClip(sigma=3, maxiters=10)
    bkg_stats = ApertureStats(image, annulus_apertures, sigma_clip=sigclip, mask=image_mask)
    aper_stats = ApertureStats(image, apertures, sigma_clip=None, mask=image_mask)
    
    # Calculate background-subtracted flux
    bkg_mean = bkg_stats.mean
    n_aper_pixels = aper_stats.sum_aper_area.value
    total_background = bkg_mean * n_aper_pixels
    net = aper_stats.sum - total_background

    # Error estimation
    bkg_std_per_pixel = bkg_stats.std
    n_bkg_pixels = bkg_stats.sum_aper_area.value
    error = np.sqrt(
        np.abs(net) +
        n_aper_pixels * bkg_std_per_pixel**2 +
        n_aper_pixels**2 * bkg_std_per_pixel**2 / n_bkg_pixels
    )

    # Create output table
    phot_table = Table()
    phot_table['id'] = np.arange(len(positions))
    phot_table['xcenter'] = positions[:, 0]
    phot_table['ycenter'] = positions[:, 1]
    phot_table['Raw'] = aper_stats.sum
    phot_table['Bkg'] = total_background
    phot_table['BkgMean'] = bkg_mean
    phot_table['BkgStd'] = bkg_std_per_pixel
    phot_table['Net'] = net
    phot_table['ErrNet'] = error
    
    # Calculate FWHM on local background-subtracted cutouts
    phot_table['FWHM'] = np.nan
    phot_table['Eccentricity'] = np.nan

    for i in range(len(positions)):
        x_int, y_int = int(positions[i, 0]), int(positions[i, 1])
        cutout_size = int(2 * b_out) + 10
        y_min = max(0, y_int - cutout_size)
        y_max = min(NAXIS2, y_int + cutout_size)
        x_min = max(0, x_int - cutout_size)
        x_max = min(NAXIS1, x_int + cutout_size)
    
        if y_max > y_min and x_max > x_min:
            cutout = image[y_min:y_max, x_min:x_max] - bkg_mean[i]
            cutout_pos = [(positions[i, 0] - x_min, positions[i, 1] - y_min)]
            cutout_aper = CircularAperture(cutout_pos, r=rstar)
            cutout_stats = ApertureStats(cutout, cutout_aper, sigma_clip=None)
        
            if np.isfinite(cutout_stats.fwhm.value):
                phot_table['FWHM'][i] = cutout_stats.fwhm.value
            if np.isfinite(cutout_stats.eccentricity):
                phot_table['Eccentricity'][i] = cutout_stats.eccentricity

    phot_table['Max'] = aper_stats.max
    phot_table['Min'] = aper_stats.min

    # Magnitude calculation (zero point = 28)
    phot_table['phot_mag'] = 28 - 2.5 * np.log10(np.abs(phot_table['Net']))
    phot_table['phot_mag_raw'] = 28 - 2.5 * np.log10(np.abs(phot_table['Raw']))

    # Make magnitudes negative for negative fluxes
    phot_table['phot_mag'] = np.select([phot_table['Net'] > 0],
                                       [phot_table['phot_mag']],
                                       default=-phot_table['phot_mag'])
    phot_table['phot_mag_raw'] = np.select([phot_table['Net'] > 0],
                                           [phot_table['phot_mag_raw']],
                                           default=-phot_table['phot_mag_raw'])

    # Add PSF quality metrics if requested
    if add_psf_metrics:
        # Signal-to-noise ratio
        phot_table['SNR'] = np.abs(phot_table['Net']) / phot_table['ErrNet']

        # Light concentration (peak/mean flux density)
        mean_flux_density = phot_table['Net'] / n_aper_pixels
        phot_table['Concentration'] = phot_table['Max'] / mean_flux_density

        # Background contamination (normalized to median)
        median_bkg_std = np.nanmedian(phot_table['BkgStd'])
        phot_table['BkgContam'] = phot_table['BkgStd'] / median_bkg_std

    # Format output
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'

    # Merge with original source table
    phot_table = hstack([phot_table, sources])

    # Add source names if not present
    if 'Source_name' not in phot_table.colnames:
        id_col = 'id' if 'id' in phot_table.colnames else 'id_1'
        if id_col in phot_table.colnames:
            names = [f'x{one[id_col]:05d}' for one in phot_table]
        else:
            names = [f'x{i:05d}' for i in range(len(phot_table))]
        phot_table['Source_name'] = names

    return phot_table


def do_one(filename='foo.fits', outroot='', nrows_max=-1,
           rstar=6, b_in=8, b_out=12, catalog='gaia'):
    """
    Process a single multi-extension FITS file.

    Performs forced photometry on all image extensions in a FITS file using
    Gaia catalog sources. Automatically retrieves Gaia sources for each
    extension's field of view.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file to process. Default: 'foo.fits'.
    outroot : str, optional
        Root name for output table. If empty, derived from filename.
        Output written to ``TabPhot/{outroot}.fits``. Default: ''.
    nrows_max : int, optional
        Maximum sources per extension. If -1, process all. Default: -1.
    rstar : float, optional
        Aperture radius in pixels. Default: 6.
    b_in : float, optional
        Inner background annulus radius in pixels. Default: 8.
    b_out : float, optional
        Outer background annulus radius in pixels. Default: 12.
    catalog : str, optional
        Reference catalog to use: ``'gaia'`` (default) or ``'smash'``.
        SMASH is only available over the Magellanic Cloud footprint but
        uses the same DECam photometric system as the images, giving a
        smaller color term in ZeroCalc.

    Returns
    -------
    Table
        Combined photometry table from all extensions with additional columns:

        * EXT : int - FITS extension number
        * CCD : str - CCD/detector name
        * Filter : str - Filter name from header
        * Exptime : float - Exposure time in seconds
        * MAGZERO : float - Magnitude zero point
        * SEEING : float - Seeing FWHM in arcseconds
        * Star_rad : float - Star radius parameter
        * Filename : str - Input filename
        * Catalog : str - Reference catalog used ('Gaia' or 'SMASH')

    Raises
    ------
    IOError
        If the file cannot be opened.

    Notes
    -----
    **Processing Steps:**
    
    1. Open FITS file and read primary header
    2. Identify all image extensions
    3. For each extension:
       
       - Determine field center and size
       - Retrieve Gaia sources in field
       - Perform forced photometry
       - Add extension metadata
    
    4. Stack results from all extensions
    5. Write to TabPhot/ directory
    
    **Output Directory:**
    
    Automatically creates ``TabPhot/`` directory if it doesn't exist.
    
    **Metadata:**
    
    Output table metadata includes processing timestamp, input filename,
    and photometry parameters (radius, background annulus).
    
    **Filter Extraction:**
    
    Tries to read FILTER keyword from header. If missing, extracts from
    filename (assumes format: ...{filter}.t{number}.fits).

    Examples
    --------
    >>> # Process with default parameters
    >>> phot = do_one('observation.fits')
    >>> print(f"Processed {len(phot)} sources")
    
    >>> # Custom aperture settings
    >>> phot = do_one('observation.fits', rstar=5, b_in=7, b_out=12)
    
    >>> # Limit sources for testing
    >>> phot = do_one('observation.fits', nrows_max=100)
    """
    try:
        x = fits.open(filename)
    except:
        print(f'Could not locate {filename}')
        raise IOError

    catalog = catalog.lower()
    print(f'do_one: Starting {filename} with radius {rstar:.1f} '
          f'and annulus {b_in:.1f} {b_out:.1f} using {catalog} catalog')

    xexptime = x['PRIMARY'].header['EXPTIME']

    try:
        magzero = x['PRIMARY'].header['MAGZERO']
    except:
        magzero = 28.

    try:
        srad = x['PRIMARY'].header['RADIUS']
    except:
        srad = -99.
        
    try:
        ssee = x['PRIMARY'].header['SEEING']
    except:
        ssee = -99.
        
    try:
        xfilter = x['PRIMARY'].header['FILTER']
    except:
        words = filename.split('.')
        xfilter = words[-3]
        print(f'Filter keyword is missing. Setting to {xfilter} for {filename}')
        
    image_extensions = ImageSum.list_image_extensions(filename)
    phot_tables = []

    for i, one_extension in enumerate(np.array(image_extensions['EXT'])):
        info = ImageSum.get_image_center_and_size_from_header(x[one_extension].header)
        ra = info['center_ra']
        dec = info['center_dec']
        width = info['width_deg']
        height = info['height_deg']
        size = np.sqrt(width*width + height*height) / 2.
        
        if catalog == 'smash':
            cat_file = get_smash(ra, dec, size)
        else:
            cat_file = get_gaia(ra, dec, size)
        phot_table = do_forced_photometry(filename, one_extension, cat_file,
                                          nrows_max, rstar, b_in, b_out)
        phot_table['EXT'] = one_extension
        phot_table['CCD'] = image_extensions['NAME'][i]
        phot_tables.append(phot_table)

    phot = vstack(phot_tables, metadata_conflicts='silent')
    phot['Filter'] = xfilter
    phot['Exptime'] = xexptime
    phot['MAGZERO'] = magzero
    phot['SEEING'] = ssee
    phot['Star_rad'] = srad
    phot['Filename'] = filename
    phot['Catalog'] = 'SMASH' if catalog == 'smash' else 'Gaia'

    # Write output
    os.makedirs('TabPhot', exist_ok=True)
    if outroot == '':
        outroot = filename.split('/')[-1]
        outroot = outroot.replace('.fz', '').replace('.fits', '')
    cat_suffix = '.smash' if catalog == 'smash' else '.gaia'
    outfile = f'TabPhot/{outroot}{cat_suffix}.fits'

    now = Time.now()

    # Extension 0: copy primary header from input MEF file
    primary_hdu = fits.PrimaryHDU(header=x['PRIMARY'].header.copy())

    # Extension 1: photometry table with extraction parameters in header
    table_hdu = fits.table_to_hdu(phot)
    table_hdu.header['DATE'] = (now.isot, 'Processing timestamp')
    table_hdu.header['FILE'] = (filename[:68], 'Input FITS file')
    table_hdu.header['RADIUS'] = (rstar, 'Aperture radius in pixels')
    table_hdu.header['B_IN'] = (b_in, 'Inner background annulus radius in pixels')
    table_hdu.header['B_OUT'] = (b_out, 'Outer background annulus radius in pixels')
    table_hdu.header['CATALOG'] = ('SMASH' if catalog == 'smash' else 'Gaia',
                                   'Reference catalog used for source positions')

    fits.HDUList([primary_hdu, table_hdu]).writeto(outfile, overwrite=True)
    x.close()
    print(f'do_one: Completed {filename} and written to {outfile}')
    return phot


def _safe_do_one_with_index(args):
    """
    Wrapper for parallel processing with error handling.
    
    Internal function used by do_many() for safe parallel execution.

    Parameters
    ----------
    args : tuple
        Tuple of (index, filename, outroot, nrows_max, rstar, b_in, b_out)
        where:
        
        * index : int - File index for output naming
        * filename : str - Path to FITS file
        * outroot : str - Output root name
        * nrows_max : int - Maximum sources to process
        * rstar : float - Aperture radius
        * b_in : float - Inner background radius
        * b_out : float - Outer background radius
        * catalog : str - Reference catalog ('gaia' or 'smash')

    Returns
    -------
    tuple
        Success case: (filename, True, None)
        
        Failure case: (filename, False, error_msg, traceback)
        
        where error_msg is a string describing the error and traceback
        is the full Python traceback as a string.

    Notes
    -----
    This function catches all exceptions to prevent multiprocessing pool
    failures. Exceptions are converted to string messages for reporting.
    """
    index, filename, outroot, nrows_max, rstar, b_in, b_out, catalog = args
    try:
        numbered_outroot = f"{outroot}_{index:03d}" if outroot else ''
        do_one(filename, outroot=numbered_outroot, nrows_max=nrows_max,
               rstar=rstar, b_in=b_in, b_out=b_out, catalog=catalog)
        return (filename, True, None)
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        tb = traceback.format_exc()
        return (filename, False, error_msg, tb)


def do_many(filenames, outroot='', nrows_max=-1, rstar=6, b_in=8, b_out=12,
            n_processes=None, logfile=None, verbose_errors=False, catalog='gaia'):
    """
    Process multiple FITS files in parallel.

    Distributes photometry tasks across multiple CPU cores for efficient
    batch processing. Includes progress tracking and error reporting.

    Parameters
    ----------
    filenames : list of str
        List of FITS files to process.
    outroot : str, optional
        Output root name. Each file gets ``outroot_NNN`` suffix where NNN
        is a zero-padded index. Default: ''.
    nrows_max : int, optional
        Maximum sources per extension. -1 for all. Default: -1.
    rstar : float, optional
        Aperture radius in pixels. Default: 6.
    b_in : float, optional
        Inner background annulus radius. Default: 8.
    b_out : float, optional
        Outer background annulus radius. Default: 12.
    n_processes : int, optional
        Number of parallel processes. If None, uses CPU count - 1. Default: None.
    logfile : str, optional
        Path to write error log. If None, prints to screen only. Default: None.
    verbose_errors : bool, optional
        If True, print full tracebacks for errors. Default: False.

    Returns
    -------
    list of tuple
        List of (filename, error_message, traceback) tuples for failed files.
        Empty list if all files processed successfully.

    Notes
    -----
    **Progress Display:**
    
    Uses tqdm to display a progress bar showing processing status.
    
    **Error Handling:**
    
    Errors in individual files do not stop the batch. Failed files are
    collected and reported at the end.
    
    **Performance:**
    
    Processing speed scales roughly linearly with number of cores up to
    the number of available CPUs. I/O bottlenecks may limit scaling for
    very fast processors.
    
    **Output Naming:**
    
    Files are numbered sequentially. For example, with outroot='phot',
    outputs will be: phot_000.fits, phot_001.fits, etc.

    Examples
    --------
    >>> # Process all files with 8 cores
    >>> files = ['obs1.fits', 'obs2.fits', 'obs3.fits']
    >>> failed = do_many(files, n_processes=8, rstar=5)
    Processing images: 100%|██████████| 3/3 [00:24<00:00,  8.2s/it]
    Successfully processed all 3 files
    
    >>> # With error logging
    >>> failed = do_many(files, n_processes=4, logfile='errors.log')
    >>> if failed:
    ...     print(f"{len(failed)} files failed - see errors.log")
    """
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)

    args_list = [(i, fname, outroot, nrows_max, rstar, b_in, b_out, catalog)
                 for i, fname in enumerate(filenames)]

    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap(_safe_do_one_with_index, args_list),
                           total=len(filenames),
                           desc="Processing images"))

    # Collect failures
    failed_files = []
    for result in results:
        if len(result) == 4:  # Failure case
            fname, success, error_msg, tb = result
            if not success:
                failed_files.append((fname, error_msg, tb))

    # Report results
    if failed_files:
        print(f"\n{len(failed_files)}/{len(filenames)} files failed to process")
        print("\nErrors:")
        for fname, error_msg, tb in failed_files:
            print(f"\n  {fname}:")
            print(f"    {error_msg}")
            if verbose_errors:
                print("    Full traceback:")
                for line in tb.split('\n'):
                    print(f"      {line}")

        if logfile:
            with open(logfile, 'w') as f:
                f.write(f"{len(failed_files)}/{len(filenames)} files failed\n\n")
                for fname, error_msg, tb in failed_files:
                    f.write(f"{fname}\n")
                    f.write(f"  Error: {error_msg}\n")
                    if verbose_errors:
                        f.write(f"  Traceback:\n")
                        for line in tb.split('\n'):
                            f.write(f"    {line}\n")
                    f.write("\n")
            print(f"\nFailed filenames and errors written to {logfile}")
    else:
        print(f"\nSuccessfully processed all {len(filenames)} files")

    return failed_files


def steer(argv):
    """
    Parse command line arguments and execute photometry.
    
    Main entry point for command-line execution. Handles argument parsing
    and delegates to do_one() or do_many() as appropriate.

    Parameters
    ----------
    argv : list
        Command line arguments (typically sys.argv).

    Returns
    -------
    None
        Results are written to files in TabPhot/ directory.

    Notes
    -----
    **Argument Parsing:**
    
    Uses a simple while loop to parse arguments. Unknown switches cause
    an error message and early return.
    
    **File List Input:**
    
    If an argument contains '.txt' or 'tab', it's treated as a file
    containing a list of FITS files to process (must have 'Filename' column).
    
    **Parallel Processing:**
    
    If more than one file and np > 1, uses do_many() for parallel processing.
    Otherwise uses sequential processing with do_one().
    
    **Parameter Validation:**
    
    Checks that rstar < b_in < b_out before processing. Prints error and
    returns if this constraint is violated.

    Examples
    --------
    From command line::
    
        python MefPhot.py -np 16 -r 5 -b 7 10 obs1.fits obs2.fits
        python MefPhot.py -out results file_list.txt
        python MefPhot.py -h
    """
    filenames = []
    np_proc = 8
    root = ''
    rstar = 6
    nrows_max = -1
    b_in = 8
    b_out = 12
    catalog = 'gaia'

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i][:3] == '-np':
            i += 1
            np_proc = int(argv[i])
        elif argv[i][:4] == '-out':
            i += 1
            root = argv[i]
        elif argv[i][:4] == '-cat':
            i += 1
            catalog = argv[i].lower()
            if catalog not in ('gaia', 'smash'):
                print(f'Error: -cat must be gaia or smash, got {argv[i]}')
                return
        elif argv[i][:2] == '-r':
            i += 1
            rstar = float(argv[i])
            b_in = rstar + 1
            b_out = rstar + 4
        elif argv[i][:2] == '-b':
            i += 1
            b_in = float(argv[i])
            i += 1
            b_out = float(argv[i])
        elif argv[i][0] == '-':
            print('Error: unknown switch:', argv[i])
            return
        elif '.txt' in argv[i] or 'tab' in argv[i]:
            xtab = ascii.read(argv[i])
            filenames = list(xtab['Filename'])
        else:
            filenames.append(argv[i])
        i += 1

    print(f'Starting with {len(filenames)} filenames and rstar of {rstar:.1f} '
          f'and background annulus of {b_in:.1f} {b_out:.1f} using {catalog} catalog')

    if rstar > b_in or b_in > b_out:
        print('UNPHYSICAL limits for photometry')
        return

    np_proc = check_smash_memory(np_proc, catalog, n_files=len(filenames))

    if len(filenames) == 1 or np_proc < 2:
        for one_file in filenames:
            do_one(filename=one_file, outroot=root, nrows_max=nrows_max,
                   rstar=rstar, b_in=b_in, b_out=b_out, catalog=catalog)
        return

    do_many(filenames, outroot=root, nrows_max=nrows_max, rstar=rstar,
            b_in=b_in, b_out=b_out, n_processes=np_proc, logfile=None,
            catalog=catalog)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
