#!/usr/bin/env python
"""
MefPhot - Multi-Extension FITS Forced Photometry
=================================================

This module performs forced photometry on multi-extension FITS files using
Gaia catalog positions as input sources.

Overview
--------
MefPhot reads FITS files with multiple image extensions, retrieves Gaia sources
in the field of view, and performs aperture photometry on those positions. The
results are saved as FITS tables in the ``TabPhot/`` directory.

Command Line Usage
------------------
.. code-block:: bash

    python MefPhot.py [-h] [-np N] [-r RADIUS] [-b INNER OUTER] [-out ROOT] file1 [file2 ...]

Arguments:
    file1, file2, ...   One or more FITS files to process, or a text/tab file 
                        containing a 'Filename' column
    -h                  Display help message and exit
    -np N               Number of parallel processes (default: 8)
    -r RADIUS           Aperture radius in pixels for flux extraction (default: 6)
                        Also sets background annulus to RADIUS+1 and RADIUS+4
    -b INNER OUTER      Inner and outer radius of background annulus in pixels
                        (default: 8, 12)
    -out ROOT           Root name for output tables (default: derived from filename)

Examples
--------
Process a single file with default parameters:

.. code-block:: bash

    python MefPhot.py observation.fits

Process multiple files in parallel with custom aperture:

.. code-block:: bash

    python MefPhot.py -np 16 -r 5 -b 7 10 file1.fits file2.fits

Process files listed in a text file:

.. code-block:: bash

    python MefPhot.py -out my_photometry file_list.txt

Output
------
Results are written to ``TabPhot/`` directory as FITS tables containing:

- Source positions (pixel and sky coordinates)
- Raw and background-subtracted fluxes
- Photometric uncertainties
- FWHM and eccentricity measurements
- Instrumental magnitudes
- Gaia catalog information

Dependencies
------------
- photutils >= 2.3.0 (earlier versions may hang)
- astropy
- numpy
- scipy
- matplotlib
- GaiaCat module (for Gaia catalog access)
- ImageSum module (for image utilities)

Notes
-----
- The Gaia catalog file ``Gaia_MagClouds.fits`` must be present either locally
  or in the ``kred/xdata`` directory
- Processing time is approximately 8 minutes per MEF file on an M1 Mac
- Background is estimated using sigma-clipped statistics in an annulus
- Sources outside image boundaries (with margin) are automatically excluded

History
-------
- 2025-11-28: Initial coding (KSL)
- 2025-12-05: Tested with photutils 2.3.0 (KSL)

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
import ImageSum


def random_rows(tab, nrows, seed=None):
    """
    Randomly select rows from an Astropy Table without duplicates.

    Parameters
    ----------
    tab : astropy.table.Table
        Input table from which to select rows.
    nrows : int
        Number of rows to randomly select. Must be <= len(tab).
    seed : int, optional
        Random seed for reproducibility. Default is None.

    Returns
    -------
    astropy.table.Table
        Table containing the randomly selected rows.

    Notes
    -----
    If nrows exceeds the table length, returns the full table with a warning.
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

    Attempts to read a table using multiple format detection strategies.

    Parameters
    ----------
    filename : str
        Path to the table file.

    Returns
    -------
    astropy.table.Table
        The loaded table.

    Raises
    ------
    IOError
        If the file does not exist or cannot be read in any supported format.

    Notes
    -----
    Tries FITS format first, then falls back to ASCII detection.
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
                         rstar=4, b_in=4, b_out=8):
    """
    Perform forced aperture photometry at specified sky positions.

    Extracts photometry for sources with known positions (typically from Gaia)
    using aperture photometry with local background subtraction.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file containing image data. Default is 'LMC_c48_T08.r.t060.fits'.
    image_ext : int, optional
        FITS extension number containing the image. Default is 1.
    object_file : str, optional
        Path to file containing source positions (RA, Dec columns required).
        Default is 'objects.txt'.
    nrows_max : int, optional
        Maximum number of sources to process. If -1, process all sources.
        Default is -1.
    rstar : float, optional
        Aperture radius in pixels for source extraction. Default is 4.
    b_in : float, optional
        Inner radius of background annulus in pixels. Default is 4.
    b_out : float, optional
        Outer radius of background annulus in pixels. Default is 8.

    Returns
    -------
    astropy.table.Table or str
        Table containing photometry results with columns:
        
        - id: Source identifier
        - xcenter, ycenter: Pixel coordinates
        - Raw: Raw aperture sum
        - Bkg: Total background in aperture
        - BkgMean: Mean background level
        - BkgStd: Background standard deviation
        - Net: Background-subtracted flux
        - ErrNet: Flux uncertainty
        - FWHM: Full-width at half-maximum
        - Eccentricity: Source eccentricity
        - Max, Min: Maximum and minimum pixel values in aperture
        - phot_mag: Instrumental magnitude (negative if Net < 0)
        - phot_mag_raw: Magnitude from raw flux
        - Original columns from object_file (RA, Dec, Gaia data, etc.)
        
        Returns 'Error' string if file cannot be opened.

    Notes
    -----
    - Background is estimated using sigma-clipped statistics in an annulus
    - Sources outside image boundaries are automatically excluded
    - FWHM is calculated on locally background-subtracted cutouts
    - Magnitude zero point is set to 28
    - Coordinate transformation failures are handled gracefully

    Examples
    --------
    >>> phot = do_forced_photometry('image.fits', image_ext=1,
    ...                             object_file='gaia_sources.fits',
    ...                             rstar=5, b_in=7, b_out=12)
    >>> print(phot['Net', 'ErrNet', 'phot_mag'])
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
        good = ~sources['R'].mask
        sources = sources[good]

    coords = SkyCoord(ra=sources['RA']*u.deg, dec=sources['Dec']*u.deg)

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

    # Format output
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'

    # Merge with original source table
    phot_table = hstack([phot_table, sources])

    # Add source names if not present
    if 'Source_name' not in phot_table.colnames:
        names = [f'x{one["id"]:05d}' for one in phot_table]
        phot_table['Source_name'] = names

    return phot_table


def do_one(filename='foo.fits', outroot='', nrows_max=-1,
           rstar=4, b_in=4, b_out=8):
    """
    Process a single multi-extension FITS file.

    Performs forced photometry on all image extensions in a FITS file using
    Gaia catalog sources.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file to process. Default is 'foo.fits'.
    outroot : str, optional
        Root name for output table. If empty, derived from filename.
        Output written to ``TabPhot/{outroot}.fits``. Default is ''.
    nrows_max : int, optional
        Maximum sources per extension. If -1, process all. Default is -1.
    rstar : float, optional
        Aperture radius in pixels. Default is 4.
    b_in : float, optional
        Inner background annulus radius in pixels. Default is 4.
    b_out : float, optional
        Outer background annulus radius in pixels. Default is 8.

    Returns
    -------
    astropy.table.Table
        Combined photometry table from all extensions with additional columns:
        
        - EXT: FITS extension number
        - CCD: CCD/detector name
        - Filter: Filter name from header
        - Exptime: Exposure time
        - Filename: Input filename

    Raises
    ------
    IOError
        If the file cannot be opened.

    Notes
    -----
    - Automatically retrieves Gaia sources for each extension's field
    - Output table saved to ``TabPhot/`` directory
    - Metadata includes processing timestamp and parameters
    - Filter keyword extracted from header or filename

    Examples
    --------
    >>> phot = do_one('observation.fits', rstar=5, b_in=7, b_out=12)
    >>> print(f"Processed {len(phot)} sources")
    """
    try:
        x = fits.open(filename)
    except:
        print(f'Could not locate {filename}')
        raise IOError

    print(f'do_one: Starting {filename} with radius {rstar:.1f} '
          f'and annulus {b_in:.1f} {b_out:.1f}')

    xexptime = x['PRIMARY'].header['EXPTIME']

    try:
        magzero=x['PRIMARY'].header['MAGZERO']
    except:
        magzero=28.

    # zpt=x['PRIMARY'].header['MAGZERO']
    try:
        srad=x['PRIMARY'].header['RADIUS']
    except:
        srad=-99.
    try:
        ssee=x['PRIMARY'].header['SEEING']
    except:
        ssee=-99.
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
        
        gaia_file = get_gaia(ra, dec, size)
        phot_table = do_forced_photometry(filename, one_extension, gaia_file,
                                          nrows_max, rstar, b_in, b_out)
        phot_table['EXT'] = one_extension
        phot_table['CCD'] = image_extensions['NAME'][i]
        phot_tables.append(phot_table)

    phot = vstack(phot_tables, metadata_conflicts='silent')
    phot['Filter'] = xfilter
    phot['Exptime'] = xexptime
    phot['MAGZERO']=magzero
    phot['SEEING']=ssee
    phot['Star_rad']=srad
    phot['Filename'] = filename

    # Write output
    os.makedirs('TabPhot', exist_ok=True)
    if outroot == '':
        outroot = filename.split('/')[-1]
        outroot = outroot.replace('.fz', '').replace('.fits', '')
    outfile = f'TabPhot/{outroot}.fits'

    # Store metadata
    now = Time.now()
    phot.meta['DATE'] = now.isot
    phot.meta['FILE'] = filename
    phot.meta['RADIUS'] = rstar
    phot.meta['B_IN'] = b_in
    phot.meta['B_OUT'] = b_out

    phot.write(outfile, format='fits', overwrite=True)
    print(f'do_one: Completed {filename} and written to {outfile}')
    return phot


def _safe_do_one_with_index(args):
    """
    Wrapper for parallel processing with error handling.

    Parameters
    ----------
    args : tuple
        (index, filename, outroot, nrows_max, rstar, b_in, b_out)

    Returns
    -------
    tuple
        (filename, success, error_msg, traceback) where success is bool
    """
    index, filename, outroot, nrows_max, rstar, b_in, b_out = args
    try:
        numbered_outroot = f"{outroot}_{index:03d}" if outroot else ''
        do_one(filename, outroot=numbered_outroot, nrows_max=nrows_max,
               rstar=rstar, b_in=b_in, b_out=b_out)
        return (filename, True, None)
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        tb = traceback.format_exc()
        return (filename, False, error_msg, tb)


def do_many(filenames, outroot='', nrows_max=-1, rstar=4, b_in=4, b_out=8,
            n_processes=None, logfile=None, verbose_errors=False):
    """
    Process multiple FITS files in parallel.

    Distributes photometry tasks across multiple CPU cores for efficient
    batch processing.

    Parameters
    ----------
    filenames : list of str
        List of FITS files to process.
    outroot : str, optional
        Output root name. Each file gets ``outroot_NNN`` suffix. Default is ''.
    nrows_max : int, optional
        Maximum sources per extension. -1 for all. Default is -1.
    rstar : float, optional
        Aperture radius in pixels. Default is 4.
    b_in : float, optional
        Inner background annulus radius. Default is 4.
    b_out : float, optional
        Outer background annulus radius. Default is 8.
    n_processes : int, optional
        Number of parallel processes. If None, uses CPU count - 1. Default is None.
    logfile : str, optional
        Path to write error log. If None, prints to screen only. Default is None.
    verbose_errors : bool, optional
        If True, print full tracebacks for errors. Default is False.

    Returns
    -------
    list of tuple
        List of (filename, error_message, traceback) for failed files.

    Examples
    --------
    >>> files = ['obs1.fits', 'obs2.fits', 'obs3.fits']
    >>> failed = do_many(files, n_processes=4, rstar=5, logfile='errors.log')
    >>> if failed:
    ...     print(f"{len(failed)} files failed")
    """
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)

    args_list = [(i, fname, outroot, nrows_max, rstar, b_in, b_out)
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

    Parameters
    ----------
    argv : list
        Command line arguments (typically sys.argv)

    Examples
    --------
    Command line usage::

        python MefPhot.py -np 16 -r 5 -b 7 10 -out results file1.fits file2.fits
    """
    filenames = []
    np_proc = 8
    root = ''
    rstar = 6
    nrows_max = -1
    b_in = 8
    b_out = 12

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
          f'and background annulus of {b_in:.1f} {b_out:.1f}')

    if rstar > b_in or b_in > b_out:
        print('UNPHYSICAL limits for photometry')
        return

    if len(filenames) == 1 or np_proc < 2:
        for one_file in filenames:
            do_one(filename=one_file, outroot=root, nrows_max=nrows_max,
                   rstar=rstar, b_in=b_in, b_out=b_out)
        return

    do_many(filenames, outroot=root, nrows_max=nrows_max, rstar=rstar,
            b_in=b_in, b_out=b_out, n_processes=np_proc, logfile=None)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
