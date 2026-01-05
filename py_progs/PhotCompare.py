#!/usr/bin/env python
# coding: utf-8

"""PhotCompare - Photometric Comparison Tool

Space Telescope Science Institute

Command Line Usage
------------------

::

    PhotCompare.py [-h] [-dir DIRNAME] [-nmax N] [-forced] [-unforced]
                   [-gcat FILE] [-out NAME] file1 file2 ...

**Operating Modes:**

There are two basic modes of operation:

1. **Directory mode** (with -dir): Process all FITS files in the specified
   directory and any subdirectories. Individual file arguments are ignored.

2. **File mode** (without -dir): Process only the specified files.

**Optional Arguments:**

-h
    Print help message and exit

-dir DIRNAME
    Process all files in DIRNAME and subdirectories. Assumes these are
    swarped versions of the original data.

-nmax N
    Limit the number of positions used for forced photometry from the
    Gaia catalog. If nmax < 0, all positions are processed. Default: 30000.

-forced
    Use forced photometry (default mode). Performs photometry at Gaia
    catalog positions.

-unforced
    Search for sources in the image, then cross-match positions to Gaia.
    This is a diagnostic mode useful for checking relative astrometry
    between Gaia and our images.

-gcat FILE
    Use specified Gaia catalog file instead of auto-generating

-out NAME
    Specify output root name for results

file1 file2 ...
    One or more FITS files to process (ignored if -dir is specified)

Processing Modes
**Forced Photometry Mode (default):**

Performs aperture photometry at positions from the Gaia catalog. This is
the standard mode for most applications.

**Unforced Mode:**

Searches for sources in the image using DAOStarFinder, then cross-matches
detected positions to Gaia. Useful for diagnosing astrometric issues.
Search results are stored in ``TabPhot/``.

Output
------
The routine generates:

* **Figures**: Saved to ``Figs_phot/`` directory showing:

  - Magnitude comparisons (Gaia vs DECam)
  - Color-magnitude diagrams
  - Residual plots

* **Tables**: Saved to ``TabPhot/`` directory containing:

  - Photometry results
  - Cross-matched catalogs
  - Source lists (unforced mode)

Performance Notes
The most time-consuming operation is Gaia catalog retrieval. To optimize:

* Catalogs are cached and reused when processing multiple files with the
  same field center and size
* Cached catalogs are stored in a ``GAIA/`` subdirectory
* If all files cover the same region, retrieval happens only once

Examples
--------
Process a single file with forced photometry::

    python PhotCompare.py image.fits

Process all files in a directory::

    python PhotCompare.py -dir DECamSWARP2/SMC_c01

Use unforced mode with limited sources::

    python PhotCompare.py -unforced -nmax 5000 image.fits

Process multiple files with custom Gaia catalog::

    python PhotCompare.py -gcat my_gaia.fits file1.fits file2.fits

Notes
-----

The most time-consuming operation is Gaia catalog retrieval. To optimize:

* Catalogs are cached and reused when processing multiple files with the
  same field center and size
* Cached catalogs are stored in a ``GAIA/`` subdirectory
* If all files cover the same region, retrieval happens only once

Version History
---------------

240318 ksl

    Coding begun

240527 ksl

    Speed up catalog matching with KDTree

251105 ksl

    Split finding sources from doing photometry

251130 ksl

    Starting cleaning

Author
Space Telescope Science Institute
"""


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
from astropy.table import Table, join, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn", force=True)
import pathlib
import os.path as path
import requests
from gaiaxpy import calibrate
from scipy.spatial import KDTree
from astropy.wcs import NoConvergence
from astropy.wcs._wcs import InvalidCoordinateError
from http.client import IncompleteRead

import ImageSum
from GaiaCat import get_gaia


#: Directory suffix for isolating different runs of PhotCompare
XDIR = ''


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
    Table
        The loaded astropy Table.
    
    Raises
    ------
    IOError
        If the file does not exist or cannot be read in any supported format.
    
    Notes
    -----
    Tries FITS format first, then falls back to ASCII detection.
    
    Examples
    --------
    >>> tab = read_table('sources.fits')
    >>> tab = read_table('sources.txt')
    """
    print('XXXX - filename ', filename)

    if not os.path.isfile(filename):
        raise IOError('read_table: %s does not appear to exist' % filename)

    try:
        xtable = Table.read(filename)
    except:
        try:
            xtable = ascii.read(filename)
        except:
            raise IOError('read_table: %s exist, but could not be read' % filename)
    return xtable


def random_rows(tab, nrows, seed=None):
    """
    Randomly select rows from an Astropy Table without duplicates.

    Parameters
    ----------
    tab : Table
        Input astropy Table.
    nrows : int
        Number of rows to randomly select (must be <= len(tab)).
    seed : int, optional
        Random seed for reproducibility. Default: None.

    Returns
    -------
    Table
        Astropy Table containing the randomly selected rows.
        
    Notes
    -----
    If nrows exceeds the table length, returns the full table with a
    warning message.
    
    Examples
    --------
    >>> from astropy.table import Table
    >>> tab = Table({'a': [1, 2, 3, 4, 5]})
    >>> subset = random_rows(tab, 3, seed=42)
    >>> len(subset)
    3
    """
    if nrows > len(tab):
        print("Requested more rows than available in table")
        return tab

    rng = np.random.default_rng(seed)
    indices = rng.choice(len(tab), size=nrows, replace=False)
    return tab[indices]


def unique_rows_within_tol(tab, tol=0.01):
    """
    Return unique rows based on approximate equality within tolerance.
    
    Identifies unique field positions (RA, Dec, Size) within a specified
    tolerance, useful for grouping observations of the same field.

    Parameters
    ----------
    tab : Table
        Table containing columns 'RA', 'Dec', and 'Size' (in degrees).
    tol : float, optional
        Matching tolerance in degrees. Default: 0.01°.

    Returns
    -------
    unique_tab : Table
        New table containing one representative row per unique group.
    mapping : ndarray
        Array of length len(tab) where mapping[i] gives the index in
        unique_tab that row i of the original table maps to.
    
    Notes
    -----
    This function is used to identify files that observe the same field,
    allowing Gaia catalogs to be reused and avoiding redundant downloads.
    
    The tolerance of 0.01° (~36 arcsec) is typically sufficient to
    identify overlapping fields while avoiding false matches.
    
    Examples
    --------
    >>> from astropy.table import Table
    >>> tab = Table({'RA': [10.0, 10.001, 20.0],
    ...              'Dec': [-30.0, -30.001, -40.0],
    ...              'Size': [0.5, 0.5, 0.5]})
    >>> unique, mapping = unique_rows_within_tol(tab, tol=0.01)
    >>> len(unique)
    2
    >>> mapping
    array([0, 0, 1])
    """
    # Stack RA, Dec, Size into a NumPy array
    data = np.vstack([tab['RA'], tab['Dec'], tab['Size']]).T

    # Initialize list of unique rows and mapping array
    unique_indices = []
    mapping = np.full(len(data), -1, dtype=int)
    used = np.zeros(len(data), dtype=bool)

    for i in range(len(data)):
        if used[i]:
            continue

        # Find all rows within tolerance of row i
        diff = np.abs(data - data[i])
        mask = np.all(diff < tol, axis=1)

        # Mark them as used and map them to the current unique group
        used[mask] = True
        group_idx = len(unique_indices)
        mapping[mask] = group_idx

        # Add representative row
        unique_indices.append(i)

    return tab[unique_indices], mapping


def do_fig(xtab, outroot=''):
    """
    Create diagnostic photometry comparison figures.
    
    Generates a 2x2 panel figure comparing Gaia and DECam photometry,
    including magnitude comparisons and residual plots.

    Parameters
    ----------
    xtab : Table
        Cross-matched table containing both Gaia and DECam photometry
        with columns: G, R (Gaia), phot_mag (DECam).
    outroot : str, optional
        Output filename root. Default: ''.

    Returns
    -------
    None
        Figure is saved to ``Figs_phot/`` directory.
    
    Notes
    -----
    **Figure Layout:**
    
    * Panel 1 (top-left): DECam vs Gaia G magnitude
    * Panel 2 (top-right): DECam vs Gaia R magnitude
    * Panel 3 (bottom-left): Residuals vs Gaia G
    * Panel 4 (bottom-right): Residuals vs Gaia R
    
    All panels use G-R color coding (plasma colormap) to show color trends.
    Negative DECam magnitudes (from negative fluxes) are plotted separately.
    
    **Output:**
    
    Saved as PNG to ``Figs_phot{XDIR}/{outroot}.png``
    
    Examples
    --------
    >>> xtab = ascii.read('cross_match.txt')
    >>> do_fig(xtab, outroot='LMC_field1')
    """
    outdir = './Figs_phot%s' % XDIR
    os.makedirs(outdir, exist_ok=True)
    
    plt.figure(1, (9, 8))
    plt.clf()
    
    # Panel 1: DECam vs Gaia G
    plt.subplot(2, 2, 1)
    if 'G' in xtab.colnames:
        sc = plt.scatter(xtab['G'], xtab['phot_mag'], marker='.', alpha=.05, 
                        c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
        sc = plt.scatter(xtab['G'], -xtab['phot_mag'], marker='.', alpha=.05, 
                        c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
        cbar = plt.colorbar(sc)
        cbar.set_label('G-R')
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0)
    else:
        plt.scatter(xtab['G'], xtab['phot_mag'], marker='.', alpha=.05)
        plt.scatter(xtab['G'], -xtab['phot_mag'], marker='.', alpha=.05)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11, 24], [11, 24], 'k-')
    plt.ylim(14, 22)
    plt.xlim(14, 22)

    # Panel 2: DECam vs Gaia R
    plt.subplot(2, 2, 2)
    if 'G' in xtab.colnames:
        sc = plt.scatter(xtab['R'], xtab['phot_mag'], marker='.', alpha=.05, 
                        c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
        sc = plt.scatter(xtab['R'], -xtab['phot_mag'], marker='.', alpha=.05, 
                        c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
        cbar = plt.colorbar(sc)
        cbar.set_label('G-R')
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0)
    else:
        plt.scatter(xtab['R'], xtab['phot_mag'], marker='.', alpha=.05)
        plt.scatter(xtab['R'], -xtab['phot_mag'], marker='.', alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11, 24], [11, 24], 'k-')
    plt.ylim(14, 22)
    plt.xlim(14, 22)

    # Panel 3: Residuals vs G
    plt.subplot(2, 2, 3)
    sc = plt.scatter(xtab['G'], xtab['phot_mag']-xtab['G'], marker='.', alpha=.01, 
                    c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
    sc = plt.scatter(xtab['G'], xtab['phot_mag']+xtab['G'], marker='.', alpha=.01, 
                    c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
    cbar = plt.colorbar(sc)
    cbar.set_label('G-R')
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam - Gaia (mag)')
    plt.plot([11, 24], [0, 0], 'k-')
    plt.ylim(-2, 2)
    plt.xlim(14, 22)

    # Panel 4: Residuals vs R
    plt.subplot(2, 2, 4)
    under = xtab[xtab['phot_mag'] > 0]
    sc = plt.scatter(xtab['R'], xtab['phot_mag']-xtab['R'], marker='.', alpha=.01, 
                    c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
    sc = plt.scatter(xtab['R'], xtab['phot_mag']+xtab['R'], marker='.', alpha=.01, 
                    c=xtab['G']-xtab['R'], cmap='plasma', vmin=-1, vmax=1)
    cbar = plt.colorbar(sc)
    cbar.set_label('G-R')
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam - Gaia (mag)')
    plt.plot([11, 24], [0, 0], 'k-')
    plt.ylim(-2, 2)
    plt.xlim(14, 22)

    plt.suptitle(outroot)
    plt.tight_layout()
    plt.savefig('%s/%s.png' % (outdir, outroot))


def do_fig_diff(xtab, outroot):
    """
    Create simplified residual-only comparison figures.
    
    Generates a 1x2 panel figure showing photometry residuals between
    Gaia and DECam in G and R bands.

    Parameters
    ----------
    xtab : Table
        Cross-matched table with Gaia and DECam photometry.
    outroot : str
        Output filename root.

    Returns
    -------
    None
        Figure is saved to ``Figs_phot/`` directory.
    
    Notes
    -----
    This is a simplified version of do_fig() showing only residuals,
    useful for quick diagnostic checks. Reports number of positive vs
    negative flux detections.
    
    Examples
    --------
    >>> xtab = ascii.read('cross_match.txt')
    >>> do_fig_diff(xtab, 'LMC_field1_diff')
    """
    outdir = './Figs_phot%s' % XDIR
    os.makedirs(outdir, exist_ok=True)

    plt.figure(1, (12, 6))
    plt.clf()
    
    plt.subplot(1, 2, 1)
    plt.plot(xtab['G'], xtab['phot_mag']-xtab['G'], '.', alpha=.01)
    plt.plot(xtab['G'], xtab['phot_mag']+xtab['G'], '.', alpha=.01)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam - Gaia (mag)')
    plt.plot([11, 24], [0, 0], 'k-')
    plt.text(13, 2, outroot)
    plt.ylim(-5, 5)
    plt.xlim(11, 22)

    under = xtab[xtab['phot_mag'] > 0]

    plt.subplot(1, 2, 2)
    plt.text(13, 4, 'Under %d Over %d' % (len(under), len(xtab)-len(under)))
    plt.plot(xtab['R'], xtab['phot_mag']-xtab['R'], '.', alpha=.01)
    plt.plot(xtab['R'], xtab['phot_mag']+xtab['R'], '.', alpha=.01)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam - Gaia (mag)')
    plt.plot([11, 24], [0, 0], 'k-')
    plt.ylim(-5, 5)
    plt.xlim(11, 22)
    
    plt.tight_layout()
    plt.savefig('%s/%s.png' % (outdir, outroot))


def get_objects_from_image(filename='LMC_c48_T08.r.t060.fits', outroot=''):
    """
    Detect sources in an image using DAOStarFinder.
    
    Performs source detection on a FITS image and saves the results
    as a FITS table with sky coordinates.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file. Default: 'LMC_c48_T08.r.t060.fits'.
    outroot : str, optional
        Output filename root. If empty, derived from filename. Default: ''.

    Returns
    -------
    str or 'Error'
        Path to output FITS table containing detected sources, or 'Error'
        if file cannot be opened.
    
    Notes
    -----
    **Detection Parameters:**
    
    * FWHM: 4.0 pixels
    * Threshold: 3.0 * background sigma
    * Background: Median-subtracted
    * Sigma estimation: MAD (median absolute deviation)
    
    **Output Table Columns:**
    
    Standard DAOStarFinder columns plus RA and Dec in degrees.
    
    **Output Location:**
    
    ``TabPhot{XDIR}/{outroot}_sources.fits``
    
    Examples
    --------
    >>> sources_file = get_objects_from_image('image.fits')
    >>> sources = Table.read(sources_file)
    >>> print(f"Detected {len(sources)} sources")
    """
    try:
        x = fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

    print('get_photometry: Beginning photometry of %s' % filename)

    xexptime = x['PRIMARY'].header['EXPTIME']
    try:
        xfilter = x['PRIMARY'].header['FILTER']
    except:
        words = filename.split('.')
        xfilter = words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter, filename))

    tab_dir = './TabPhot%s' % XDIR
    os.makedirs(tab_dir, exist_ok=True)

    if outroot == '':
        words = filename.split('/')
        outroot = words[-1].replace('.fits', '')

    # Determine which extension contains the image
    if x[0].data is not None:
        image_wcs = WCS(x[0].header)
        image = x[0].data
    elif x[1].data is not None:
        image_wcs = WCS(x[1].header)
        image = x[1].data

    image -= np.median(image)
    bkg_sigma = mad_std(image)

    daofind = DAOStarFinder(fwhm=4.0, threshold=3.0 * bkg_sigma)
    sources = daofind(image)

    pos = image_wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    sources['RA'] = pos.ra.degree
    sources['Dec'] = pos.dec.degree

    for col in sources.colnames:
        sources[col].info.format = '%.8g'

    outname = '%s/%s_sources.fits' % (tab_dir, outroot)
    sources.write(outname, format='fits', overwrite=True)

    return outname


def locate_first_image_extension(xx):
    """
    Find the first FITS extension containing image data.
    
    Searches through a FITS file to locate the first extension with
    actual image data, handling both primary HDUs and image extensions.

    Parameters
    ----------
    xx : HDUList
        Opened FITS file (from fits.open()).

    Returns
    -------
    int
        Index of the first image extension, or -1 if no image found.
    
    Notes
    -----
    This handles the case where CCD images may have the image in
    extension 1 instead of the primary HDU (extension 0).
    
    Checks for:
    
    * PrimaryHDU
    * ImageHDU
    * CompImageHDU
    
    And verifies that data is not None.
    
    Examples
    --------
    >>> x = fits.open('image.fits')
    >>> ext = locate_first_image_extension(x)
    >>> if ext >= 0:
    ...     image = x[ext].data
    """
    i = 0
    while i < len(xx):
        if isinstance(xx[i], (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)) and xx[i].data is not None:
            return i
        i += 1
    return -1


def do_forced_photometry(filename='LMC_c48_T08.r.t060.fits', object_file='objects.txt',
                         nrows_max=-1, outroot='', rstar=6, b_in=8, b_out=12):
    """
    Perform forced photometry at catalog positions.
    
    Extracts aperture photometry at specified sky positions (typically
    from Gaia catalog) with local background subtraction.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file. Default: 'LMC_c48_T08.r.t060.fits'.
    object_file : str, optional
        Path to file containing source positions (RA, Dec columns required).
        Default: 'objects.txt'.
    nrows_max : int, optional
        Maximum sources to process. If -1, process all. Default: -1.
    outroot : str, optional
        Output filename root. If empty, derived from filename. Default: ''.
    rstar : float, optional
        Aperture radius in pixels. Default: 6.
    b_in : float, optional
        Inner background annulus radius in pixels. Default: 8.
    b_out : float, optional
        Outer background annulus radius in pixels. Default: 12.

    Returns
    -------
    str or 'Error'
        Path to output photometry table, or 'Error' if file cannot be opened.
    
    Notes
    -----
    **NOTE:** This version should be replaced by MefPhot.do_forced_photometry()
    which has been better tested. This version is maintained for compatibility
    but writes output directly within the routine.
    
    **Processing:**
    
    1. Load image and source catalog
    2. Transform sky coordinates to pixel coordinates
    3. Filter sources within detector boundaries
    4. Perform aperture photometry with local background
    5. Calculate magnitudes (zero point = 28)
    6. Write results to TabPhot directory
    
    **Output Table:**
    
    Written to ``TabPhot{XDIR}/{outroot}_phot.txt``
    
    Examples
    --------
    >>> phot_file = do_forced_photometry('image.fits', 'gaia_sources.fits')
    >>> phot = ascii.read(phot_file)
    >>> print(f"Measured {len(phot)} sources")
    """
    try:
        x = fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

    image_ext = locate_first_image_extension(x)
    if image_ext < 0:
        raise IOError('No Image extension in %s' % filename)

    image_wcs = WCS(x[image_ext].header)
    image = x[image_ext].data
    image -= np.median(image)
    image_mask = (image == 0) | ~np.isfinite(image)
    NAXIS1 = x[image_ext].header['NAXIS1']
    NAXIS2 = x[image_ext].header['NAXIS2']

    xexptime = x['PRIMARY'].header['EXPTIME']
    try:
        xfilter = x['PRIMARY'].header['FILTER']
    except:
        words = filename.split('.')
        xfilter = words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter, filename))

    sources = read_table(object_file)
    if 'G' in sources.colnames:
        good = ~sources['R'].mask
        sources = sources[good]

    coords = SkyCoord(ra=sources['RA']*u.deg, dec=sources['Dec']*u.deg)

    # Initialize with NaNs
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

    mask = (
        np.isfinite(sources['xcentroid']) &
        np.isfinite(sources['ycentroid']) &
        (sources['xcentroid'] >= 0) &
        (sources['xcentroid'] < NAXIS1) &
        (sources['ycentroid'] >= 0) &
        (sources['ycentroid'] < NAXIS2)
    )

    sources = sources[mask]
    print(f"Kept {np.sum(mask)} sources on detector out of {len(mask)}")

    npossible = len(sources)

    if nrows_max > 0 and len(sources) > nrows_max:
        sources = random_rows(sources, nrows=nrows_max, seed=None)

    print('Forced photometry of %d of %d possible sources' % (len(sources), npossible))

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=rstar)
    annulus_apertures = CircularAnnulus(positions, r_in=b_in, r_out=b_out)

    phot_table = aperture_photometry(image, apertures)
    aper_stats = ApertureStats(image, apertures, sigma_clip=None, mask=image_mask)
    sigclip = SigmaClip(sigma=3, maxiters=10)
    bkg_stats = ApertureStats(image, annulus_apertures, sigma_clip=sigclip, mask=image_mask)
    total_background = bkg_stats.mean * aper_stats.sum_aper_area.value
    net = aper_stats.sum - total_background

    # Error estimation
    bkg_std_per_pixel = bkg_stats.std
    n_aper_pixels = aper_stats.sum_aper_area.value
    n_bkg_pixels = bkg_stats.sum_aper_area.value

    error = np.sqrt(
        np.abs(net) +
        n_aper_pixels * bkg_std_per_pixel**2 +
        n_aper_pixels**2 * bkg_std_per_pixel**2 / n_bkg_pixels
    )

    phot_table['Raw'] = aper_stats.sum
    phot_table['Bkg'] = total_background
    phot_table['Net'] = net
    phot_table['ErrNet'] = error

    phot_table['phot_mag'] = 28 - 2.5*np.log10(np.fabs(phot_table['Net']))
    phot_table['phot_mag_simple'] = 28 - 2.5*np.log10(np.fabs(phot_table['aperture_sum']))

    phot_table['phot_mag'] = np.select([phot_table['Net'] > 0], [phot_table['phot_mag']],
                                       default=-phot_table['phot_mag'])
    phot_table['phot_mag_simple'] = np.select([phot_table['Net'] > 0], [phot_table['phot_mag_simple']],
                                              default=-phot_table['phot_mag_simple'])

    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'

    pos = image_wcs.pixel_to_world(phot_table['xcenter'], phot_table['ycenter'])
    names = []
    for one in phot_table:
        names.append('x%05d' % one['id'])
    phot_table['Source_name'] = names
    phot_table['RA'] = pos.ra.degree
    phot_table['Dec'] = pos.dec.degree
    phot_table['File'] = outroot
    phot_table['Filter'] = xfilter
    phot_table['Exptime'] = xexptime

    tab_dir = './TabPhot%s' % XDIR
    os.makedirs(tab_dir, exist_ok=True)

    if outroot == '':
        words = filename.split('/')
        outroot = words[-1].replace('.fits', '')

    outfile = '%s/%s_phot.txt' % (tab_dir, outroot)
    phot_table.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)
    print('Wrote %s with %d objects' % (outfile, len(phot_table)))
    return outfile


def find_closest_objects(table1_path, table2_path, max_sep=0.5):
    """
    Cross-match two catalogs finding closest matches.
    
    Finds objects within a specified distance between two astropy tables
    using efficient KDTree algorithm. Returns only the closest match for
    each object in table1.

    Parameters
    ----------
    table1_path : str
        Path to first table (typically Gaia catalog).
    table2_path : str
        Path to second table (typically photometry results).
    max_sep : float, optional
        Maximum separation in arcseconds. Default: 0.5.

    Returns
    -------
    Table or empty list
        Cross-matched table combining columns from both inputs, with added
        'Sep' column giving separation in arcseconds. Returns empty list if
        no matches found.
    
    Notes
    -----
    **Algorithm:**
    
    Uses KDTree in Cartesian coordinates for efficient matching. This is
    much faster than direct spherical distance calculations for large
    catalogs (>1000 sources).
    
    **Optimization (240527):**
    
    Replaced previous implementation with KDTree version, providing
    significant speedup for large catalogs.
    
    **Output:**
    
    Written to ``TabPhot{XDIR}/{table2_name}_x_{table1_name}.txt``
    
    Removes duplicate columns (Source_name, RA, Dec) from table2 to
    avoid conflicts.
    
    Examples
    --------
    >>> xtab = find_closest_objects('gaia.fits', 'phot.txt', max_sep=1.0)
    >>> print(f"Matched {len(xtab)} sources")
    >>> median_sep = np.median(xtab['Sep'])
    >>> print(f"Median separation: {median_sep:.3f} arcsec")
    """
    table1 = read_table(table1_path)
    table2 = read_table(table2_path)

    print('get_closest_objects: Beginning x-match of %s and %s' % (table1_path, table2_path))

    # Convert RA and Dec to SkyCoord objects
    coords1 = SkyCoord(ra=table1['RA'] * u.degree, dec=table1['Dec'] * u.degree)
    coords2 = SkyCoord(ra=table2['RA'] * u.degree, dec=table2['Dec'] * u.degree)

    # Convert to Cartesian coordinates
    cartesian_coords1 = np.array([coords1.cartesian.x.value, coords1.cartesian.y.value, 
                                  coords1.cartesian.z.value]).T
    cartesian_coords2 = np.array([coords2.cartesian.x.value, coords2.cartesian.y.value, 
                                  coords2.cartesian.z.value]).T

    # Build KDTree for efficient matching
    tree = KDTree(cartesian_coords2)

    # Query for closest neighbor
    distances, indices = tree.query(cartesian_coords1)

    # Extract matching rows (resorted to match table1 order)
    closest_matches = table2[indices]

    # Compute separations
    closest_coords = SkyCoord(ra=closest_matches['RA']*u.degree, dec=closest_matches['Dec']*u.degree)
    separations = coords1.separation(closest_coords).arcsecond

    # Create output table
    table1['Sep'] = separations
    del closest_matches['Source_name']
    del closest_matches['RA']
    del closest_matches['Dec']
    table1['Sep'].format = '.3f'
    table1['RA'].format = '.6f'
    table1['Dec'].format = '.6f'
    xtab = hstack([table1, closest_matches])

    xtab = xtab[xtab['Sep'] < max_sep]

    print('Of %d objects in %s and %d objects in %s, found %d matches' % 
          (len(table1), table1_path, len(table2), table2_path, len(xtab)))

    tab_dir = 'TabPhot%s' % XDIR

    if len(xtab):
        words = table1_path.split('/')
        one = words[-1].replace('.txt', '').replace('.fits', '')
        words = table2_path.split('/')
        two = words[-1].replace('.txt', '').replace('.fits', '')
        outfile = '%s/%s_x_%s.txt' % (tab_dir, two, one)
        xtab.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)
    else:
        print('Error: There are no objects that are closer than %f arcsec' % max_sep)
        return []

    return xtab


def get_size(filename='LMC_c48_T08.r.t060.fits'):
    """
    Calculate image center and field size from WCS.
    
    Determines the RA, Dec, and angular size of a FITS image from its
    WCS information.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file. Default: 'LMC_c48_T08.r.t060.fits'.

    Returns
    -------
    ra : float
        Right ascension of image center in degrees.
    dec : float
        Declination of image center in degrees.
    size_deg : float
        Field size in degrees (diagonal from center to corner).
    
    Raises
    ------
    IOError
        If file cannot be opened or WCS information is missing.
    
    Notes
    -----
    Size is calculated as the angular separation from image center to
    corner, providing a conservative estimate of field coverage for
    catalog queries.
    
    Tries extension 0 first, then extension 1 if needed.
    
    Examples
    --------
    >>> ra, dec, size = get_size('image.fits')
    >>> print(f"Field center: RA={ra:.3f}, Dec={dec:.3f}")
    >>> print(f"Field size: {size:.3f} degrees")
    """
    try:
        x = fits.open(filename)
    except:
        print('get_size: Could not open %s' % filename)
        raise IOError('get_size: Could not open %s' % filename)

    try:
        wcs = WCS(x[0].header)
        naxis1 = x[0].header['NAXIS1']
        naxis2 = x[0].header['NAXIS2']
    except:
        try:
            wcs = WCS(x[1].header)
            naxis1 = x[1].header['NAXIS1']
            naxis2 = x[1].header['NAXIS2']
        except:
            raise IOError('get_size: Could not get info for %s' % filename)

    # Calculate center position
    center_pixel = (naxis1 / 2, naxis2 / 2)
    center_ra_dec = wcs.pixel_to_world(center_pixel[0], center_pixel[1])

    # Calculate field size (diagonal)
    corner_pixel = (0, 0)
    corner_ra_dec = wcs.pixel_to_world(corner_pixel[0], corner_pixel[1])
    size_deg = center_ra_dec.separation(corner_ra_dec).to(u.degree).value
    
    ra = center_ra_dec.ra.deg
    dec = center_ra_dec.dec.deg
    
    return ra, dec, size_deg


def do_xphot(filename, gaia_file, forced, nrows_max, outroot):
    """
    Execute photometry and cross-matching pipeline.
    
    Performs photometry (forced or unforced), cross-matches with Gaia,
    and generates diagnostic figures.

    Parameters
    ----------
    filename : str
        Path to FITS image file.
    gaia_file : str
        Path to Gaia catalog file.
    forced : bool
        If True, use forced photometry at Gaia positions. If False,
        detect sources then cross-match.
    nrows_max : int
        Maximum sources to process in forced mode. Ignored for unforced.
    outroot : str
        Output filename root for results.

    Returns
    -------
    None
        Results are written to files and figures are saved.
    
    Notes
    -----
    This is the main pipeline orchestrator that ties together:
    
    1. Photometry (forced or unforced mode)
    2. Cross-matching with Gaia
    3. Figure generation
    
    If cross-matching fails (no matches), prints error and returns without
    generating figures.
    
    Examples
    --------
    >>> do_xphot('image.fits', 'gaia.fits', forced=True, 
    ...          nrows_max=10000, outroot='field1')
    """
    print('XXX - do_xphot %s gaia %s' % (filename, gaia_file))

    if forced:
        object_file = gaia_file
        phot_file = do_forced_photometry(filename, object_file, nrows_max, outroot)
    else:
        object_file = get_objects_from_image(filename, outroot)
        phot_file = do_forced_photometry(filename, object_file, nrows_max=-1, outroot=outroot)

    closest_objects_table = find_closest_objects(gaia_file, phot_file)
    if len(closest_objects_table) == 0:
        print('Error: There are no objects that were xmatched')
        return

    if outroot == '':
        word = filename.split('/')
        outroot = word[-1].replace('.fits', '')

    do_fig(closest_objects_table, outroot)


def do_one(filename='LMC_c48_T08.r.t060.fits', gaia_cat_file='', forced=False, 
           nrows_max=-1, outroot=''):
    """
    Process a single image for photometric comparison.
    
    Complete pipeline for comparing photometry in a single image to Gaia
    catalog, including catalog retrieval, photometry, and figure generation.

    Parameters
    ----------
    filename : str, optional
        Path to FITS file. Default: 'LMC_c48_T08.r.t060.fits'.
    gaia_cat_file : str, optional
        Path to existing Gaia catalog. If empty or non-existent, will
        retrieve new catalog. Default: ''.
    forced : bool, optional
        Photometry mode (True=forced, False=unforced). Default: False.
    nrows_max : int, optional
        Maximum sources for forced photometry. Default: -1 (all).
    outroot : str, optional
        Output filename root. Default: ''.

    Returns
    -------
    None
        Results written to TabPhot/ and figures to Figs_phot/.
    
    Raises
    ------
    ValueError
        If FITS file cannot be opened.
    
    Notes
    -----
    **Gaia Catalog Handling:**
    
    * If gaia_cat_file exists: uses it
    * Otherwise: calculates field center/size and retrieves new catalog
    
    Catalog is cached for reuse in subsequent calls with the same field.
    
    Examples
    --------
    >>> # Use existing Gaia catalog
    >>> do_one('image.fits', gaia_cat_file='gaia.fits', forced=True)
    
    >>> # Auto-retrieve Gaia catalog
    >>> do_one('image.fits', forced=True, nrows_max=5000)
    """
    try:
        x = fits.open(filename)
    except:
        print('Could not open %s' % filename)
        raise ValueError

    if gaia_cat_file != '' and os.path.isfile(gaia_cat_file) == True:
        gaia_file = gaia_cat_file
        print('Using existing GaiaCat file: %s' % gaia_cat_file)
    else:
        ra, dec, size_deg = get_size(filename)
        print('Making new GaiaCat file - %.2f %.2f %.2f' % (ra, dec, size_deg))
        print('do_one - RA, Dec, size: ', ra, dec, size_deg)
        gaia_file = get_gaia(ra, dec, size_deg, outroot)

    do_xphot(filename, gaia_file, forced, nrows_max, outroot)
    return


def do_many(filenames=['LMC_c48_T08.r.t060.fits'], gaia_cat_file='', forced=True, 
            nrows_max=10000, outroot=''):
    """
    Process multiple images with optimized Gaia catalog retrieval.
    
    Efficiently processes multiple images by identifying unique field
    positions and reusing Gaia catalogs for overlapping fields.

    Parameters
    ----------
    filenames : list of str, optional
        List of FITS files to process. Default: ['LMC_c48_T08.r.t060.fits'].
    gaia_cat_file : str, optional
        Ignored (kept for API compatibility). Default: ''.
    forced : bool, optional
        Photometry mode. Default: True.
    nrows_max : int, optional
        Maximum sources for forced photometry. Default: 10000.
    outroot : str, optional
        Output filename root. Default: ''.

    Returns
    -------
    None
        Results written to files.
    
    Raises
    ------
    IOError
        If any file cannot be opened.
    
    Notes
    -----
    **Optimization Strategy:**
    
    1. Calculate field centers and sizes for all files
    2. Identify unique fields (within 0.01° tolerance)
    3. Retrieve Gaia catalogs only for unique fields
    4. Map each file to its Gaia catalog
    5. Process all files using cached catalogs
    
    This dramatically reduces Gaia query time when processing many images
    of the same field (e.g., different filters or epochs).
    
    **Intermediate Files:**
    
    * xpos.txt - All file positions
    * zpos.txt - Unique field positions
    * xxpos.txt - Files with assigned Gaia catalogs
    
    Examples
    --------
    >>> files = ['field1_r.fits', 'field1_g.fits', 'field1_i.fits']
    >>> do_many(files, forced=True, nrows_max=5000)
    Finished getting gaia tables for 1 files
    Processing images...
    """
    xra = []
    xdec = []
    xsize = []

    for filename in filenames:
        try:
            x = fits.open(filename)
        except:
            print('do_many: Could not open %s' % filename)
            raise IOError
        ra, dec, size = get_size(filename)
        xra.append(ra)
        xdec.append(dec)
        xsize.append(size)

    xpos = Table([filenames, xra, xdec, xsize], names=['filename', 'RA', 'Dec', 'Size'])
    zpos, mapping = unique_rows_within_tol(xpos, tol=0.01)

    print("Finished getting positions")

    gaia_files = []
    for one in zpos:
        gaia_file = get_gaia(one['RA'], one['Dec'], one['Size'], outroot='')
        gaia_files.append(gaia_file)
    zpos['gaia_file'] = gaia_files

    print('Finished getting gaia tables for %d files' % len(zpos))

    xpos.write('xpos.txt', format='ascii.fixed_width_two_line', overwrite=True)
    zpos.write('zpos.txt', format='ascii.fixed_width_two_line', overwrite=True)

    xpos['gaia_file'] = zpos['gaia_file'][mapping]
    xpos.write('xxpos.txt', format='ascii.fixed_width_two_line', overwrite=True)

    # Process all files using cached Gaia catalogs
    for one in xpos:
        print('ZZZ', one)
        do_xphot(one['filename'], one['gaia_file'], forced, nrows_max, outroot)

    return


def do_dir(xdir='DECam_SWARP2/LMC_c37/T16', nrows_max=30000, forced=True):
    """
    Process all images in a directory and subdirectories.
    
    Recursively finds all FITS files in a directory tree and processes
    them with optimized Gaia catalog retrieval.

    Parameters
    ----------
    xdir : str, optional
        Directory path to process. Default: 'DECam_SWARP2/LMC_c37/T16'.
    nrows_max : int, optional
        Maximum sources for forced photometry. Default: 30000.
    forced : bool, optional
        Photometry mode. Default: True.

    Returns
    -------
    None
        Results written to files.
    
    Notes
    -----
    Uses ImageSum.table_create() to recursively find all FITS files.
    Then calls do_many() to process with optimized Gaia catalog caching.
    
    This is the recommended approach for processing large datasets where
    multiple images cover the same fields.
    
    Examples
    --------
    >>> do_dir('DECamSWARP2/SMC_c01', nrows_max=20000, forced=True)
    Starting 145 files
    Finished getting gaia tables for 12 files
    Processing images...
    """