#!/usr/bin/env python
# coding: utf-8

"""Smash - SMASH DR2 Catalog Retrieval

Space Telescope Science Institute

Synopsis
--------

Retrieve stars from the SMASH (Survey of the MAgellanic Stellar History) DR2
catalog for a specified region on the sky.

Command Line Usage
------------------

::

    Smash.py [-h] [-rad size_deg] [-rmag max] [-kfrac frac] [-plot] [-ex] [-out root] ra dec

**Required Arguments:**

ra
    Right Ascension of cone center in degrees.

dec
    Declination of cone center in degrees.

**Optional Arguments:**

-h
    Print this help message and exit.

-rad size
    Set the cone search radius in degrees. Default: 0.5.

-rmag max
    Set the maximum r-band magnitude. Default: 22.0.

-kfrac frac
    Set the fraction of objects to keep. Default: 0.5.

-plot
    Create a scatter plot of selected positions.

-ex
    Run example query (SMC region) without further inputs.

-out root
    Set the output filename root. Default: smash_cat.

Description
-----------

This module retrieves photometric data from the SMASH DR2 catalog using
the NOAO Data Lab query client. It performs a cone search around a specified
position and returns filtered sources with high probability stellar
classifications.

The output table includes:

- Position (RA, Dec)
- Photometry in U, G, R, Z bands
- Probability of being a star (prob)

Sources are filtered using select_fraction_smash_dr2, which applies a
PSF-likeness score based on prob, chi, and sharp parameters to select
a consistent fraction of high-quality stars across different sky regions.

Query results are cached in the local 'Smash' directory to avoid redundant
queries to the SMASH archive. Cached files include header keywords describing
the query parameters for verification.

Primary Routines
----------------

get_smash
    Retrieve a SMASH catalog for a field and return the file path.
    Tries a pre-assembled local file first, falls back to the archive.
    This is the main entry point for use by MefPhot and other pipeline
    tools; it mirrors the interface of GaiaCat.get_gaia.

get_smash_from_file
    Extract sources from a pre-assembled local SMASH catalog file
    (default: Smash_MagClouds.fits in the current directory or
    $KRED/xdata/), mirroring GaiaCat.get_gaia_from_file.

smash_cone_search
    Perform a cone search on SMASH DR2 catalog (with caching).

select_fraction_smash_dr2
    Select a fixed fraction of objects using PSF-likeness scoring.

plot_positions
    Create a scatter plot of selected objects vs RA and Dec with histograms.

do_one
    Retrieve and filter SMASH catalog for a single sky position and write
    the result to a FITS file.

steer
    Command line interface for catalog retrieval.

Caching Functions
-----------------

get_cache_filename
    Generate cache filename from query parameters.

check_cache
    Check if cached query result exists.

save_to_cache
    Save query result to cache with header keywords.

Output Column Names
-------------------

After filtering, column names are standardised to match the Gaia
convention used throughout the pipeline:

* RA, Dec  -- sky coordinates (degrees)
* U, G, R, Z -- DECam ugriz photometry (magnitudes)

The G and R columns are in the DECam photometric system (not Gaia's
broadband G and R), so the color term c_1 in ZeroCalc is expected to
be small for DECam r-band data.

Notes
-----

* Requires the NOAO Data Lab client library (``dl``).
* SMASH DR2 covers the Magellanic Cloud region only; use GaiaCat for
  fields outside this footprint.
* Raw query results are cached in ``Smash/`` to avoid redundant archive
  queries.  The filtered output used by MefPhot is also cached in
  ``Smash/`` with a filename that encodes all query parameters.

Version History
---------------

260108 ksl
    Coding begun

260414 ksl
    Added get_smash() as the pipeline integration entry point.

"""

import sys
from astropy.io import ascii, fits
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import os

from dl import queryClient as qc

# Cache directory for raw SMASH query results
SMASH_CACHE_DIR = 'Smash'

# Module-level cache for the pre-assembled catalog (populated on first use)
_smash_preassembled_table = None
_smash_preassembled_path  = None


def get_cache_filename(ra, dec, radius):
    """
    Generate a cache filename based on query parameters.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in degrees.

    Returns
    -------
    str
        Path to the cache file.
    """
    # Round to avoid floating point issues in filenames
    filename = f"smash_ra{ra:.4f}_dec{dec:.4f}_r{radius:.3f}.fits"
    return os.path.join(SMASH_CACHE_DIR, filename)


def check_cache(ra, dec, radius, tol=1e-4):
    """
    Check if a cached query result exists for the given parameters.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in degrees.
    tol : float, optional
        Tolerance for matching query parameters (default: 1e-4).

    Returns
    -------
    astropy.table.Table or None
        Cached table if found and parameters match, None otherwise.
    """
    cache_file = get_cache_filename(ra, dec, radius)

    if not os.path.exists(cache_file):
        return None

    try:
        # Read the cached file and check header keywords
        with fits.open(cache_file) as hdu:
            header = hdu[1].header
            cached_ra = header.get('QUERY_RA', None)
            cached_dec = header.get('QUERY_DE', None)
            cached_rad = header.get('QUERY_R', None)

            if cached_ra is None or cached_dec is None or cached_rad is None:
                print(f"Cache file {cache_file} missing query keywords, will re-query")
                return None

            # Check if parameters match within tolerance
            if (abs(cached_ra - ra) < tol and
                abs(cached_dec - dec) < tol and
                abs(cached_rad - radius) < tol):
                print(f"Using cached query result from {cache_file}")
                return Table.read(cache_file)
            else:
                print(f"Cache file parameters don't match, will re-query")
                return None
    except Exception as e:
        print(f"Error reading cache file {cache_file}: {e}")
        return None


def save_to_cache(table, ra, dec, radius):
    """
    Save a query result to the cache with header keywords.

    Parameters
    ----------
    table : astropy.table.Table
        The query result table.
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in degrees.

    Returns
    -------
    str
        Path to the saved cache file.
    """
    # Create cache directory if needed
    if not os.path.exists(SMASH_CACHE_DIR):
        os.makedirs(SMASH_CACHE_DIR)
        print(f"Created cache directory: {SMASH_CACHE_DIR}")

    cache_file = get_cache_filename(ra, dec, radius)

    # Write table to FITS
    table.write(cache_file, format='fits', overwrite=True)

    # Add header keywords describing the query
    # Coerce to plain Python floats in case numpy scalars were passed
    with fits.open(cache_file, mode='update') as hdu:
        hdu[1].header['QUERY_RA'] = (float(ra), 'Query center RA (deg)')
        hdu[1].header['QUERY_DE'] = (float(dec), 'Query center Dec (deg)')
        hdu[1].header['QUERY_R'] = (float(radius), 'Query radius (deg)')
        hdu[1].header['CATALOG'] = ('SMASH_DR2', 'Source catalog')
        hdu[1].header['NOBJ_RAW'] = (len(table), 'Number of objects in raw query')
        hdu.flush()

    print(f"Saved raw query to cache: {cache_file}")
    return cache_file


def smash_cone_search(ra, dec, radius_deg, limit=10000000, verbose=True, use_cache=True):
    """
    Perform a cone search on SMASH DR2 using Data Lab query client.

    Queries the SMASH DR2 object catalog for all sources within a
    specified angular radius of the given sky position. Results are
    cached locally to avoid redundant queries.

    Parameters
    ----------
    ra : float
        Right Ascension of the cone center in degrees.
    dec : float
        Declination of the cone center in degrees.
    radius_deg : float
        Radius of the cone search in degrees.
    limit : int, optional
        Maximum number of rows to return (default: 10000000).
    verbose : bool, optional
        If True, write raw query to foo.fits (default: True).
    use_cache : bool, optional
        If True, check for and use cached results (default: True).

    Returns
    -------
    astropy.table.Table
        Table of SMASH DR2 objects within the search cone.
        Columns include position, photometry, and classification info.

    Examples
    --------
    >>> table = smash_cone_search(13.1867, -72.8286, 0.5)
    >>> print(f"Found {len(table)} sources")
    """
    # Check cache first
    if use_cache:
        cached = check_cache(ra, dec, radius_deg)
        if cached is not None:
            print(f"Retrieved {len(cached)} objects from cache")
            return cached

    sql = f"""
    SELECT *
    FROM smash_dr2.object
    WHERE q3c_radial_query(
        ra, dec,
        {ra}, {dec}, {radius_deg}
    )
    LIMIT {limit}
    """
    # Query returns an Astropy Table directly
    xtable = qc.query(sql=sql, fmt="table")
    print("The cone search returned %d objects compared to the limit of %d or  %.2f percent of the max" % (len(xtable),limit,100.*len(xtable)/limit))
    if verbose:
        xtable.write('foo.fits',format='fits',overwrite=True)

    # Save to cache
    if use_cache:
        save_to_cache(xtable, ra, dec, radius_deg)

    return xtable



def select_fraction_smash_dr2(
    tab,
    rmag_max,
    keep_frac,
    rmag_ref=21.0,
    sigma_chi=0.5,
    sharp0=0.05,
    eps=1e-3
):
    """
    Select a fixed fraction of a SMASH DR2 object table for artificial star subtraction.

    Parameters
    ----------
    tab : astropy.table.Table
        SMASH DR2 object table.
    rmag_max : float
        Hard faint-end magnitude limit.
    keep_frac : float
        Fraction (0 < keep_frac <= 1) of objects to keep after magnitude cut.
    rmag_ref : float
        Reference magnitude for sharpness scaling.
    sigma_chi : float
        Width of acceptable chi distribution.
    sharp0 : float
        Base sharpness scatter at bright magnitudes.
    eps : float
        Small number to stabilize log(prob).

    Returns
    -------
    astropy.table.Table
        Subset of input table containing exactly keep_frac of eligible objects.
    """

    if not (0 < keep_frac <= 1):
        raise ValueError("keep_frac must be in (0, 1].")

    # ------------------------------------------------------------
    # 1. Hard magnitude cut + sanity filters
    # ------------------------------------------------------------
    good = (
        np.isfinite(tab['rmag']) &
        np.isfinite(tab['prob']) &
        np.isfinite(tab['chi']) &
        np.isfinite(tab['sharp']) &
        (tab['rmag'] <= rmag_max) &
        (tab['ndetr'] > 0)
    )

    t = tab[good].copy()
    n_total = len(t)
    if n_total == 0:
        return t

    # ------------------------------------------------------------
    # 2. PSF-likeness score
    # ------------------------------------------------------------
    dm = np.clip(t['rmag'] - rmag_ref, 0, None)

    sigma_sharp = np.sqrt(sharp0**2 + (0.02 * dm)**2)

    score = (
        np.log(t['prob'] + eps)
        - 0.5 * ((t['chi'] - 1.0) / sigma_chi)**2
        - 0.5 * (t['sharp'] / sigma_sharp)**2
    )

    # ------------------------------------------------------------
    # 3. Quantile selection
    # ------------------------------------------------------------
    n_keep = int(np.floor(keep_frac * n_total))
    if n_keep < 1:
        return t[:0]

    thresh = np.partition(score, -n_keep)[-n_keep]
    keep = score >= thresh

    return t[keep]


def plot_positions(table, outfile='smash_positions.png', title=None,
                   rmag_max=None, keep_frac=None, n_original=None):
    """
    Create a scatter plot of selected objects as a function of RA and Dec.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing RA and Dec columns (accepts 'RA'/'Dec' or 'ra'/'dec').
    outfile : str, optional
        Output filename for the plot (default: 'smash_positions.png').
    title : str, optional
        Title for the plot. If None, shows number of objects.
    rmag_max : float, optional
        Maximum r-band magnitude used for filtering (for annotation).
    keep_frac : float, optional
        Fraction of objects kept after filtering (for annotation).
    n_original : int, optional
        Number of objects in original query before filtering (for annotation).

    Returns
    -------
    None
        Saves plot to outfile.

    Notes
    -----
    Uses small markers ('.') and low alpha to visualize density variations
    when plotting large numbers of objects. Includes marginal histograms
    showing the distribution of sources in RA and Dec. The RA axis is
    scaled by cos(dec) to account for spherical projection.
    """
    # Handle both column naming conventions
    if 'RA' in table.colnames:
        ra_col, dec_col = 'RA', 'Dec'
    else:
        ra_col, dec_col = 'ra', 'dec'

    ra = np.array(table[ra_col])
    dec = np.array(table[dec_col])

    # Adjust alpha based on number of objects
    n_obj = len(table)
    if n_obj > 100000:
        alpha = 0.05
    elif n_obj > 50000:
        alpha = 0.1
    elif n_obj > 10000:
        alpha = 0.2
    elif n_obj > 1000:
        alpha = 0.3
    else:
        alpha = 0.5

    # Create figure with marginal histograms using gridspec
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                          hspace=0.05, wspace=0.05)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)

    # Main scatter plot
    ax_main.scatter(ra, dec, marker='.', s=1, alpha=alpha, c='blue')
    ax_main.set_xlabel('RA (deg)')
    ax_main.set_ylabel('Dec (deg)')
    ax_main.invert_xaxis()  # RA increases to the left

    # Apply cos(dec) correction for aspect ratio
    dec_mean = np.mean(dec)
    cos_dec = np.cos(np.radians(dec_mean))
    ax_main.set_aspect(1.0 / cos_dec, adjustable='box')

    # RA histogram (bottom, but placed at top in this layout)
    ax_histx.hist(ra, bins=50, color='blue', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histx.set_ylabel('N')
    ax_histx.invert_xaxis()  # Match main plot RA direction

    # Dec histogram (right side)
    ax_histy.hist(dec, bins=50, orientation='horizontal', color='blue', alpha=0.7,
                  edgecolor='black', linewidth=0.5)
    ax_histy.tick_params(axis='y', labelleft=False)
    ax_histy.set_xlabel('N')

    if title is None:
        title = f'SMASH DR2 Selected Objects (N={n_obj:,})'
    ax_histx.set_title(title)

    # Add selection criteria info box
    info_lines = []
    if n_original is not None:
        actual_frac = n_obj / n_original
        info_lines.append(f'Original: {n_original:,}')
        info_lines.append(f'Selected: {n_obj:,} ({100*actual_frac:.1f}%)')
    if rmag_max is not None:
        info_lines.append(f'rmag_max: {rmag_max}')
    if keep_frac is not None:
        info_lines.append(f'keep_frac: {keep_frac}')

    if info_lines:
        info_text = '\n'.join(info_lines)
        ax_main.text(0.02, 0.98, info_text, transform=ax_main.transAxes,
                     fontsize=9, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved position plot to {outfile}")


def do_one(ra, dec, radius=0.5, outroot='smash_cat', rmag_max=22.0, keep_frac=0.5, plot=False):
    """
    Retrieve and filter SMASH DR2 catalog for a single sky position.

    Performs a cone search on the SMASH DR2 catalog, applies quality
    filtering using select_fraction_smash_dr2, standardizes column names,
    and writes the result to a FITS file.

    Parameters
    ----------
    ra : float
        Right Ascension of the cone center in degrees.
    dec : float
        Declination of the cone center in degrees.
    radius : float, optional
        Radius of the cone search in degrees (default: 0.5).
    outroot : str, optional
        Output filename root (default: 'smash_cat').
    rmag_max : float, optional
        Maximum r-band magnitude for filtering (default: 22.0).
    keep_frac : float, optional
        Fraction of objects to keep after quality filtering (default: 0.5).
    plot : bool, optional
        If True, create a scatter plot of selected positions (default: False).

    Returns
    -------
    astropy.table.Table
        The filtered and column-renamed table that was written to file.

    Notes
    -----
    Uses select_fraction_smash_dr2 to select a consistent fraction of
    high-quality stellar sources across different sky regions.

    Column names are standardized: ra->RA, dec->Dec, umag->U, etc.
    """
    print("Querying SMASH DR2...")
    table = smash_cone_search(ra, dec, radius)
    n_original = len(table)
    print(f"Retrieved {n_original} rows")

    # Use select_fraction_smash_dr2 for consistent filtering across sky regions
    table = select_fraction_smash_dr2(table, rmag_max=rmag_max, keep_frac=keep_frac)
    print(f"After select_fraction_smash_dr2 (rmag_max={rmag_max}, keep_frac={keep_frac}): {len(table)} rows")

    table.pprint(max_lines=10)

    # Standardize column names
    table.rename_column('ra', 'RA')
    table.rename_column('dec', 'Dec')
    table.rename_column('umag', 'U')
    table.rename_column('gmag', 'G')
    table.rename_column('rmag', 'R')
    table.rename_column('zmag', 'Z')

    if outroot == '':
        outroot = 'Smash_%06.2f_%06.2f_%5.2f' % (ra, dec, radius)
    outfile = outroot + '.fits'
    table.write(outfile, format='fits', overwrite=True)

    # Add header keywords describing the query and filtering
    with fits.open(outfile, mode='update') as hdu:
        hdu[1].header['QUERY_RA'] = (float(ra), 'Query center RA (deg)')
        hdu[1].header['QUERY_DE'] = (float(dec), 'Query center Dec (deg)')
        hdu[1].header['QUERY_R'] = (float(radius), 'Query radius (deg)')
        hdu[1].header['CATALOG'] = ('SMASH_DR2', 'Source catalog')
        hdu[1].header['RMAG_MAX'] = (rmag_max, 'Max r-band magnitude filter')
        hdu[1].header['KEEPFRAC'] = (keep_frac, 'Fraction of objects kept')
        hdu[1].header['NOBJ_RAW'] = (n_original, 'Number of objects in raw query')
        hdu[1].header['NOBJ_SEL'] = (len(table), 'Number of selected objects')
        hdu.flush()

    print(f"Wrote {len(table)} sources to {outfile}")

    if plot:
        plot_positions(table, outfile=outroot + '_positions.png',
                       rmag_max=rmag_max, keep_frac=keep_frac, n_original=n_original)

    return table


def get_smash_from_file(ra, dec, size_deg, filename='Smash_MagClouds.fits',
                        outroot='', rmag_max=22.0, keep_frac=0.5):
    """
    Extract SMASH DR2 sources from a local pre-assembled catalog file.

    Performs a rectangular RA/Dec selection from a pre-downloaded SMASH
    catalog, applies quality filtering, standardises column names, and
    writes the subset to the ``Smash/`` cache directory.  This mirrors
    the behaviour of ``GaiaCat.get_gaia_from_file()``.

    Parameters
    ----------
    ra : float
        Right Ascension of the field centre in degrees.
    dec : float
        Declination of the field centre in degrees.
    size_deg : float
        Half-width of the extraction region in degrees.  The RA range is
        widened by ``1/cos(dec)`` to preserve equal angular coverage in
        both dimensions.
    filename : str, optional
        Name of the pre-assembled SMASH catalog.  Searched first in the
        current directory, then in ``$KRED/xdata/``.
        Default: ``'Smash_MagClouds.fits'``.
    outroot : str, optional
        Root name for the output file.  If empty, constructed from RA and
        Dec.
    rmag_max : float, optional
        Faint-end r-band magnitude limit passed to
        ``select_fraction_smash_dr2`` (default: 22.0).
    keep_frac : float, optional
        Fraction of quality-selected stars to retain (default: 0.5).

    Returns
    -------
    str
        Path to the output FITS file containing the filtered sources with
        standardised column names (RA, Dec, U, G, R, Z).

    Raises
    ------
    IOError
        If ``filename`` cannot be located in the current directory or in
        ``$KRED/xdata/``.
    """
    # Locate the pre-assembled file
    xfilename = ''
    if os.path.isfile(filename):
        xfilename = filename
    else:
        KRED = os.environ.get('KRED')
        if KRED is not None:
            candidate = os.path.join(KRED, 'xdata', filename)
            if os.path.isfile(candidate):
                xfilename = candidate
            else:
                raise IOError(f'Could not locate {filename} locally or in $KRED/xdata/')
        else:
            raise IOError(
                f'Could not locate {filename} locally and KRED environment variable is not set'
            )

    global _smash_preassembled_table, _smash_preassembled_path
    if _smash_preassembled_path == xfilename and _smash_preassembled_table is not None:
        xtab = _smash_preassembled_table
    else:
        size_gb = os.path.getsize(xfilename) / 1e9
        print(f'get_smash_from_file: loading {xfilename} ({size_gb:.1f} GB)...')
        xtab = Table.read(xfilename)
        _smash_preassembled_table = xtab
        _smash_preassembled_path  = xfilename
        print(f'get_smash_from_file: loaded {len(xtab):,} rows into memory')

    # Detect column naming convention:
    #   raw archive  -> lowercase 'ra', 'dec', 'rmag', 'umag', 'gmag', 'zmag'
    #   pre-assembled -> standardised 'RA', 'Dec', 'R', 'U', 'G', 'Z'
    raw_names = 'ra' in xtab.colnames
    ra_col  = 'ra'  if raw_names else 'RA'
    dec_col = 'dec' if raw_names else 'Dec'

    # Rectangular sky selection (RA range scaled by cos(dec))
    dec_min = dec - 0.5 * size_deg
    dec_max = dec + 0.5 * size_deg
    xscale = np.cos(np.radians(dec))
    factor = 0.5 * size_deg / xscale
    ra_min = ra - factor
    ra_max = ra + factor

    mask = (
        (xtab[ra_col]  > ra_min) & (xtab[ra_col]  < ra_max) &
        (xtab[dec_col] > dec_min) & (xtab[dec_col] < dec_max)
    )
    ftab = xtab[mask]

    # select_fraction_smash_dr2 expects the raw column name 'rmag'.
    # Pre-assembled files already call it 'R', so rename temporarily.
    if not raw_names:
        ftab.rename_column('R', 'rmag')

    ftab = select_fraction_smash_dr2(ftab, rmag_max=rmag_max, keep_frac=keep_frac)

    # Standardise column names to match pipeline convention
    if raw_names:
        ftab.rename_column('ra',   'RA')
        ftab.rename_column('dec',  'Dec')
        ftab.rename_column('umag', 'U')
        ftab.rename_column('gmag', 'G')
        ftab.rename_column('rmag', 'R')
        ftab.rename_column('zmag', 'Z')
    else:
        # Undo the temporary rename used for filtering
        ftab.rename_column('rmag', 'R')

    os.makedirs(SMASH_CACHE_DIR, exist_ok=True)
    if outroot == '':
        outroot = f'{ra:.4f}_{dec:+.4f}'
    outfile = os.path.join(SMASH_CACHE_DIR, f'Smash.{outroot}.fits')

    ftab.write(outfile, format='fits', overwrite=True)
    return outfile


def get_smash(ra, dec, size, rmag_max=22.0, keep_frac=0.5):
    """
    Retrieve a SMASH DR2 catalog for the given field, using a local cache
    if available.  Returns the path to the catalog FITS file.

    This is the SMASH analogue of GaiaCat.get_gaia(ra, dec, size) and
    provides the same calling interface: given a field centre and a
    half-diagonal radius in degrees it returns a file path that can be
    passed directly to MefPhot.do_forced_photometry as ``object_file``.

    The returned table has standardised column names RA, Dec, U, G, R, Z
    and is compatible with ZeroCalc without further processing.

    The function first attempts to extract sources from a pre-assembled
    local file (``Smash_MagClouds.fits``, searched in the current directory
    then in ``$KRED/xdata/``).  If the file is not available it falls back
    to a live cone search against the NOAO Data Lab archive.

    Parameters
    ----------
    ra : float
        Right Ascension of the field centre in degrees.
    dec : float
        Declination of the field centre in degrees.
    size : float
        Half-diagonal search radius in degrees (same convention as
        get_gaia).
    rmag_max : float, optional
        Faint-end r-band magnitude limit passed to do_one (default: 22.0).
    keep_frac : float, optional
        Fraction of quality-selected stars to retain (default: 0.5).

    Returns
    -------
    str
        Path to a FITS file containing the filtered SMASH DR2 catalog.

    Notes
    -----
    When using the pre-assembled file the output is cached in
    ``Smash/Smash.<ra>_<dec>.fits``.  When falling back to the archive
    the output is cached using a filename that also encodes rmag_max and
    keep_frac so that queries with different parameters produce separate
    cache entries.
    """
    # Try the pre-assembled local file first
    try:
        outroot = f'{ra:.4f}_{dec:+.4f}_r{rmag_max:.1f}_k{keep_frac:.2f}'
        outfile = os.path.join(SMASH_CACHE_DIR, f'Smash.{outroot}.fits')
        if os.path.exists(outfile):
            return outfile
        return get_smash_from_file(ra, dec, size,
                                   rmag_max=rmag_max, keep_frac=keep_frac,
                                   outroot=outroot)
    except IOError:
        print('get_smash: pre-assembled file not available, falling back to archive')

    # Fall back to live archive query
    os.makedirs(SMASH_CACHE_DIR, exist_ok=True)
    outroot = os.path.join(
        SMASH_CACHE_DIR,
        f'smash_{ra:.4f}_{dec:+.4f}_{size:.3f}_r{rmag_max:.1f}_k{keep_frac:.2f}'
    )
    outfile = outroot + '.fits'
    if os.path.exists(outfile):
        return outfile
    do_one(ra, dec, radius=size, outroot=outroot,
           rmag_max=rmag_max, keep_frac=keep_frac, plot=False)
    return outfile


def steer(argv):
    """
    Command line interface for SMASH catalog retrieval.

    Parses command line arguments and calls do_one to perform the
    cone search and filtering.

    Parameters
    ----------
    argv : list
        Command line arguments (typically sys.argv).
    """
    ra = -1.
    dec = -1.
    radius = 0.5  # degrees
    outroot = 'smash_cat'
    rmag_max = 22.0
    keep_frac = 0.5
    plot = False

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i][:4] == '-rad':
            i += 1
            radius = eval(argv[i])
        elif argv[i][:4] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i][:5] == '-rmag':
            i += 1
            rmag_max = eval(argv[i])
        elif argv[i][:6] == '-kfrac':
            i += 1
            keep_frac = eval(argv[i])
        elif argv[i][:5] == '-plot':
            plot = True
        elif argv[i][:3] == '-ex':
            ra = 13.1867
            dec = -72.8286
            break
        elif argv[i][0] == '-' and ra == -1:
            print('Error: Unknown option :', argv)
            return
        elif ra == -1.:
            ra = eval(argv[i])
        elif dec == -1.:
            dec = eval(argv[i])
        else:
            print('Error: Cannot parse command line: ', argv)
            return
        i += 1

    if ra == -1 or dec == -1:
        print('Error: RA and Dec are required')
        print(__doc__)
        return

    print("RA:", ra)
    print("Dec:", dec)
    print("Radius (deg):", radius)
    print("R_mag:", rmag_max)
    print("Fraction:", keep_frac)
    print("Out root:", outroot)
    print("Plot:", plot)


    do_one(ra, dec, radius=radius, outroot=outroot, rmag_max=rmag_max, keep_frac=keep_frac, plot=plot)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)

