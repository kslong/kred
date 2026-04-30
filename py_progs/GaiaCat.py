#!/usr/bin/env python
# coding: utf-8
"""GaiaCat - GAIA Database Interaction Module.

Space Telescope Science Institute

Synopsis
--------

Retrieve GAIA catalog data either from the archive or from a local file.
Can determine the search center from a FITS image WCS or from explicit
RA/Dec coordinates.

Command Line Usage
------------------

::

    GaiaCat.py [-h] [-archive] [-redo] [-gfile FILENAME] [-rad DEGREES] [-out OUTROOT]
               input.fits or RA Dec

Arguments
---------

input
    Either a FITS file (to extract center from WCS), or RA in degrees.
    If RA is given, dec must also be provided.

RA, Dec  RA and DEC of field center

or

whatever.fits  a fits file with a WCS,  note that size is not taken from WCS

Options
-------

-h
    Print this help message and exit

-archive
    Retrieve from GAIA archive instead of local file. By default,
    tries local file first, then falls back to archive.

-redo
    Force re-retrieval from the archive even if the output file already
    exists. Implies -archive.

-gfile FILENAME
    Name of local GAIA catalog file. Searches locally first, then
    in $KRED/xdata/. Default: Gaia_MagClouds.fits

-rad DEGREES
    Search radius in degrees. Default: 0.5

-out OUTROOT
    Output filename root. Default: derived from RA/Dec or FITS filename

Description
-----------

This module provides functions to retrieve GAIA DR3 photometric data.
The primary functions are:

- get_gaia_from_archive(): Query the GAIA archive with cone search
- get_gaia_from_file(): Extract from a local pre-downloaded catalog

Notes
-----

This module requires the astroquery package for archive access.
All archive-dependent functions will provide clear error messages if
astroquery is not available or if GAIA services are unreachable.

History:

240318 ksl Coding begun
240527 ksl Speed up the catalog matching
251105 ksl Split finding sources in an image from doing photometry
251130 ksl Cleaned up to focus on GAIA catalog interaction
251130 ksl Robust handling of astroquery import vs service availability
251211 ksl Handle gaiaxpy version compatibility (2.1.1 vs 2.1.2)
250116 ksl Added command-line steering with FITS/WCS support

Example Usage
-------------

Command line with FITS file::

    $ GaiaCat.py myimage.fits -rad 0.3

Command line with RA/Dec::

    $ GaiaCat.py 84.925 -66.274 -rad 0.5

Force archive retrieval (skip local file)::

    $ GaiaCat.py 84.925 -66.274 -archive -rad 0.3

Force re-retrieval from archive even if output file exists::

    $ GaiaCat.py 84.925 -66.274 -redo -rad 0.3

Use a different local catalog file::

    $ GaiaCat.py 84.925 -66.274 -gfile my_gaia_catalog.fits

Python usage::

    >>> from GaiaCat import get_gaia_from_archive
    >>> outfile = get_gaia_from_archive(ra=84.925, dec=-66.274, rad_deg=0.3)
    >>> print(f"Catalog saved to: {outfile}")
"""


import os
import time
import timeit
import pathlib
import os.path as path

import numpy as np
from http.client import IncompleteRead

from astropy.io import fits, ascii
from astropy.table import Table, join, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs


# --------------------------------------------------------------------------------
# Gaia loader (astroquery) — distinguishes dependency vs asset/service problems
# --------------------------------------------------------------------------------
def load_Gaia(probe_service=True, credentials_file=None):
    """Return the astroquery.gaia.Gaia class with service validation.

    This function provides a controlled import of the GAIA query interface,
    distinguishing between installation issues and service availability problems.

    Parameters
    ----------
    probe_service : bool, optional
        If True (default), perform a minimal network check to verify that
        external GAIA services are reachable. If False, only import the
        Gaia class without network validation.
    credentials_file : str, optional
        Path to file containing Gaia credentials (username on line 1,
        password on line 2). If None, defaults to ~/.gaia_credentials
        if it exists.

    Returns
    -------
    Gaia : class
        The Gaia class from astroquery.gaia, ready for use in queries.

    Raises
    ------
    RuntimeError
        If astroquery is not installed, or if ``probe_service=True`` and
        the GAIA external services are unavailable or unreachable.

    Notes
    -----
    This function also configures astropy logging to reduce verbosity during
    GAIA operations by setting the log level to ERROR and disabling IERS
    auto-downloads.

    If a credentials file exists, the function will automatically login
    to the Gaia archive.

    Examples
    --------
    Load GAIA with service validation::

        >>> Gaia = load_Gaia(probe_service=True)
        >>> # Now safe to use Gaia for queries

    Load GAIA without immediate service check::

        >>> Gaia = load_Gaia(probe_service=False)
        >>> # Service errors will occur at first query attempt

    """
    import os

    # 1) Import-time: distinguish "not installed"
    try:
        from astroquery.gaia import Gaia
    except ImportError as e:
        raise RuntimeError(
            "astroquery (or 'astroquery.gaia') is not installed in this environment."
        ) from e

    # Optional: reduce astropy logger noise during later operations
    try:
        from astropy.logger import log
        log.setLevel('ERROR')
        from astropy.utils import iers
        iers.conf.auto_download = False
        iers.conf.auto_max_age = None
    except Exception:
        # If astropy settings change or aren't present, just continue.
        pass

    # 2) Login using credentials file if available
    if credentials_file is None:
        credentials_file = os.path.expanduser('~/.gaia_credentials')

    if os.path.isfile(credentials_file):
        try:
            Gaia.login(credentials_file=credentials_file)
        except Exception as e:
            print(f"Warning: Gaia login failed: {e}")
    else:
        print(f"Note: No Gaia credentials file found at {credentials_file}")
        print("      Create one with username on line 1, password on line 2")

    # 3) Service probe: distinguish "assets unavailable"
    if probe_service:
        try:
            # Minimal, fast probe (hits TAP briefly):
            Gaia.launch_job("SELECT TOP 1 source_id FROM gaiadr3.gaia_source", dump_to_file=False)
            # Alternatively, if you want absolutely minimal probing:
            # _ = Gaia.tap
        except Exception as e:
            raise RuntimeError(
                "Gaia external assets/services appear unavailable or unreachable."
            ) from e

    return Gaia


# --------------------------------------------------------------------------------
# Utilities that do not require astroquery
# --------------------------------------------------------------------------------
def random_rows(tab, nrows, seed=None):
    """Randomly select rows from an Astropy Table without replacement.

    Parameters
    ----------
    tab : astropy.table.Table
        Input table from which to select rows.
    nrows : int
        Number of rows to randomly select. Must be less than or equal to
        the length of the table.
    seed : int, optional
        Random seed for reproducibility. If None, uses system entropy.

    Returns
    -------
    subtab : astropy.table.Table
        New table containing the randomly selected rows.

    Warnings
    --------
    If ``nrows`` exceeds the table length, a warning is printed and the
    entire table is returned unchanged.

    Examples
    --------
    Select 10 random rows from a catalog::

        >>> from astropy.table import Table
        >>> catalog = Table.read('gaia_stars.fits')
        >>> random_subset = random_rows(catalog, 10, seed=42)
        >>> len(random_subset)
        10

    """
    if nrows > len(tab):
        print("Requested more rows than available in table")
        return tab
    rng = np.random.default_rng(seed)
    indices = rng.choice(len(tab), size=nrows, replace=False)
    return tab[indices]


def unique_rows_within_tol(tab, tol=0.01):
    """Return unique rows based on approximate equality of position and size.

    This function identifies and removes duplicate rows where RA, Dec, and Size
    are within a specified tolerance, keeping one representative row per group.

    Parameters
    ----------
    tab : astropy.table.Table
        Table containing at minimum the columns 'RA', 'Dec', and 'Size',
        all in degrees.
    tol : float, optional
        Matching tolerance in degrees. Default is 0.01° (~36 arcseconds).

    Returns
    -------
    unique_tab : astropy.table.Table
        New table containing one representative row per unique group.

    Notes
    -----
    The algorithm uses a greedy approach: rows are processed sequentially,
    and the first row in each group is kept as the representative.

    Examples
    --------
    Remove near-duplicate observations::

        >>> catalog = Table.read('observations.fits')
        >>> unique_catalog = unique_rows_within_tol(catalog, tol=0.001)
        >>> print(f"Reduced from {len(catalog)} to {len(unique_catalog)} rows")

    """
    # Stack RA, Dec, Size into a NumPy array
    data = np.vstack([tab['RA'], tab['Dec'], tab['Size']]).T

    # Initialize list of unique rows
    unique_indices = []
    used = np.zeros(len(data), dtype=bool)
    for i in range(len(data)):
        if used[i]:
            continue
        diff = np.abs(data - data[i])
        mask = np.all(diff < tol, axis=1)
        used[mask] = True
        unique_indices.append(i)

    return tab[unique_indices]


# --------------------------------------------------------------------------------
# Gaia XP spectrum helper (does not use astroquery)
# --------------------------------------------------------------------------------
def old_get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./GaiaSpec', redo=False):
    """Load or download a GAIA XP spectrum and convert to physical units.

    This function retrieves the XP continuous spectrum for a given GAIA source,
    caching it locally for future use. Spectra are converted from the archive
    format to erg/s/cm²/Å.

    Parameters
    ----------
    gaiaID : int or str
        GAIA DR3 source identifier.
    GAIA_CACHE_DIR : str, optional
        Directory path for caching spectrum files. Default is './GaiaSpec'.
        Will be created if it doesn't exist.
    redo : bool, optional
        If True, re-download the spectrum even if cached. Default is False.

    Returns
    -------
    wave : numpy.ndarray
        Wavelength array in Angstroms.
    flux : numpy.ndarray
        Flux array in erg/s/cm²/Å, or empty list ``[]`` if the spectrum
        cannot be retrieved.

    Notes
    -----
    This function was adapted from the lvmdrp package. It handles two different
    CSV formats produced by different versions of gaiaxpy (2.1.1 vs 2.1.2):
    - Newer versions use numpy representation strings
    - Older versions use comma-separated values

    The function requires the ``requests`` and ``gaiaxpy`` packages to download
    new spectra, but can read cached spectra without these dependencies.

    .. warning::
       This is the legacy spectrum retrieval function. Consider using
       :func:`get_gaia_spec` for the current recommended approach.

    Examples
    --------
    Retrieve and plot a GAIA spectrum::

        >>> import matplotlib.pyplot as plt
        >>> wave, flux = old_get_gaia_spec(4658615927801509760)
        >>> plt.plot(wave, flux)
        >>> plt.xlabel('Wavelength (Å)')
        >>> plt.ylabel('Flux (erg/s/cm²/Å)')
        >>> plt.show()

    """
    # create cache dir if it does not exist
    pathlib.Path(GAIA_CACHE_DIR).mkdir(parents=True, exist_ok=True)

    flux_path = f"{GAIA_CACHE_DIR}/gaia_spec_{gaiaID}.csv"
    wave_path = f"{GAIA_CACHE_DIR}/gaia_spec_{gaiaID}_sampling.csv"

    if path.exists(flux_path) and path.exists(wave_path) and redo == False:
        print('Star is in cache')
        gaiaflux = Table.read(flux_path, format="csv")
        gaiawave = Table.read(wave_path, format="csv")
    else:
        print('Star must be retrieved')
        # Deferred imports to keep module import-safe
        import requests
        from gaiaxpy import calibrate

        # need to download from Gaia archive
        CSV_URL = (
            "https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=XP_CONTINUOUS&ID=Gaia+DR3+"
            + str(gaiaID)
            + "&format=CSV&DATA_STRUCTURE=RAW"
        )
        FILE = f"{GAIA_CACHE_DIR}/XP_{gaiaID}_RAW.csv"

        with requests.get(CSV_URL, stream=True) as r:
            r.raise_for_status()
            if len(r.content) < 2:
                return []
            with open(FILE, "w") as f:
                f.write(r.content.decode("utf-8"))

        # convert coefficients to sampled spectrum
        _, _ = calibrate(
            FILE,
            output_path=GAIA_CACHE_DIR,
            output_file=f"gaia_spec_{gaiaID}",
            output_format="csv",
        )

        # read the flux and wavelength tables
        gaiaflux = Table.read(flux_path, format="csv")
        gaiawave = Table.read(wave_path, format="csv")

    # make numpy arrays from gaia tables
    # Handle both old CSV format (comma-separated) and new format (numpy repr strings)
    wave_str = gaiawave["pos"][0]
    flux_str = gaiaflux["flux"][0]
    
    # Try to parse as numpy array representation first (new format)
    try:
        # New format has strings like "(np.float64(1.23), np.float64(4.56), ...)"
        # Extract numeric values from inside np.float64() calls
        import re
        wave_numbers = re.findall(r'np\.float64\(([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\)', wave_str)
        flux_numbers = re.findall(r'np\.float64\(([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\)', flux_str)
        
        if len(wave_numbers) > 10 and len(flux_numbers) > 10:  # Sanity check
            wave = np.array([float(x) for x in wave_numbers]) * 10  # Angstrom
            flux = 1e4 * np.array([float(x) for x in flux_numbers])  # W/s/nm -> erg/s/cm^2/Å
        else:
            raise ValueError("Could not extract enough numeric values")
    except (ValueError, AttributeError):
        # Fall back to old CSV format (comma-separated)
        wave = np.fromstring(wave_str.strip("()"), sep=",") * 10  # Angstrom
        flux = 1e4 * np.fromstring(flux_str.strip("()"), sep=",")  # W/s/nm -> erg/s/cm^2/Å

    return wave, flux


def get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./GaiaSpec', redo=False):
    """Retrieve a GAIA XP spectrum with robust version handling.

    This is the recommended function for retrieving GAIA XP spectra. It provides
    improved error handling and compatibility with different gaiaxpy versions.

    Parameters
    ----------
    gaiaID : int or str
        GAIA DR3 source identifier.
    GAIA_CACHE_DIR : str, optional
        Directory path for caching spectrum files. Default is './GaiaSpec'.
        Will be created if it doesn't exist.
    redo : bool, optional
        If True, re-download the spectrum even if cached. Default is False.

    Returns
    -------
    wave : numpy.ndarray
        Wavelength array in Angstroms.
    flux : numpy.ndarray
        Flux array in erg/s/cm²/Å, or empty list if retrieval fails.

    Notes
    -----
    This function wraps :func:`old_get_gaia_spec` with additional error handling
    to gracefully handle failures in spectrum retrieval or conversion.

    Examples
    --------
    Retrieve a spectrum with error handling::

        >>> result = get_gaia_spec(4658615927801509760)
        >>> if isinstance(result, list):
        ...     print("Spectrum retrieval failed")
        ... else:
        ...     wave, flux = result
        ...     print(f"Retrieved spectrum with {len(wave)} points")

    """
    try:
        wave, flux = old_get_gaia_spec(gaiaID, GAIA_CACHE_DIR=GAIA_CACHE_DIR, redo=redo)
        return wave, flux
    except Exception as e:
        print(f'Error retrieving spectrum for GAIA ID {gaiaID}: {e}')
        return []


def get_gaia_spectra_batch(source_ids, GAIA_CACHE_DIR='./GaiaSpec', redo=False):
    """Retrieve XP spectra for multiple Gaia sources in a single archive call.

    Uses gaiaxpy.calibrate() on the full list, then splits and caches results
    in the same per-star format used by get_gaia_spec(), so subsequent calls
    to get_gaia_spec() for the same IDs will be served from cache.

    Parameters
    ----------
    source_ids : list of int
        Gaia DR3 source identifiers.
    GAIA_CACHE_DIR : str, optional
        Cache directory. Default is './GaiaSpec'.
    redo : bool, optional
        If False (default), skip IDs already in cache.

    Returns
    -------
    n_ok : int
        Number of spectra successfully retrieved and cached.
    n_fail : int
        Number of IDs for which retrieval failed.
    """
    import warnings
    import pandas as pd
    from gaiaxpy import calibrate

    pathlib.Path(GAIA_CACHE_DIR).mkdir(parents=True, exist_ok=True)

    if not redo:
        pending = [sid for sid in source_ids
                   if not (path.exists(f'{GAIA_CACHE_DIR}/gaia_spec_{sid}.csv') and
                           path.exists(f'{GAIA_CACHE_DIR}/gaia_spec_{sid}_sampling.csv'))]
        n_cached = len(source_ids) - len(pending)
        if n_cached:
            print(f'get_gaia_spectra_batch: {n_cached} already cached, fetching {len(pending)}')
    else:
        pending = list(source_ids)

    if not pending:
        print('get_gaia_spectra_batch: all spectra already cached')
        return len(source_ids), 0

    # Read credentials if available
    username, password = None, None
    cred_file = os.path.expanduser('~/.gaia_credentials')
    if os.path.isfile(cred_file):
        with open(cred_file) as f:
            lines = f.read().splitlines()
        if len(lines) >= 2:
            username, password = lines[0].strip(), lines[1].strip()

    print(f'get_gaia_spectra_batch: retrieving {len(pending)} spectra in one call')
    spectra = None
    for attempt in range(3):
        try:
            if attempt > 0:
                print(f'get_gaia_spectra_batch: retry {attempt}/2...')
                time.sleep(5)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='.*archive is unstable.*')
                spectra, sampling = calibrate(pending, save_file=False,
                                              username=username, password=password)
            break
        except Exception as e:
            print(f'get_gaia_spectra_batch: attempt {attempt+1} failed: {e}')
    if spectra is None:
        return 0, len(pending)

    # sampling is a numpy array of wavelength positions (same for all stars)
    import numpy as np
    if isinstance(sampling, pd.DataFrame):
        wave_arr = sampling['pos'].iloc[0] if 'pos' in sampling.columns else sampling.iloc[0].values
    else:
        wave_arr = np.asarray(sampling)
    wave_str = str(tuple(wave_arr))

    # identify flux column name
    flux_col = 'flux' if 'flux' in spectra.columns else [c for c in spectra.columns if c != 'source_id'][0]

    n_ok = 0
    n_fail = 0
    retrieved_ids = set(spectra['source_id'].values)

    for sid in pending:
        if sid not in retrieved_ids:
            n_fail += 1
            continue
        row = spectra[spectra['source_id'] == sid].iloc[0]
        flux_path = f'{GAIA_CACHE_DIR}/gaia_spec_{sid}.csv'
        wave_path = f'{GAIA_CACHE_DIR}/gaia_spec_{sid}_sampling.csv'
        try:
            pd.DataFrame({'flux': [str(tuple(row[flux_col]))]}).to_csv(flux_path, index=False)
            pd.DataFrame({'pos': [wave_str]}).to_csv(wave_path, index=False)
            n_ok += 1
        except Exception as e:
            print(f'get_gaia_spectra_batch: failed caching {sid}: {e}')
            n_fail += 1

    print(f'get_gaia_spectra_batch: {n_ok} cached, {n_fail} failed')
    return n_ok, n_fail


def get_gaia_from_file(ra=84.92500000000001, dec=-66.27416666666667,
                       size_deg=0.3, filename='Gaia_MagClouds.fits', outroot=''):
    """Extract GAIA sources from a local catalog file within a sky region.

    This function performs a rectangular selection from a pre-downloaded GAIA
    catalog file, useful when working offline or with large local catalogs.

    Parameters
    ----------
    ra : float, optional
        Right Ascension of the field center in degrees (J2000).
        Default is 84.925 (approximately LMC).
    dec : float, optional
        Declination of the field center in degrees (J2000).
        Default is -66.274 (approximately LMC).
    size_deg : float, optional
        Size of the extraction region in degrees. Default is 0.3°.
    filename : str, optional
        Name of the local GAIA catalog file. Searches locally first, then
        in $KRED/xdata/. Default is 'Gaia_MagClouds.fits'.
    outroot : str, optional
        Root name for the output file. If empty, constructs from RA and Dec.

    Returns
    -------
    outfile : str
        Path to the output FITS file containing the extracted sources.

    Raises
    ------
    IOError
        If the input file cannot be located in either the current directory
        or $KRED/xdata/.

    Notes
    -----
    The function performs a rectangular (RA, Dec) selection rather than a
    cone search. The RA range is adjusted for declination to approximate
    equal angular sizes in both dimensions.

    Output files are written to a ``Gaia/`` subdirectory in FITS format.

    Examples
    --------
    Extract sources from a local catalog::

        >>> outfile = get_gaia_from_file(
        ...     ra=150.0, dec=-30.0,
        ...     size_deg=0.5,
        ...     filename='gaia_dr3_subset.fits'
        ... )
        >>> print(f"Extracted catalog: {outfile}")

    """
    # 1) Determine input file path - search locally first, then $KRED/xdata/
    xfilename = ''
    if os.path.isfile(filename):
        xfilename = filename
        print(f'get_gaia_from_file: Using local file {xfilename}')
    else:
        KRED = os.environ.get("KRED")
        if KRED is not None:
            candidate = f"{KRED}/xdata/{filename}"
            if os.path.isfile(candidate):
                xfilename = candidate
                print(f'get_gaia_from_file: Using {xfilename}')
            else:
                raise IOError(f'Could not locate {filename} locally or in $KRED/xdata/')
        else:
            raise IOError(f'Could not locate {filename} locally and KRED environment variable is not set')

    # read the local table
    if xfilename.lower().endswith('.fits'):
        xtab = Table.read(xfilename)
    else:
        xtab = ascii.read(xfilename)

    # rectangular selection around (ra, dec)
    dec_min = dec - 0.5 * size_deg
    dec_max = dec + 0.5 * size_deg
    xscale = np.cos(dec / 57.29578)
    factor = 0.5 * size_deg / xscale
    ra_min = ra - factor
    ra_max = ra + factor

    mask = ((ra_min < xtab['RA']) & (xtab['RA'] < ra_max) &
            (dec_min < xtab['Dec']) & (xtab['Dec'] < dec_max))
    ftab = xtab[mask]

    if outroot == '':
        outroot = '%06.2f_%06.2f' % (ra, dec)
    os.makedirs('Gaia', exist_ok=True)
    outfile = f'Gaia/Gaia.{outroot}.fits'

    ftab.write(outfile, format='fits', overwrite=True)
    print(f'Wrote {outfile} with {len(ftab)} objects')
    return outfile


# --------------------------------------------------------------------------------
# Gaia archive access (astroquery used only inside these functions via load_Gaia)
# --------------------------------------------------------------------------------
def get_gaia_from_archive(ra=84.92500000000001, dec=-66.27416666666667,
                          rad_deg=0.3, outroot='', nmax=-1, redo=False,
                          max_retries=3, retry_delay=5):
    """Query the GAIA archive with cone search and automatic retry on errors.

    This function queries the GAIA DR3 archive for photometric data within a
    specified cone, with robust handling of network errors through automatic
    retries.

    Parameters
    ----------
    ra : float, optional
        Right Ascension of the cone center in degrees (J2000).
        Default is 84.925.
    dec : float, optional
        Declination of the cone center in degrees (J2000).
        Default is -66.274.
    rad_deg : float, optional
        Cone search radius in degrees. Default is 0.3°.
    outroot : str, optional
        Root name for the output file. If empty, constructs from RA and Dec.
    nmax : int, optional
        Maximum number of rows to return. If -1 (default), returns all
        matching sources.
    redo : bool, optional
        If True, re-query even if the output file exists. Default is False.
    max_retries : int, optional
        Maximum number of retry attempts for failed queries. Default is 3.
    retry_delay : float, optional
        Delay in seconds between retry attempts. Default is 5.

    Returns
    -------
    outfile : str
        Path to the output file containing the catalog, or empty list if
        no sources were retrieved.

    Raises
    ------
    RuntimeError
        If astroquery is not installed or GAIA services are unavailable.
    IncompleteRead
        If max retries are exceeded due to persistent network errors.

    Notes
    -----
    The function automatically renames columns to a simplified convention:
    - ``ra`` → ``RA``
    - ``dec`` → ``Dec``
    - ``source_id`` → ``Source_name``
    - ``phot_g_mean_mag`` → ``G``
    - ``phot_bp_mean_mag`` → ``B``
    - ``phot_rp_mean_mag`` → ``R``
    - ``teff_gspphot`` → ``teff``
    - ``logg_gspphot`` → ``log_g``
    - ``distance_gspphot`` → ``D``
    - ``has_xp_continuous`` → ``xp_spec_exists`` (boolean flag)

    Output is written to ``Gaia/Gaia.<outroot>.fits`` in FITS table format.

    Examples
    --------
    Basic cone search with default retry behavior::

        >>> outfile = get_gaia_from_archive(
        ...     ra=180.0, dec=45.0,
        ...     rad_deg=0.5,
        ...     outroot='ngc1234'
        ... )
        >>> print(f"Catalog written to: {outfile}")

    Query with custom retry parameters::

        >>> outfile = get_gaia_from_archive(
        ...     ra=10.0, dec=-5.0,
        ...     rad_deg=0.1,
        ...     max_retries=5,
        ...     retry_delay=10
        ... )

    """
    if outroot == '':
        outroot = '%05.1f_%05.1f' % (ra, dec)
    os.makedirs('Gaia', exist_ok=True)
    outfile = 'Gaia/Gaia.%s.fits' % outroot

    if not redo and os.path.isfile(outfile):
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile

    print('get_gaia: Getting data for RA Dec of %.5f %.5f and size of %.2f' % (ra, dec, rad_deg))

    r = None
    start_time = time.time()

    # Try MAST first — synchronous, avoids ESA async job storage bugs
    try:
        print('get_gaia: Querying MAST GaiaDR3...')
        coord = SkyCoord(ra=ra, dec=dec, unit=u.degree, frame='icrs')
        r = Catalogs.query_region(coord, radius='%f deg' % rad_deg, catalog='GaiaDR3')
        if nmax > 0:
            r = r[:nmax]
        print('get_gaia: MAST query returned %d rows.' % len(r))
    except Exception as e:
        print(f'get_gaia: MAST query failed ({type(e).__name__}: {e}), falling back to ESA archive...')
        r = None

    # Fall back to ESA archive with retry loop
    if r is None:
        try:
            Gaia = load_Gaia(probe_service=False)
        except RuntimeError as err:
            raise RuntimeError(f"Gaia archive access failed: {err}") from err

        Gaia.ROW_LIMIT = nmax
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

        for attempt in range(max_retries):
            try:
                if attempt > 0:
                    print(f'get_gaia: Retry attempt {attempt + 1}/{max_retries}...')
                print('get_gaia: Submitting query to ESA Gaia archive...')
                j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))
                print('get_gaia: Waiting for results...')
                r = j.get_results()
                print('get_gaia: Results received.')
                break
            except IncompleteRead as e:
                print(f'get_gaia: IncompleteRead error on attempt {attempt + 1}: {e}')
                if attempt < max_retries - 1:
                    print(f'get_gaia: Retrying in {retry_delay} seconds...')
                    time.sleep(retry_delay)
                else:
                    print('get_gaia: Max retries reached. Query failed.')
                    raise
            except Exception as e:
                print(f'get_gaia: Unexpected error on attempt {attempt + 1}: {type(e).__name__}: {e}')
                if attempt < max_retries - 1:
                    print(f'get_gaia: Retrying in {retry_delay} seconds...')
                    time.sleep(retry_delay)
                else:
                    print('get_gaia: Max retries reached. Query failed.')
                    raise

    elapsed_time = time.time() - start_time
    print(f'get_gaia: Query completed in {elapsed_time:.1f} seconds')

    if r is None or len(r) == 0:
        print('Error: get_gaia: No objects were retrieved')
        return []

    # Process and rename columns
    r.rename_column('ra', 'RA')
    r.rename_column('dec', 'Dec')
    r['RA'].unit = None
    r['Dec'].unit = None
    try:
        r.rename_column('source_id', 'Source_name')
    except Exception:
        r.rename_column('SOURCE_ID', 'Source_name')
    r['Source_name'] = np.array(r['Source_name'], dtype='int64')
    r.rename_column('phot_g_mean_mag', 'G')
    r.rename_column('phot_bp_mean_mag', 'B')
    r.rename_column('phot_rp_mean_mag', 'R')
    r.rename_column('teff_gspphot', 'teff')
    r.rename_column('logg_gspphot', 'log_g')
    r.rename_column('distance_gspphot', 'D')

    r['log_g'].unit = None  # strip unparseable 'log(cm.s**-2)' unit
    cols = ['Source_name', 'RA', 'Dec', 'B', 'G', 'R', 'teff', 'log_g', 'D']
    if 'has_xp_continuous' in r.colnames:
        r.rename_column('has_xp_continuous', 'xp_spec_exists')
        cols.append('xp_spec_exists')
    r[cols].write(outfile, format='fits', overwrite=True)
    n_spec = int(np.sum(r['xp_spec_exists'])) if 'xp_spec_exists' in r.colnames else 0
    print('Wrote %s with %d objects (%d with XP spectra)' % (outfile, len(r), n_spec))
    return outfile


def get_gaia(ra,dec,size):
    """
    Retrieve data from a file if possible, but if it is not possible try
    to get the data from the archive
    """

    try:
        outfile=get_gaia_from_file(ra,dec,size)
    except:
        print('Could not get data locally, so trying the GaiaArchive')
        outfile=get_gaiat_from_archive(ra,dec,size)

    return outfile
    
        



def get_gaia_from_archive_old(ra=84.92500000000001, dec=-66.27416666666667,
                               rad_deg=0.3, outroot='', nmax=-1, redo=False):
    """Query the GAIA archive with cone search (legacy version without retry).

    .. deprecated:: 251211
       Use :func:`get_gaia_from_archive` instead, which includes retry logic
       for improved reliability.

    This is the original cone search function without automatic retry handling.
    It is retained for backward compatibility but is not recommended for new code.

    Parameters
    ----------
    ra : float, optional
        Right Ascension of the cone center in degrees (J2000).
        Default is 84.925.
    dec : float, optional
        Declination of the cone center in degrees (J2000).
        Default is -66.274.
    rad_deg : float, optional
        Cone search radius in degrees. Default is 0.3°.
    outroot : str, optional
        Root name for the output file. If empty, constructs from RA and Dec.
    nmax : int, optional
        Maximum number of rows to return. If -1 (default), returns all sources.
    redo : bool, optional
        If True, re-query even if the output file exists. Default is False.

    Returns
    -------
    outfile : str
        Path to the output file, or empty list if no sources retrieved.

    Raises
    ------
    RuntimeError
        If astroquery is not installed or GAIA services are unavailable.

    See Also
    --------
    get_gaia_from_archive : Recommended function with retry logic

    """
    try:
        Gaia = load_Gaia(probe_service=True)
    except RuntimeError as err:
        raise RuntimeError(f"Gaia archive access failed: {err}") from err

    if outroot == '':
        outroot = '%06.2f_%06.2f' % (ra, dec)
    os.makedirs('Gaia', exist_ok=True)
    outfile = 'Gaia/Gaia.%s.txt' % outroot

    if not redo and os.path.isfile(outfile):
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile

    print('get_gaia: Getting data for RA Dec of %.5f %.5f and size of %.2f' % (ra, dec, rad_deg))
    Gaia.ROW_LIMIT = nmax  # Ensure the default row limit.
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))
    r = j.get_results()

    if len(r) == 0:
        print('Error: get_gaia: No objects were retrieved')
        return []

    r.rename_column('ra', 'RA')
    r.rename_column('dec', 'Dec')
    try:
        r.rename_column('source_id', 'Source_name')
    except Exception:
        r.rename_column('SOURCE_ID', 'Source_name')
    r.rename_column('phot_g_mean_mag', 'G')
    r.rename_column('phot_bp_mean_mag', 'B')
    r.rename_column('phot_rp_mean_mag', 'R')
    r.rename_column('teff_gspphot', 'teff')
    r.rename_column('logg_gspphot', 'log_g')
    r.rename_column('distance_gspphot', 'D')

    r['Source_name', 'RA', 'Dec', 'B', 'G', 'R', 'teff', 'log_g', 'D'].write(
        outfile, format='ascii.fixed_width_two_line', overwrite=True
    )
    print('Wrote %s with %d objects' % (outfile, len(r)))
    return outfile


def simple_test():
    """Run basic functionality tests for GaiaCat module.

    This function performs a simple test of core functionality:
    1. Retrieves a small catalog from the GAIA archive
    2. Retrieves a spectrum for a known source

    Returns
    -------
    None

    Notes
    -----
    This test requires network access and will create files in the
    ``Gaia/`` and ``GaiaSpec/`` directories.

    Examples
    --------
    Run the test suite::

        >>> simple_test()
        Here we just check if a small amount of the SW works
        Can we retrieve a catalog of stars
        ...
        If there were no errors we have success

    """
    print('Here we just check if a small amount of the SW works')
    print('Can we retrieve a catalog of stars')
    if os.path.isfile('Gaia/Gaia.foo.txt'):
        os.remove('Gaia/Gaia.foo.txt')
        print('Removed Gaia/Gaia.foo.txt, before retrieving catalog from archive')
    get_gaia_from_archive(rad_deg=0.1, outroot='foo')
    print('Can we retrieve a spectrum')
    get_gaia_spec(4658615927801509760, redo=True)
    print('\nIf there were no errors we have success\n')
    return


def get_wcs_center(fitsfile):
    """Extract the center RA/Dec from a FITS file WCS.

    Parameters
    ----------
    fitsfile : str
        Path to FITS file with valid WCS in header

    Returns
    -------
    ra : float
        Right Ascension of image center in degrees
    dec : float
        Declination of image center in degrees

    Raises
    ------
    ValueError
        If WCS cannot be extracted from the FITS file
    """
    from astropy.wcs import WCS

    with fits.open(fitsfile) as hdul:
        header = hdul[0].header
        data_shape = hdul[0].data.shape

        try:
            wcs = WCS(header)
        except Exception as e:
            raise ValueError(f"Could not extract WCS from {fitsfile}: {e}")

        # Get center pixel
        ny, nx = data_shape
        cx, cy = nx / 2.0, ny / 2.0

        # Convert to sky coordinates
        ra, dec = wcs.wcs_pix2world(cx, cy, 0)

    return float(ra), float(dec)


# --------------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------------
def steer(argv):
    """Execute command-line interface for GaiaCat module.

    Parses command-line arguments and retrieves GAIA catalog data.
    By default, tries local file first, then falls back to archive.

    Parameters
    ----------
    argv : list of str
        Command-line arguments

    Returns
    -------
    str
        Path to output file, or None if error
    """
    force_archive = False
    redo = False
    rad_deg = 0.5
    outroot = ''
    ra = None
    dec = None
    fitsfile = None
    gaia_file = 'Gaia_MagClouds.fits'

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-archive':
            force_archive = True
        elif argv[i] == '-redo':
            redo = True
            force_archive = True
        elif argv[i][:6] == '-gfile':
            i += 1
            gaia_file = argv[i]
        elif argv[i][:4] == '-rad':
            i += 1
            rad_deg = float(argv[i])
        elif argv[i][:4] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i][0] == '-' and ra is not None and dec is None:
            try:
                dec = float(argv[i])
            except ValueError:
                print('Error: Unknown option:', argv[i])
                return
        elif argv[i][0] == '-':
            print('Error: Unknown option:', argv[i])
            return
        elif ra is None:
            # First positional argument - could be FITS file or RA
            if argv[i].endswith('.fits') or argv[i].endswith('.fits.gz'):
                fitsfile = argv[i]
            else:
                try:
                    ra = float(argv[i])
                except ValueError:
                    # Assume it's a FITS file without .fits extension
                    if os.path.isfile(argv[i]):
                        fitsfile = argv[i]
                    else:
                        print(f'Error: Cannot parse {argv[i]} as RA or find as file')
                        return
        elif dec is None:
            try:
                dec = float(argv[i])
            except ValueError:
                print(f'Error: Cannot parse {argv[i]} as Dec')
                return
        else:
            print('Error: Too many arguments:', argv[i])
            return
        i += 1

    # If FITS file provided, extract RA/Dec from WCS
    if fitsfile is not None:
        if not os.path.isfile(fitsfile):
            print(f'Error: FITS file not found: {fitsfile}')
            return
        try:
            ra, dec = get_wcs_center(fitsfile)
            print(f'Extracted center from {fitsfile}: RA={ra:.5f}, Dec={dec:.5f}')
        except Exception as e:
            print(f'Error extracting WCS: {e}')
            return

        # Default outroot from FITS filename
        if outroot == '':
            outroot = os.path.basename(fitsfile).replace('.fits.gz', '').replace('.fits', '')

    # Validate we have coordinates
    if ra is None or dec is None:
        print('Error: Must provide either a FITS file or RA and Dec')
        print(__doc__)
        return

    # Default outroot from coordinates
    if outroot == '':
        outroot = '%.2f_%.2f' % (ra, dec)

    # Print parameters
    print('         RA : %.5f' % ra)
    print('        Dec : %.5f' % dec)
    print('     Radius : %.3f deg' % rad_deg)
    print('    Outroot : %s' % outroot)
    print('  Gaia file : %s' % gaia_file)

    # Try local file first (unless -archive specified), then fall back to archive
    if force_archive:
        print('Source mode : archive (forced)')
        outfile = get_gaia_from_archive(ra=ra, dec=dec, rad_deg=rad_deg, outroot=outroot, redo=redo)
    else:
        print('Source mode : local file (with archive fallback)')
        try:
            outfile = get_gaia_from_file(ra=ra, dec=dec, size_deg=rad_deg * 2,
                                         filename=gaia_file, outroot=outroot)
        except (IOError, FileNotFoundError) as e:
            print(f'Local file not available: {e}')
            print('Falling back to archive...')
            outfile = get_gaia_from_archive(ra=ra, dec=dec, rad_deg=rad_deg, outroot=outroot, redo=redo)

    return outfile


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
