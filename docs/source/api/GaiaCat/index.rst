GaiaCat
=======

.. py:module:: GaiaCat

.. autoapi-nested-parse::

   GaiaCat - GAIA Database Interaction Module.

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

   This module provides functions to retrieve GAIA DR3 photometric and
   spectroscopic data.  Catalog cone searches query MAST first
   (``astroquery.mast.Catalogs``) because it is synchronous and avoids the
   async-job storage bugs on the ESA server; ESA (``astroquery.gaia``) is
   retained as a fallback.  The primary functions are:

   - get_gaia_from_archive(): Query MAST (primary) or ESA (fallback) with
     cone search
   - get_gaia_from_file(): Extract from a local pre-downloaded catalog
   - get_gaia_spectra_batch(): Bulk-retrieve Gaia XP spectra for a list of
     source IDs, pre-populating the per-star cache used by get_gaia_spec()

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
   260429 ksl Switch catalog queries to MAST as primary source; ESA async
              endpoint retained as fallback.  Add get_gaia_spectra_batch()
              for bulk XP spectrum retrieval.
   260502 ksl Add redo parameter to get_gaia_from_file; skip re-extraction
              when output file already exists.

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



Functions
---------

.. autoapisummary::

   GaiaCat.get_gaia
   GaiaCat.get_gaia_from_archive
   GaiaCat.get_gaia_from_archive_old
   GaiaCat.get_gaia_from_file
   GaiaCat.get_gaia_spec
   GaiaCat.get_gaia_spectra_batch
   GaiaCat.get_wcs_center
   GaiaCat.load_Gaia
   GaiaCat.old_get_gaia_spec
   GaiaCat.random_rows
   GaiaCat.simple_test
   GaiaCat.steer
   GaiaCat.unique_rows_within_tol


Module Contents
---------------

.. py:function:: get_gaia(ra, dec, size)

   Retrieve data from a file if possible, but if it is not possible try
   to get the data from the archive


.. py:function:: get_gaia_from_archive(ra=84.92500000000001, dec=-66.27416666666667, rad_deg=0.3, outroot='', nmax=-1, redo=False, max_retries=3, retry_delay=5)

   Query the GAIA archive with cone search and automatic retry on errors.

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



.. py:function:: get_gaia_from_archive_old(ra=84.92500000000001, dec=-66.27416666666667, rad_deg=0.3, outroot='', nmax=-1, redo=False)

   Query the GAIA archive with cone search (legacy version without retry).

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



.. py:function:: get_gaia_from_file(ra=84.92500000000001, dec=-66.27416666666667, size_deg=0.3, filename='Gaia_MagClouds.fits', outroot='', redo=False)

   Extract GAIA sources from a local catalog file within a sky region.

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
   redo : bool, optional
       If True, re-extract and overwrite the output file even if it already
       exists. Default is False.

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



.. py:function:: get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./GaiaSpec', redo=False)

   Retrieve a GAIA XP spectrum with robust version handling.

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



.. py:function:: get_gaia_spectra_batch(source_ids, GAIA_CACHE_DIR='./GaiaSpec', redo=False)

   Retrieve XP spectra for multiple Gaia sources in a single archive call.

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


.. py:function:: get_wcs_center(fitsfile)

   Extract the center RA/Dec from a FITS file WCS.

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


.. py:function:: load_Gaia(probe_service=True, credentials_file=None)

   Return the astroquery.gaia.Gaia class with service validation.

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



.. py:function:: old_get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./GaiaSpec', redo=False)

   Load or download a GAIA XP spectrum and convert to physical units.

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



.. py:function:: random_rows(tab, nrows, seed=None)

   Randomly select rows from an Astropy Table without replacement.

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



.. py:function:: simple_test()

   Run basic functionality tests for GaiaCat module.

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



.. py:function:: steer(argv)

   Execute command-line interface for GaiaCat module.

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


.. py:function:: unique_rows_within_tol(tab, tol=0.01)

   Return unique rows based on approximate equality of position and size.

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



