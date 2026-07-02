Smash
=====

.. py:module:: Smash

.. autoapi-nested-parse::

   Smash - SMASH DR2 Catalog Retrieval

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



Attributes
----------

.. autoapisummary::

   Smash.SMASH_CACHE_DIR


Functions
---------

.. autoapisummary::

   Smash.check_cache
   Smash.clear_preassembled_cache
   Smash.do_one
   Smash.get_cache_filename
   Smash.get_smash
   Smash.get_smash_from_file
   Smash.plot_positions
   Smash.save_to_cache
   Smash.select_fraction_smash_dr2
   Smash.smash_cone_search
   Smash.steer


Module Contents
---------------

.. py:data:: SMASH_CACHE_DIR
   :value: 'Smash'


.. py:function:: check_cache(ra, dec, radius, tol=0.0001)

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


.. py:function:: clear_preassembled_cache()

   Release the module-level SMASH catalog table from memory.

   Call this after pre-caching all tile files and before launching a
   multiprocessing pool so forked workers inherit a lean parent process
   rather than a copy of the 15 GB table.


.. py:function:: do_one(ra, dec, radius=0.5, outroot='smash_cat', rmag_max=22.0, keep_frac=0.5, plot=False)

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


.. py:function:: get_cache_filename(ra, dec, radius)

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


.. py:function:: get_smash(ra, dec, size, rmag_max=22.0, keep_frac=0.5)

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


.. py:function:: get_smash_from_file(ra, dec, size_deg, filename='Smash_MagClouds.fits', outroot='', rmag_max=22.0, keep_frac=0.5)

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


.. py:function:: plot_positions(table, outfile='smash_positions.png', title=None, rmag_max=None, keep_frac=None, n_original=None)

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


.. py:function:: save_to_cache(table, ra, dec, radius)

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


.. py:function:: select_fraction_smash_dr2(tab, rmag_max, keep_frac, rmag_ref=21.0, sigma_chi=0.5, sharp0=0.05, eps=0.001)

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


.. py:function:: smash_cone_search(ra, dec, radius_deg, limit=10000000, verbose=True, use_cache=True)

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


.. py:function:: steer(argv)

   Command line interface for SMASH catalog retrieval.

   Parses command line arguments and calls do_one to perform the
   cone search and filtering.

   Parameters
   ----------
   argv : list
       Command line arguments (typically sys.argv).


