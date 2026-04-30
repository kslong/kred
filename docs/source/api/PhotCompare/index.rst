PhotCompare
===========

.. py:module:: PhotCompare

.. autoapi-nested-parse::

   PhotCompare - Photometric Comparison Tool

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
   ----------------

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

     - Photometry results (``TabPhot/xmatch_*.txt``)
     - Cross-matched catalogs
     - Source lists (unforced mode)

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

   The xmatch output files written to ``TabPhot/`` are the primary input for
   ``ZeroPoint.py``.

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

   260430 ksl
       Handle photutils >= 2.x returning shaped arrays from
       ``ApertureStats.fwhm`` and ``ApertureStats.eccentricity``; use
       ``.flat[0]`` instead of ``float()`` conversion.

   Author
   ------
   Space Telescope Science Institute



Attributes
----------

.. autoapisummary::

   PhotCompare.XDIR


Functions
---------

.. autoapisummary::

   PhotCompare.do_dir
   PhotCompare.do_fig
   PhotCompare.do_fig_diff
   PhotCompare.do_forced_photometry
   PhotCompare.do_many
   PhotCompare.do_one
   PhotCompare.do_xphot
   PhotCompare.find_closest_objects
   PhotCompare.get_objects_from_image
   PhotCompare.get_size
   PhotCompare.locate_first_image_extension
   PhotCompare.random_rows
   PhotCompare.read_table
   PhotCompare.steer
   PhotCompare.unique_rows_within_tol


Module Contents
---------------

.. py:data:: XDIR
   :value: ''


.. py:function:: do_dir(xdir='DECam_SWARP2/LMC_c37/T16', nrows_max=30000, forced=True)

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


.. py:function:: do_fig(xtab, outroot='')

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


.. py:function:: do_fig_diff(xtab, outroot)

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


.. py:function:: do_forced_photometry(filename='LMC_c48_T08.r.t060.fits', object_file='objects.txt', nrows_max=-1, outroot='', rstar=6, b_in=8, b_out=12, add_psf_metrics=True)

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
   add_psf_metrics : bool, optional
       If True, adds columns useful for PSF star selection (SNR, Concentration,
       BkgContam). Default: False.

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

   Written to ``TabPhot{XDIR}/{outroot}_phot.fits``

   Examples
   --------
   >>> phot_file = do_forced_photometry('image.fits', 'gaia_sources.fits')
   >>> phot = ascii.read(phot_file)
   >>> print(f"Measured {len(phot)} sources")


.. py:function:: do_many(filenames=['LMC_c48_T08.r.t060.fits'], gaia_cat_file='', forced=True, nrows_max=10000, outroot='')

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


.. py:function:: do_one(filename='LMC_c48_T08.r.t060.fits', gaia_cat_file='', forced=False, nrows_max=-1, outroot='')

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


.. py:function:: do_xphot(filename, gaia_file, forced, nrows_max, outroot)

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


.. py:function:: find_closest_objects(table1_path, table2_path, max_sep=0.5)

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


.. py:function:: get_objects_from_image(filename='LMC_c48_T08.r.t060.fits', outroot='')

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


.. py:function:: get_size(filename='LMC_c48_T08.r.t060.fits')

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


.. py:function:: locate_first_image_extension(xx)

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


.. py:function:: random_rows(tab, nrows, seed=None)

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


.. py:function:: read_table(filename)

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


.. py:function:: steer(argv)

   Parse command-line arguments and run PhotCompare.

   Usage: PhotCompare.py [-h] [-dir DIRNAME] [-nmax N] [-forced] [-unforced]
                         [-gcat FILE] [-out NAME] file1 file2 ...


.. py:function:: unique_rows_within_tol(tab, tol=0.01)

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


