MefPhot
=======

.. py:module:: MefPhot

.. autoapi-nested-parse::

   MefPhot - Multi-Extension FITS Forced Photometry

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
   B_IN, B_OUT, CATALOG, FIELD, ROOT, and DETNAME, and columns including:

   * Source positions (pixel and sky coordinates)
   * Raw and background-subtracted fluxes with uncertainties
   * FWHM and eccentricity measurements
   * Instrumental magnitudes (zero point = 28)
   * Reference catalog photometry (RA, Dec, G, R, ...)
   * EXT, CCD (per-row provenance columns)

     File-level metadata is in the extension 1 header, not repeated per row:
     FILTER, EXPTIME, MAGZERO, SEEING, CATALOG, RADIUS, B_IN, B_OUT,
     FILE, FIELD, ROOT, DETNAME.

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

   2026-05-02 ksl
       Fix photutils >= 2.x compatibility: ApertureStats.fwhm and eccentricity
       now return shaped arrays; use .flat[0] to extract scalar values.

   2026-06-08 ksl
       Add FIELD, ROOT, DETNAME keywords to extension 1 header, parsed from
       the input file path, so downstream tools can identify the source MEF
       without reparsing filenames.

   Author
   ------
   Space Telescope Science Institute

   .. moduleauthor:: KSL



Functions
---------

.. autoapisummary::

   MefPhot.check_smash_memory
   MefPhot.do_forced_photometry
   MefPhot.do_many
   MefPhot.do_one
   MefPhot.get_available_memory
   MefPhot.get_total_memory
   MefPhot.precache_smash_tiles
   MefPhot.random_rows
   MefPhot.read_table
   MefPhot.steer


Module Contents
---------------

.. py:function:: check_smash_memory(n_processes, catalog, n_files=None)

   Warn and prompt when launching multiple SMASH workers would exhaust RAM.

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


.. py:function:: do_forced_photometry(filename='LMC_c48_T08.r.t060.fits', image_ext=1, object_file='objects.txt', nrows_max=-1, rstar=6, b_in=8, b_out=12, add_psf_metrics=True, mag_bright=14.0, mag_faint=22.0)

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


.. py:function:: do_many(filenames, outroot='', nrows_max=-1, rstar=6, b_in=8, b_out=12, n_processes=None, logfile='ErrorsMefPhot.txt', verbose_errors=False, catalog='gaia', redo=False, mag_bright=14.0, mag_faint=22.0)

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


.. py:function:: do_one(filename='foo.fits', outroot='', nrows_max=-1, rstar=6, b_in=8, b_out=12, catalog='gaia', verbose=True, redo=False, mag_bright=14.0, mag_faint=22.0)

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


.. py:function:: get_available_memory()

   Return available (free + reclaimable) memory in bytes, cross-platform.

   Tries psutil first (works on Linux and macOS), then falls back to
   platform-specific methods so the function works without psutil.

   Returns
   -------
   int or None
       Available memory in bytes, or None if it cannot be determined.


.. py:function:: get_total_memory()

   Return total installed physical RAM in bytes, cross-platform.

   Returns
   -------
   int or None
       Total RAM in bytes, or None if it cannot be determined.


.. py:function:: precache_smash_tiles(filenames)

   Load Smash_MagClouds.fits once and write per-tile cache files for all extensions.

   Call this before do_many() when using the SMASH catalog.  All workers
   will find their tile already cached on disk and never load the large
   pre-assembled file, keeping per-worker memory low.  After this function
   returns, call Smash.clear_preassembled_cache() to free the big table
   before the multiprocessing pool is created.

   Parameters
   ----------
   filenames : list of str
       FITS files that will be passed to do_many().


.. py:function:: random_rows(tab, nrows, seed=None)

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


.. py:function:: read_table(filename)

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


.. py:function:: steer(argv)

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


