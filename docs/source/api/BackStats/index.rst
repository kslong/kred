BackStats
=========

.. py:module:: BackStats

.. autoapi-nested-parse::

   BackStats - Background Statistics Calculator

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       BackStats.py [-np 8] [-all] [-rm] [-sigma_clipped] field T01 T02 ...

   **Required Arguments:**

   field
       Field name (e.g., 'LMC_c42')

   T01 T02 ...
       One or more tile identifiers (unless -all is specified)

   **Optional Arguments:**

   -all
       Process all tiles (T01 through T16)

   -rm
       Remove files created by BackPrep after successful completion. Only
       occurs at the end of the program so that if something fails, no
       background files will have been deleted.

   -np N
       Perform calculations in parallel using N threads (default: 1)

   -sigma_clipped
       Diagnostic mode that calculates sigma-clipped versions of mean,
       median, and std. Greatly increases run time but provides robust
       statistics.

   Output
   ------
   Results written to ``Summary/{field}_{tile}_xxx.txt`` containing:

   * Overlapping file pairs
   * Number of overlapping pixels
   * Mean, median, std of differences
   * Sigma-clipped statistics (if requested)
   * Skewness and kurtosis
   * Mode (half-sample algorithm)
   * Filter and exposure time metadata

   Performance Notes
   This routine can take considerable time because:

   * Each tile involves many files (often 100-1000)
   * Must check all pairwise combinations (N × (N-1) matches)
   * For 100 files: ~10,000 comparisons
   * For 1000 files: ~1,000,000 comparisons

   **Optimization:**

   Running with maximum available threads is strongly recommended. Use ``-np 0``
   for auto-detection or explicitly set thread count based on your system.

   Background Level Determination
   The mode (calculated using half-sample algorithm) is used for setting
   background levels between images. This is more robust than mean or median
   for skewed distributions common in astronomical images.

   Notes
   -----

   This routine can take considerable time because:

   * Each tile involves many files (often 100-1000)
   * Must check all pairwise combinations (N × (N-1) matches)
   * For 100 files: ~10,000 comparisons
   * For 1000 files: ~1,000,000 comparisons

   **Optimization:**

   Running with maximum available threads is strongly recommended. Use ``-np 0``
   for auto-detection or explicitly set thread count based on your system.

   Background Level Determination
   The mode (calculated using half-sample algorithm) is used for setting
   background levels between images. This is more robust than mean or median
   for skewed distributions common in astronomical images.

   Version History
   ---------------

   230505 ksl

       Coding begun and routine parallelized

   230617 ksl

       Adapted to new version of routines

   230704 ksl

       Mode and additional statistics added. Mode now used for background levels.

   Author
   Space Telescope Science Institute



Attributes
----------

.. autoapisummary::

   BackStats.BACKDIR
   BackStats.ierror


Functions
---------

.. autoapisummary::

   BackStats.calculate_one
   BackStats.do_one_tile
   BackStats.halfsamplemode
   BackStats.quartile_stats
   BackStats.steer


Module Contents
---------------

.. py:data:: BACKDIR
   :value: 'DECam_BACK'


.. py:function:: calculate_one(file, xmatch_files, indir, calc_sigma_clipped=False, npix_min=100)

   Calculate background statistics for one file vs multiple matches.

   Computes statistics on the difference in overlap regions between one
   base file and a list of comparison files.

   Parameters
   ----------
   file : str
       Base FITS filename.
   xmatch_files : list of str
       List of FITS filenames to compare against base file.
   indir : str
       Directory containing the FITS files.
   calc_sigma_clipped : bool, optional
       If True, calculate sigma-clipped statistics. Increases runtime
       significantly. Default: False.
   npix_min : int, optional
       Minimum number of overlapping pixels required. Pairs with fewer
       pixels are skipped. Default: 100.

   Returns
   -------
   list of list
       List of records, where each record is a list containing:
       
       * file1, file2 : str - Filenames
       * npix : int - Number of overlapping pixels
       * mean, med, std : float - Unbiased statistics
       * mean_clipped, med_clipped, std_clipped : float - Sigma-clipped stats
       * skew, kurt : float - Higher moments
       * Delta : float - Mode of difference

   Notes
   -----
   **Processing Steps:**

   1. Open base file
   2. For each comparison file:
      
      - Find overlapping pixels (both non-zero)
      - Calculate difference in overlap region
      - Compute statistics (mean, median, std)
      - Optionally compute sigma-clipped stats
      - Calculate skewness and kurtosis
      - Calculate mode using half-sample algorithm

   3. Close files and clean up memory

   **Memory Management:**

   Uses explicit garbage collection and file closure to handle large
   datasets without memory leaks.

   **Error Handling:**

   * Sets global ierror flag on failures
   * Logs errors for problematic file pairs
   * Skips pairs with insufficient overlap or NaN results

   Examples
   --------
   >>> records = calculate_one('base.fits', ['comp1.fits', 'comp2.fits'],
   ...                         'DECam_BACK/LMC_c42/T07')
   >>> print(f"Processed {len(records)} comparisons")


.. py:function:: do_one_tile(field='LMC_c42', tile='T07', nproc=1, calc_sigma_clipped=False)

   Process all overlapping file pairs for one tile.

   Orchestrates background statistics calculation for all file pairs in
   a tile, with optional parallel processing.

   Parameters
   ----------
   field : str, optional
       Field name. Default: 'LMC_c42'.
   tile : str, optional
       Tile identifier. Default: 'T07'.
   nproc : int, optional
       Number of parallel processes. If 1, uses serial processing.
       Default: 1.
   calc_sigma_clipped : bool, optional
       If True, calculate sigma-clipped statistics. Default: False.

   Returns
   -------
   Table
       Astropy Table containing background statistics for all file pairs.

   Raises
   ------
   IOError
       If required input files or directories cannot be found.

   Notes
   -----
   **Input Files:**

   * ``Summary/{field}_{tile}_overlap.txt`` - Overlap file pairs
   * ``Summary/{field}_{tile}.txt`` - Image metadata
   * ``DECam_BACK/{field}/{tile}/`` - Background FITS files

   **Processing:**

   1. Read overlap file listing all file pairs
   2. Create work queue for all base files
   3. Process in parallel (if nproc > 1) or serially
   4. Combine results into single table
   5. Rescale by pixel area (0.2631² → 2² arcsec²)
   6. Join with metadata (filter, exposure time)
   7. Write to Summary directory

   **Memory Management:**

   Uses explicit garbage collection and careful table handling to avoid
   memory leaks in astropy join operations.

   **Output Scaling:**

   Results are scaled by (0.2631/2)² to account for pixel size difference
   between native and binned images.

   Examples
   --------
   >>> # Serial processing
   >>> tab = do_one_tile('LMC_c42', 'T07')

   >>> # Parallel processing with 8 threads
   >>> tab = do_one_tile('LMC_c42', 'T07', nproc=8)

   >>> # With sigma-clipped statistics
   >>> tab = do_one_tile('LMC_c42', 'T07', nproc=8, calc_sigma_clipped=True)


.. py:function:: halfsamplemode(inputData, axis=None)

   Compute mode of an array using the half-sample algorithm.

   This is a robust estimator of the mode that works well for
   approximately unimodal distributions.

   Parameters
   ----------
   inputData : ndarray
       Input array of values.
   axis : int, optional
       Axis along which to compute the mode. If None, computes mode
       of flattened array. Default: None.

   Returns
   -------
   float
       Mode estimate.

   Notes
   -----
   **Algorithm:**

   1. Sort the data
   2. Find the narrowest half-sample (shortest interval containing
      half the data)
   3. Recursively apply to the narrowest half
   4. Continue until 1-2 values remain

   This algorithm is more robust than simple binning approaches and
   works well for skewed distributions common in background differences.

   **Computational Complexity:**

   O(n log n) due to initial sort, then O(log n) iterations.

   Examples
   --------
   >>> data = np.random.normal(5, 1, 1000)
   >>> mode = halfsamplemode(data)
   >>> print(f"Mode: {mode:.2f}")


.. py:data:: ierror
   :value: False


.. py:function:: quartile_stats(xdata)

   Compute skewness and kurtosis with outlier rejection.

   Uses quartile-based outlier rejection before computing higher-order
   moments, providing robust statistics for distributions with outliers.

   Parameters
   ----------
   xdata : ndarray
       Input data array.

   Returns
   -------
   skew : float
       Skewness (third standardized moment).
   kurt : float
       Kurtosis (fourth standardized moment).

   Notes
   -----
   **Outlier Rejection:**

   Uses the IQR (interquartile range) method:

   * Q1 = 25th percentile
   * Q3 = 75th percentile
   * IQR = Q3 - Q1
   * Outliers: values < Q1 - 1.5×IQR or > Q3 + 1.5×IQR

   After outlier removal, computes standard moments.

   **Reference:**

   See: https://towardsdatascience.com/skewness-and-kurtosis-with-outliers-f43167532c69

   **Interpretation:**

   * Skewness = 0: symmetric distribution
   * Skewness > 0: right-skewed (tail extends right)
   * Skewness < 0: left-skewed (tail extends left)
   * Kurtosis = 3: normal distribution
   * Kurtosis > 3: heavy tails
   * Kurtosis < 3: light tails

   Examples
   --------
   >>> data = np.random.normal(0, 1, 10000)
   >>> skew, kurt = quartile_stats(data)
   >>> print(f"Skewness: {skew:.3f}, Kurtosis: {kurt:.3f}")


.. py:function:: steer(argv)

   Parse command-line arguments and execute background statistics.

   Main entry point for command-line execution. Handles argument parsing
   and orchestrates tile processing.

   Parameters
   ----------
   argv : list
       Command-line arguments (typically sys.argv).

   Returns
   -------
   None
       Results written to files and logged.

   Notes
   -----
   **Processing:**

   1. Parse arguments
   2. Expand -all to tiles T01-T16 if specified
   3. Open log file
   4. Process each tile with timing
   5. Optionally remove BackPrep files on success
   6. Close log

   **Error Handling:**

   If global ierror flag is set, does not remove BackPrep files even
   with -rm option.

   Examples
   --------
   From command line::

       python BackStats.py LMC_c42 T07 T08
       python BackStats.py -all -np 8 LMC_c42
       python BackStats.py -rm -sigma_clipped LMC_c42 T07


