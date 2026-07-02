PsfPhot
=======

.. py:module:: PsfPhot

.. autoapi-nested-parse::

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       ::

       Usage: PsfPhot.py [-out whatever] [-np N] image.fits psf.fits stars.fits

       where:
       -out : Output filename root (default: 'psf_phot')
       -np  : Number of parallel processes (default: 1, 0=auto-detect)

   Notes
   -----

   This parallelized version handles large images efficiently and includes:

   * Saturation masking to prevent flux overestimation
   * Post-fit peak consistency checks
   * Robust source grouping with size limits
   * Parallel processing with joblib
   * Progress monitoring with heartbeat status messages
   * Timing information for each processing chunk

   The algorithm performs post-processing validation by comparing predicted
   peak values with actual image values in each region.



Attributes
----------

.. autoapisummary::

   PsfPhot.BACKGROUND_CONFIG
   PsfPhot.HAS_JOBLIB
   PsfPhot.PSF_CONFIG
   PsfPhot.READNOISE
   PsfPhot.SATURATION_LIMIT


Classes
-------

.. autoapisummary::

   PsfPhot.HeartbeatPrinter
   PsfPhot.RobustSourceGrouper


Functions
---------

.. autoapisummary::

   PsfPhot.create_psf_photometry_objects
   PsfPhot.do_crowded
   PsfPhot.process_single_chunk
   PsfPhot.steer


Module Contents
---------------

.. py:data:: BACKGROUND_CONFIG

.. py:data:: HAS_JOBLIB
   :value: True


.. py:class:: HeartbeatPrinter(chunk_id, interval=900)

   Bases: :py:obj:`threading.Thread`


   Background thread that prints periodic status messages.

   Provides visual confirmation that long-running processes haven't stalled
   by printing timestamped heartbeat messages at regular intervals.

   Parameters
   ----------
   chunk_id : int
       Identifier for the processing chunk being monitored.
   interval : int, optional
       Time between heartbeat messages in seconds (default: 900 = 15 minutes).

   Attributes
   ----------
   chunk_id : int
       The chunk identifier.
   interval : int
       Heartbeat interval in seconds.
   stop_event : threading.Event
       Event used to signal thread termination.

   Examples
   --------
   >>> heartbeat = HeartbeatPrinter(chunk_id=1, interval=600)
   >>> heartbeat.start()
   >>> # ... do work ...
   >>> heartbeat.stop()


   .. py:attribute:: chunk_id


   .. py:attribute:: daemon
      :value: True


      A boolean value indicating whether this thread is a daemon thread.

      This must be set before start() is called, otherwise RuntimeError is
      raised. Its initial value is inherited from the creating thread; the
      main thread is not a daemon thread and therefore all threads created in
      the main thread default to daemon = False.

      The entire Python program exits when only daemon threads are left.




   .. py:attribute:: interval
      :value: 900



   .. py:method:: run()

      Execute the heartbeat loop until stopped.



   .. py:method:: stop()

      Stop the heartbeat thread.



   .. py:attribute:: stop_event


.. py:data:: PSF_CONFIG

.. py:data:: READNOISE
   :value: 5.0


.. py:class:: RobustSourceGrouper(min_separation, max_group_size=50)

   Bases: :py:obj:`photutils.psf.SourceGrouper`


   Enhanced SourceGrouper with maximum group size enforcement.

   Extends photutils SourceGrouper to prevent exponential slowdown
   in dense stellar clusters by splitting oversized groups.

   Parameters
   ----------
   min_separation : float
       Minimum separation distance for grouping sources (pixels).
   max_group_size : int, optional
       Maximum number of sources allowed in a single group (default: 50).
       Groups exceeding this size are split into smaller subgroups.

   Attributes
   ----------
   max_group_size : int
       Maximum allowed group size.

   Examples
   --------
   >>> grouper = RobustSourceGrouper(min_separation=5.0, max_group_size=30)
   >>> groups = grouper(x_positions, y_positions)


   .. py:method:: __call__(x, y)

      Assign group IDs to sources, splitting oversized groups.

      Parameters
      ----------
      x : array_like
          X coordinates of sources.
      y : array_like
          Y coordinates of sources.

      Returns
      -------
      ndarray
          Group ID for each source.



   .. py:attribute:: max_group_size
      :value: 50



.. py:data:: SATURATION_LIMIT
   :value: 3000


.. py:function:: create_psf_photometry_objects(image_data, psf_data, global_bkg_median, readnoise=None)

   Create and configure PSF photometry objects.

   Sets up the PSF model, NDData structure, background estimator, and
   source grouper needed for PSF photometry.

   Parameters
   ----------
   image_data : ndarray
       2D image data array.
   psf_data : ndarray
       2D PSF model array (should be normalized).
   global_bkg_median : float
       Median background level for Poisson uncertainty estimation.
   readnoise : float, optional
       Detector read noise in electrons/ADU (default: uses READNOISE config).

   Returns
   -------
   psf_phot : PSFPhotometry
       Configured PSF photometry object.
   nddata : NDData
       Image data wrapped with uncertainty information.

   Notes
   -----
   The function configures:

   * ImagePSF model with specified oversampling
   * Local background estimation using MMM algorithm
   * Robust source grouping with separation constraints
   * Levenberg-Marquardt least-squares fitter

   **Uncertainty Calculation**: Uses Poisson statistics based on the data
   values plus read noise: uncertainty = sqrt(data + readnoise²). This is
   appropriate when local background subtraction is performed, as it avoids
   double-counting background variations that would lead to unrealistically
   low reduced chi-squared values.


.. py:function:: do_crowded(image, psf, source_table, outroot='psf_phot', n_processes=1, readnoise=None)

   Perform PSF photometry on a crowded field image.

   Main processing function that orchestrates parallel PSF photometry
   on large source catalogs, with chunked processing and optional
   multi-core parallelization.

   Parameters
   ----------
   image : str
       Path to FITS image file.
   psf : str
       Path to FITS PSF model file.
   source_table : str
       Path to FITS table containing source positions.
       Must have 'xcentroid' and 'ycentroid' columns.
   outroot : str, optional
       Output filename root (default: 'psf_phot').
       Results saved as ``{outroot}.fits``.
   n_processes : int, optional
       Number of parallel processes (default: 1).
       Set to 0 for auto-detection (uses CPU count - 1).
       Requires joblib for n_processes > 1.
   readnoise : float, optional
       Detector read noise in electrons/ADU (default: uses READNOISE config).
       Used for uncertainty calculation.

   Returns
   -------
   None
       Results are written to disk.

   Notes
   -----
   **Processing Steps**:

   1. Load image, PSF, and source catalog
   2. Estimate global background median (for uncertainty calculation)
   3. Split sources into chunks of 5000
   4. Process chunks in parallel (if enabled)
   5. Validate results with peak consistency checks
   6. Write combined results to FITS table

   **Uncertainty Calculation**: Uses Poisson statistics (sqrt(counts)) plus
   read noise, rather than background RMS. This is appropriate because 
   PSFPhotometry performs local background subtraction, which removes the
   spatial background variations that would otherwise be included in the RMS.
   Using background RMS would lead to unrealistically low reduced chi-squared
   values due to overestimated uncertainties.

   **Output Columns**:

   Standard PSFPhotometry columns plus:

   * model_peak : Model-predicted peak value
   * data_peak : Measured peak value
   * peak_ratio : Quality metric (model/data)

   **Performance**: Chunk size of 5000 balances memory usage
   and parallelization efficiency. Smaller chunks increase overhead;
   larger chunks increase memory footprint.

   Examples
   --------
   >>> # Single-threaded processing
   >>> do_crowded('image.fits', 'psf.fits', 'sources.fits')

   >>> # Parallel processing with 4 cores
   >>> do_crowded('image.fits', 'psf.fits', 'sources.fits', 
   ...            outroot='my_phot', n_processes=4)

   >>> # Auto-detect CPU count with custom read noise
   >>> do_crowded('image.fits', 'psf.fits', 'sources.fits', 
   ...            n_processes=0, readnoise=7.0)


.. py:function:: process_single_chunk(chunk_info, image_data, psf_data, sources_x, sources_y, global_bkg_median, readnoise=None)

   Process a single chunk of sources with PSF photometry.

   Performs PSF fitting on a subset of sources, including saturation masking
   and post-fit validation of peak consistency. Reports timing information.

   Parameters
   ----------
   chunk_info : tuple
       Tuple of (chunk_index, start_idx, end_idx, n_chunks) defining the chunk.
   image_data : ndarray
       2D image data array.
   psf_data : ndarray
       Normalized 2D PSF model array.
   sources_x : ndarray
       X coordinates of all sources.
   sources_y : ndarray
       Y coordinates of all sources.
   global_bkg_median : float
       Median background level for uncertainty calculation.
   readnoise : float, optional
       Detector read noise in electrons/ADU (default: uses READNOISE config).

   Returns
   -------
   Table or None
       Photometry results table with additional columns:
       
       * model_peak : Predicted peak pixel value from fitted flux
       * data_peak : Actual peak pixel value in image
       * peak_ratio : Ratio of model_peak to data_peak
       
       Returns None if processing fails.

   Notes
   -----
   **Timing Information**: The function prints:

   * Start time and source index range
   * End time and elapsed duration
   * Timing breakdown: initialization, fitting, and validation

   **Saturation Handling**: Pixels above SATURATION_LIMIT are masked
   during fitting to prevent flux overestimation from saturated cores.

   **Peak Validation**: Post-fit check compares the model-predicted peak
   with the actual data peak to identify problematic fits.


.. py:function:: steer(argv)

   Parse command-line arguments and execute photometry.

   Parameters
   ----------
   argv : list
       Command-line argument list (typically sys.argv).

   Command-Line Options
   --------------------
   -out <root>
       Output filename root (default: 'psf_phot')
   -np <N>
       Number of parallel processes (default: 1, 0=auto)

   Positional Arguments
   --------------------
   image.fits
       Input image FITS file
   psf.fits
       PSF model FITS file
   stars.fits
       Source catalog FITS table

   Examples
   --------
   Command line usage::

       $ python PsfPhot.py image.fits psf.fits stars.fits
       $ python PsfPhot.py -out myresults -np 4 image.fits psf.fits stars.fits
       $ python PsfPhot.py -np 0 image.fits psf.fits stars.fits  # auto-detect CPUs


