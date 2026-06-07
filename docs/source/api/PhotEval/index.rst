PhotEval
========

.. py:module:: PhotEval

.. autoapi-nested-parse::

   PhotEval - Photometric consistency evaluation across overlapping frames

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       PhotEval.py [-h] [-dir DIR] [-filter FILTER] [-exp EXPTIME]
                   [-snr MIN_SNR] [-prob MIN_PROB] [-min_n MIN_N] [-o OUTPUT]

   Reads all ``*<FILTER>*.smash.fits`` catalogs from ``DIR`` (default:
   ``TabPhot/``), filters on exposure time and quality, groups detections by
   ``Source_name``, and writes a summary table with per-source photometric
   statistics.

   Two sets of magnitude statistics are produced:

   * **Raw** (``mag_*``): based on ``phot_mag``, which uses a fixed instrumental
     zero point of 28.
   * **Corrected** (``magc_*``): ``phot_mag + MAGZERO - 28``, i.e. the true
     calibrated magnitude using the per-image zero point from the ``MAGZERO``
     header keyword.  When images have different zero points these corrections
     matter; the corrected scatter (``magc_std``) is the meaningful
     repeatability metric.

   The key diagnostic is the corrected scatter (``magc_std``) relative to the
   expected photon-noise floor (``mag_err_mean``) and the reduced chi-squared
   (``chi2_nu_c``).

   Optional Arguments
   ------------------

   -h
       Print this help and exit.

   -dir DIR
       Directory containing the ``*.smash.fits`` catalogs (default: ``TabPhot``).
       Run from the working directory containing ``TabPhot/``.

   -filter FILTER
       Filter string used to select files, e.g. ``r``, ``N662``, ``N673``.
       Files matching ``*<FILTER>*req.smash.fits`` are read.
       Default: ``r``.

   -exp EXPTIME
       Exposure time to select, in seconds (e.g. ``30``).  If omitted all
       exposure times are used.

   -snr MIN_SNR
       Minimum SNR for a detection to be included (default: 10).

   -prob MIN_PROB
       Minimum SMASH star probability for inclusion (default: 0.5).
       Detections with ``prob >= 99`` are treated as sentinels and excluded.

   -min_n MIN_N
       Minimum number of qualifying detections required to include a source
       in the output (default: 3).

   -o OUTPUT
       Output FITS filename (default: ``phot_eval_<FILTER>_<EXPTIME>.fits``
       in the current directory).

   Output Columns
   --------------

   +-----------------+------------------------------------------------------+
   | Column          | Description                                          |
   +=================+======================================================+
   | Source_name     | SMASH source identifier                              |
   +-----------------+------------------------------------------------------+
   | RA, Dec         | Mean sky position across detections                 |
   +-----------------+------------------------------------------------------+
   | n_detect        | Number of qualifying detections                     |
   +-----------------+------------------------------------------------------+
   | mag_mean        | Simple mean of ``phot_mag`` (ZP = 28)               |
   +-----------------+------------------------------------------------------+
   | mag_wmean       | Inverse-variance-weighted mean (ZP = 28)            |
   +-----------------+------------------------------------------------------+
   | mag_median      | Median of ``phot_mag`` (robust centre)              |
   +-----------------+------------------------------------------------------+
   | mag_std         | Scatter of ``phot_mag`` across frames               |
   +-----------------+------------------------------------------------------+
   | magc_mean       | Simple mean of ZP-corrected magnitudes              |
   +-----------------+------------------------------------------------------+
   | magc_wmean      | Inverse-variance-weighted mean (corrected)          |
   +-----------------+------------------------------------------------------+
   | magc_median     | Median of ZP-corrected magnitudes                   |
   +-----------------+------------------------------------------------------+
   | magc_std        | Scatter of corrected magnitudes (key metric)        |
   +-----------------+------------------------------------------------------+
   | mag_err_mean    | Mean photon-noise magnitude error                   |
   |                 | (= 2.5/ln10 × ErrNet/Net)                           |
   +-----------------+------------------------------------------------------+
   | snr_mean        | Mean signal-to-noise ratio across detections        |
   +-----------------+------------------------------------------------------+
   | chi2_nu         | Reduced chi-squared of raw ``phot_mag`` about       |
   |                 | weighted mean                                        |
   +-----------------+------------------------------------------------------+
   | chi2_nu_c       | Reduced chi-squared of corrected magnitudes         |
   |                 | (should be ~1 for photon-noise-limited data)        |
   +-----------------+------------------------------------------------------+

   Performance
   -----------

   The key algorithmic improvement over a source-by-source loop is to use a
   single vectorised ``pandas.groupby().agg()`` call to aggregate all
   detections at once.  This reduces complexity from O(n × m) (where n is
   the number of unique sources and m the number of detections) to
   O(m log m), enabling ~60 million rows from 40 r-band frames to be
   processed in under a minute.

   Relationship to Other Tools
   ---------------------------

   * **Input**: ``TabPhot/*.smash.fits`` files produced by :doc:`MefPhot
     </api/MefPhot/index>`.  Each file contains forced photometry for one
     MEF image; ``PhotEval`` stacks them all and groups by ``Source_name``.
   * **Zero-point correction**: The ``MAGZERO`` column written by ``MefPhot``
     (from the MEF header) is used to compute corrected magnitudes.  For an
     independent check of those zero points see :doc:`CalcZeroPoint
     </api/CalcZeroPoint/index>` (quick summary) or :doc:`ZeroCalc
     </api/ZeroCalc/index>` (regression fit with optional color term).
   * **PhotAnal**: provides a related per-filter consistency check across a
     smaller set of files; ``PhotEval`` is preferred when many overlapping
     frames are available.

   Examples
   --------

   Evaluate r-band 30-second exposures::

       PhotEval.py -filter r -exp 30 -snr 10 -min_n 5

   Evaluate Hα (N662) band, all exposure times::

       PhotEval.py -filter N662 -min_n 3 -o ha_eval.fits

   Version History
   ---------------

   2026-06-07 ksl
       Initial coding.

   Author
   ------
   Space Telescope Science Institute

   .. moduleauthor:: KSL


Functions
---------

.. autoapisummary::

   PhotEval.steer
   PhotEval.do_eval


Module Contents
---------------

.. py:function:: steer(argv)

   Parse command-line arguments and invoke :func:`do_eval`.

.. py:function:: do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile)

   Read, filter, and group MefPhot catalogs; write per-source statistics.

   Parameters
   ----------
   tabphot_dir : str
       Directory containing ``*.smash.fits`` catalogs.
   filter_str : str
       Filter substring used to select files (e.g. ``'r'``, ``'N662'``).
   exptime : float or None
       Exposure time to select (seconds), or ``None`` to use all.
   snr_min : float
       Minimum SNR threshold.
   prob_min : float
       Minimum SMASH star probability.
   min_n : int
       Minimum number of detections to include a source.
   outfile : str
       Output FITS filename.
