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
     calibrated magnitude using the per-image zero point measured by MagZero.
     When images have different zero points these corrections matter; the
     corrected scatter (``magc_std``) is the meaningful repeatability metric.

   The key diagnostic is the scatter (``magc_std``) relative to the expected
   photon-noise floor (``mag_err_mean``) and the reduced chi^2 (``chi2_nu_c``).

   Optional Arguments
   ------------------

   -h
       Print this help and exit.

   -dir DIR
       Directory containing the ``*.smash.fits`` catalogs (default: ``TabPhot``).

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

   -min_n MIN_N
       Minimum number of qualifying detections required to include a source
       in the output (default: 3).

   -o OUTPUT
       Output FITS filename (default: ``phot_eval_<FILTER>_<EXPTIME>.fits``
       in the current directory).

   -zp_table FILE
       FITS table of empirical zero points (e.g. ``zeropoints.fits`` produced
       by CalcZeroPoint).  Must contain ``Filename`` and ``zp_calc`` columns.
       When supplied, a second set of corrected magnitudes is computed using
       ``zp_calc`` in place of the header ``MAGZERO``, giving additional output
       columns ``magc_emp_*`` and ``chi2_nu_emp`` for direct comparison.

   Output Columns
   --------------

   Source_name, RA, Dec, n_detect,
   mag_mean, mag_wmean, mag_median, mag_std,           (raw, ZP = 28)
   magc_mean, magc_wmean, magc_median, magc_std,       (corrected, ZP = MAGZERO)
   magc_emp_mean, magc_emp_wmean, magc_emp_median, magc_emp_std,  (empirical ZP, if -zp_table given)
   mag_err_mean, snr_mean, chi2_nu, chi2_nu_c, chi2_nu_emp

   Examples
   --------

   Evaluate r-band 30-second exposures::

       PhotEval.py -filter r -exp 30 -snr 10 -min_n 5

   Evaluate Ha (N662) band, all exposure times::

       PhotEval.py -filter N662 -min_n 3 -o ha_eval.fits

   Compare header MAGZERO vs empirical ZP from CalcZeroPoint::

       PhotEval.py -filter r -zp_table zeropoints.fits -o phot_eval_r_emp.fits



Functions
---------

.. autoapisummary::

   PhotEval.do_eval
   PhotEval.steer


Module Contents
---------------

.. py:function:: do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile, zp_lookup=None)

.. py:function:: steer(argv)

