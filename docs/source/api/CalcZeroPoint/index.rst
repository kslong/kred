CalcZeroPoint
=============

.. py:module:: CalcZeroPoint

.. autoapi-nested-parse::

   CalcZeroPoint - Calculate photometric zero points from TabPhot catalogs

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       CalcZeroPoint.py [-h] [-dir DIR] [-filter FILTER] [-snr MIN_SNR]
                        [-ref COL] [-sigma N] [-mag_lo LO] [-mag_hi HI]
                        [-o OUTPUT]

   For each ``*.smash.fits`` or ``*.gaia.fits`` catalog in ``DIR`` (default:
   ``TabPhot/``), computes a best-fit photometric zero point by comparing
   instrumental magnitudes (``phot_mag``, which use a fixed ZP = 28) to
   reference catalog magnitudes (``R`` for SMASH, ``G`` for Gaia).

   The zero point estimate for a single star is::

       zp_i = 28 + (ref_mag_i - phot_mag_i)

   After quality cuts, the distribution of ``zp_i`` values is sigma-clipped
   and summarised.  The output table has one row per input file and is
   designed to be compared directly against the ``MAGZERO`` keyword carried
   in each MEF file header.

   Running with no arguments processes all catalog files found in ``TabPhot/``
   using default settings.  Use ``-h`` to print this help.

   Optional Arguments
   ------------------

   -h
       Print this help and exit.

   -dir DIR
       Directory containing the ``*.smash.fits`` / ``*.gaia.fits`` catalogs
       (default: ``TabPhot``).  Run from the working directory containing
       ``TabPhot/``.

   -filter FILTER
       Only process files whose ``Filter`` column value starts with FILTER
       (e.g. ``r``, ``N662``, ``N673``).  The filter name is read from
       inside each file, not from the filename.
       Default: process all catalog files found.

   -snr MIN_SNR
       Minimum SNR for a detection to be used (default: 20).

   -ref COL
       Reference magnitude column to use (overrides the per-file default of
       ``R`` for SMASH and ``G`` for Gaia).

   -sigma N
       Number of sigma for iterative sigma-clipping (default: 3).

   -mag_lo LO
       Faint magnitude limit for reference stars (default: 21).
       Excludes faint stars where photon noise dominates.

   -mag_hi HI
       Bright magnitude limit for reference stars (default: 15).
       Excludes stars that may be saturated.

   -o OUTPUT
       Output FITS filename (default: ``zeropoints.fits`` in the current
       directory).

   Output Columns
   --------------

   +------------+------------------------------------------------------------+
   | Column     | Description                                                |
   +============+============================================================+
   | Filename   | Input TabPhot catalog filename (basename)                  |
   +------------+------------------------------------------------------------+
   | Filter     | Filter name as stored in the catalog                       |
   +------------+------------------------------------------------------------+
   | Exptime    | Exposure time in seconds                                   |
   +------------+------------------------------------------------------------+
   | Catalog    | Reference catalog used (``SMASH`` or ``Gaia``)             |
   +------------+------------------------------------------------------------+
   | ref_col    | Reference magnitude column used (``R`` or ``G``)           |
   +------------+------------------------------------------------------------+
   | MAGZERO    | Zero point from the MEF image header                       |
   +------------+------------------------------------------------------------+
   | zp_calc    | Our derived zero point (sigma-clipped median)              |
   +------------+------------------------------------------------------------+
   | zp_wmean   | Inverse-variance-weighted mean zero point                  |
   +------------+------------------------------------------------------------+
   | zp_std     | Scatter of individual per-star ZP estimates after clipping |
   +------------+------------------------------------------------------------+
   | zp_err     | Standard error on ``zp_calc`` (= zp_std / sqrt(n_stars))  |
   +------------+------------------------------------------------------------+
   | zp_mad     | Median absolute deviation (robust scatter measure)         |
   +------------+------------------------------------------------------------+
   | n_stars    | Stars used after sigma-clipping                            |
   +------------+------------------------------------------------------------+
   | n_total    | Stars passing initial quality cuts (before clipping)       |
   +------------+------------------------------------------------------------+
   | delta_zp   | ``zp_calc − MAGZERO`` (offset from pipeline value)         |
   +------------+------------------------------------------------------------+

   The ``delta_zp`` column is the primary diagnostic: values consistently
   offset from zero across many files indicate a systematic difference between
   the reference catalog photometric system and the pipeline zero point.

   Relationship to Other Tools
   ---------------------------

   * **Input**: ``TabPhot/*.smash.fits`` or ``TabPhot/*.gaia.fits`` files
     produced by :doc:`MefPhot </api/MefPhot/index>`.
   * **ZeroCalc**: a complementary zero-point tool that uses a weighted linear
     regression (with optional color term) instead of a sigma-clipped median.
     ``ZeroCalc`` is preferred when a color term is scientifically important;
     ``CalcZeroPoint`` is preferred for a quick multi-file summary.  Both
     tools read the same ``TabPhot/`` files and can be run independently.
   * **PhotEval**: uses the ``MAGZERO`` values from the MEF headers (not
     the derived ``zp_calc`` from this script) to apply zero-point corrections
     in its per-source scatter analysis.

   Examples
   --------

   Process all catalog files with default settings::

       CalcZeroPoint.py

   Process r-band files only::

       CalcZeroPoint.py -filter r -snr 20

   Process all filters with a brighter magnitude limit::

       CalcZeroPoint.py -mag_hi 14 -mag_lo 19 -o zp_bright.fits

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

   CalcZeroPoint.steer
   CalcZeroPoint.do_one


Module Contents
---------------

.. py:function:: steer(argv)

   Parse command-line arguments, find catalog files, and call
   :func:`do_one` for each.

.. py:function:: do_one(filepath, snr_min, ref_col_override, n_sigma, mag_hi, mag_lo)

   Compute the zero point for a single TabPhot catalog file.

   Parameters
   ----------
   filepath : str
       Path to a ``*.smash.fits`` or ``*.gaia.fits`` catalog file.
   snr_min : float
       Minimum SNR cut applied before computing ZP estimates.
   ref_col_override : str or None
       If given, use this column as the reference magnitude; otherwise
       ``'R'`` for SMASH and ``'G'`` for Gaia.
   n_sigma : float
       Sigma threshold for iterative sigma-clipping.
   mag_hi : float
       Bright magnitude limit (exclude brighter stars).
   mag_lo : float
       Faint magnitude limit (exclude fainter stars).

   Returns
   -------
   dict or None
       Dictionary of output columns (one row of the summary table), or
       ``None`` if the file cannot be read or has too few usable stars.
