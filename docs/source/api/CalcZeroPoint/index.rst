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
   ``TabPhot``), computes a best-fit photometric zero point by comparing
   instrumental magnitudes (``phot_mag``, which use a fixed ZP = 28) to
   reference catalog magnitudes (``R`` for SMASH, ``G`` for Gaia).

   The zero point for a single star is::

       zp_i = 28 + (ref_mag_i - phot_mag_i)

   After quality cuts, the distribution of ``zp_i`` is sigma-clipped and
   summarised.  The output table has one row per input file.

   Optional Arguments
   ------------------

   -h
       Print this help and exit.

   -dir DIR
       Directory containing the ``*.smash.fits`` / ``*.gaia.fits`` catalogs
       (default: ``TabPhot``).

   -filter FILTER
       Only process files whose Filter header value starts with FILTER
       (e.g. ``r``, ``N662``, ``N673``).  The filter name is read from
       the ``Filter`` column inside each file, not from the filename.
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
       Output FITS filename (default: ``zeropoints.fits``).

   Output Columns
   --------------

   Filename, Filter, Exptime, Catalog, ref_col, MAGZERO,
   zp_calc, zp_wmean, zp_std, zp_err, zp_mad,
   n_stars, n_total, delta_zp

   where:

   * ``zp_calc``  -- sigma-clipped median of individual ZP estimates
   * ``zp_wmean`` -- inverse-variance-weighted mean ZP
   * ``zp_std``   -- std of clipped ZP estimates (scatter / repeatability)
   * ``zp_err``   -- standard error on zp_calc (zp_std / sqrt(n_stars))
   * ``zp_mad``   -- median absolute deviation of clipped estimates
   * ``n_stars``  -- number of stars used after sigma-clipping
   * ``n_total``  -- number of stars passing initial quality cuts
   * ``delta_zp`` -- zp_calc - MAGZERO (should be ~0 if header ZP is good)

   Examples
   --------

   Process all SMASH r-band files::

       CalcZeroPoint.py -filter r -snr 20

   Process all filters, bright stars only::

       CalcZeroPoint.py -mag_hi 13 -mag_lo 18 -o zp_bright.fits



Functions
---------

.. autoapisummary::

   CalcZeroPoint.do_one
   CalcZeroPoint.steer


Module Contents
---------------

.. py:function:: do_one(filepath, snr_min, ref_col_override, n_sigma, mag_hi, mag_lo)

.. py:function:: steer(argv)

