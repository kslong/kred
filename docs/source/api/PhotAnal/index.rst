PhotAnal
========

.. py:module:: PhotAnal

.. autoapi-nested-parse::

   Compare stellar photometry between two filters to assess inter-filter flux consistency

   Space Telescope Science Institute

   Synopsis
   --------

   PhotAnal measures the same stars in two different filters and reports how
   consistently they are scaled relative to each other.  This is the key
   diagnostic for CleanStars continuum subtraction: if a star has DN_r in the
   r-band image and DN_N662 in the Hα image, CleanStars computes

       ha_clean = N662 - r_pure

   and a star cancels only if DN_N662 = DN_r.  After MefPrep has placed all
   images on the mag-28 scale, this condition is equivalent to asking whether
   the star has the same magc (= phot_mag + MAGZERO − 28) in both filters.
   PhotAnal measures that difference and its scatter.

   Command Line Usage
   ------------------

   ::

       PhotAnal.py -filter1 r -filter2 N662 TabPhot/*.gaia.fits
       PhotAnal.py -filter1 r -filter2 N708 TabPhot/*.gaia.fits
       PhotAnal.py -filter1 r -filter2 N662 -o result.fits -no_plot TabPhot/*.gaia.fits

   Flags
   -----

   -filter1 NAME   Reference filter (e.g. r)
   -filter2 NAME   Comparison filter (e.g. N662, N673, N708)
   -o FILE         Output FITS table (default: filter_compare_<f1>_<f2>.fits)
   -no_plot        Skip the diagnostic plot
   -h              Print this help

   Description
   -----------

   Both filter sets must be Gaia-matched TabPhot files (.gaia.fits).  Run::

       MefPhot.py -r 6 DECam_MEF/LMC_c42/*.fits

   for all filters before running PhotAnal.

   For each star detected in both filter1 and filter2:

       delta_mag = magc_filter1 − magc_filter2
                 = (phot_mag_1 + MAGZERO_1 − 28) − (phot_mag_2 + MAGZERO_2 − 28)

   This is the quantity that determines the stellar residual after continuum
   subtraction.  A non-zero mean indicates a systematic inter-filter offset
   (correctable by adjusting MAGZERO for one filter).  The scatter around
   the best-fit color correction is the irreducible floor:

       typical stellar residual ≈ σ_residual × ln(10) / 2.5
                                 ≈ σ_residual × 0.92    (in fractional DN units)

   Output
   ------

   Printed summary: N_stars, mean Δmag, total scatter, color term b, residual
   scatter, and typical stellar residual percentage.

   filter_compare_<f1>_<f2>.fits: per-star table with Source_name, magc in
   each filter, delta_mag, Gaia G−R color, delta_mag after color correction.

   filter_compare_<f1>_<f2>.png: Δmag vs G−R color with best-fit line.

   Notes
   -----

   The comparison uses per-source median magc aggregated across all files,
   so it is not sensitive to intra-filter frame-to-frame scatter (that is
   diagnosed separately by PhotEval).

   History
   -------

   251123 ksl  Prototype coded
   260701 ksl  Rewritten as inter-filter consistency tool



Functions
---------

.. autoapisummary::

   PhotAnal.do_compare
   PhotAnal.load_filter_data
   PhotAnal.steer


Module Contents
---------------

.. py:function:: do_compare(filenames, filter1, filter2, outfile=None, do_plot=True, extra_zp_lookups=None, _preloaded=None)

   Cross-match filter1 and filter2 photometry and report inter-filter scatter.

   Parameters
   ----------
   filenames : list of str
       Paths to .gaia.fits TabPhot files (both filters mixed together).
   filter1 : str
       Reference filter name (e.g. 'r').
   filter2 : str
       Comparison filter name (e.g. 'N662').
   outfile : str, optional
       Output FITS table path.  Default: filter_compare_<f1>_<f2>.fits.
   do_plot : bool
       If True, save a PNG diagnostic plot.
   extra_zp_lookups : dict of {name: {basename: zp_calc}}, optional
       Additional ZP sources to evaluate in parallel (e.g. from CheckPhot's
       CalcZeroPoint runs).  For each name the same sigma metrics are
       computed and stored in the output meta as STD_DM_{name},
       STD_RE_{name}, N_PR_{name}, STD_PD_{name}.
   _preloaded : tuple of two (Table, list) pairs, optional
       Pre-loaded (tab, chunks) from load_filter_data() for (filter1, filter2).
       When supplied, load_filter_data() is skipped (avoids redundant I/O when
       a filter appears in multiple pairs).


.. py:function:: load_filter_data(filenames, filter_name, extra_zp_lookups=None)

   Load all .gaia.fits TabPhot files matching filter_name.

   Returns (aggregate, chunks) where aggregate is an Astropy Table with one
   row per star: Source_name, magc (median across files using header MAGZERO),
   G, R (Gaia magnitudes), and optionally magc_{name} for each entry in
   extra_zp_lookups.  chunks is the list of per-file tables (before
   aggregation) needed for per-image-pair statistics.

   Parameters
   ----------
   extra_zp_lookups : dict of {name: {basename: zp_calc}}, optional
       Additional ZP sources (e.g. {'gaia_g': {...}, 'smash_r': {...}}).
       Files not in a lookup get NaN for that ZP source.


.. py:function:: steer(argv)

