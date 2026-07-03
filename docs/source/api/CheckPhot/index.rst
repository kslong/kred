CheckPhot
=========

.. py:module:: CheckPhot

.. autoapi-nested-parse::

   CheckPhot - Photometric consistency check for a DECam field

   Space Telescope Science Institute

   Synopsis
   --------

   Runs the complete photometric consistency workflow for one field and
   produces two concise summary tables: one for intra-filter consistency
   (are all exposures of the same filter on the same flux scale?) and one
   for inter-filter consistency (does the same star have the same DN in the
   r-band image as in the emission-line images?).

   Command Line Usage
   ------------------

   ::

       CheckPhot.py LMC_c42
       CheckPhot.py -np 16 LMC_c42
       CheckPhot.py -no_mefphot LMC_c42        # skip MefPhot if TabPhot/ already populated
       CheckPhot.py -redo LMC_c42              # force reprocessing of all TabPhot files
       CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits   # reprint summary only
       CheckPhot.py -config myconfig.txt LMC_c42   # alternate image-group config

   Flags
   -----

   ``-np N``
       Number of parallel MefPhot workers (default 8).

   ``-no_mefphot``
       Skip the MefPhot step; use existing ``TabPhot/`` files.

   ``-redo``
       Reprocess already-existing CalcZeroPoint/PhotEval outputs.

   ``-o DIR``
       Write all outputs to DIR (default: ``CheckPhot``).

   ``-config FILE``
       Alternate image-group config file (default: ``DeMCELS_images.txt``
       in current directory, then ``$KRED/config/``).

   ``-summary FITS``
       Read an existing ``*_phot_check.fits`` and reprint the formatted
       summary without running any photometry.

   ``-ref COL``
       Additional ZP reference column for CalcZeroPoint (may be repeated).

   ``-h``
       Print this help.

   Description
   -----------

   Exposures are grouped by *Image label* read from ``config/DeMCELS_images.txt``
   (e.g. ``N662`` = 400-800 s N662 images, ``N662_s`` = 30-60 s N662 images).
   Files in ``TabPhot/`` are assigned to groups by matching the ``FILTER``
   header keyword and ``EXPTIME`` (within 5 s).  Four steps then run in sequence:

   1. **MefPhot** — forced aperture photometry at SMASH and Gaia positions.
      Already-processed files are skipped unless ``-redo`` is given.

   2. **CalcZeroPoint** — per-file zero-point fit vs catalog for each Image group.
      Measures how far each file's header MAGZERO deviates from the measured ZP.

   3. **PhotEval** — groups detections of the same star across overlapping frames
      and computes ``chi2_nu_c``: the reduced chi-squared of the MAGZERO-corrected
      magnitudes relative to photon noise.

   4. **PhotAnal** — cross-matches the same star between two Image groups and
      measures the scatter in ``(magc_1 - magc_2)``, which sets the minimum
      stellar residual in CleanStars output.  Pairs are formed automatically
      within each tier (long vs short exposures).

   Output
   ------

   All outputs are written to ``CheckPhot/`` (or the directory given via ``-o``).
   All filenames include the field name so multiple fields share the directory
   without collision.

   ``{field}_phot_check.fits`` — single FITS file with two binary table extensions:

   **INTRA** (ext 1) — one row per Image group.  Zero-point comparison columns
   are named ``{stat}_{ref}`` where ``{ref}`` encodes the catalog and reference
   band (e.g. ``smash_r``, ``gaia_g``).  Consistency columns: ``N_sources``,
   ``chi2_nu_c_med``, ``chi2_nu_c_emp_med``, ``sigma_c_magzero``,
   ``sigma_c_gaia_g``, ``sigma_c_smash_r``.

   **INTER** (ext 2) — one row per Image-group pair: ``Image1``, ``Image2``,
   ``Filter1``, ``Filter2``, ``N_stars``, ``mean_delta_mag``, ``sigma_c_magzero``,
   ``sigma_c_gaia_g``, ``sigma_c_smash_r``, ``sigma_residual``, ``residual_pct``,
   ``color_term_b``, ``N_pairs``, ``std_dzp_magzero``, ``std_dzp_gaia_g``,
   ``std_dzp_smash_r``.

   ``sigma_c_magzero`` and ``sigma_c_gaia_g`` / ``sigma_c_smash_r`` name the
   same concept in both extensions — per-star scatter after applying a given
   ZP reference — so the INTRA and INTER tables can be compared directly.

   Intermediate files in ``CheckPhot/``:

   - ``{field}_zeropoints_{label}.fits``         — CalcZeroPoint per-file ZP table
   - ``{field}_phot_eval_{label}.fits``          — PhotEval per-source chi2 table
   - ``{field}_filter_compare_{l1}_{l2}.fits``   — PhotAnal matched-star table

   Interpreting the tables
   -----------------------

   Intra-filter:
     ``chi2_nu_c_emp_med ≈ 1`` means frames are mutually consistent at the
     photon-noise level.  ``std_dzp_{ref}`` is the MAD scatter of per-file
     MAGZERO deviations (robust to outliers); > 0.06 mag suggests the
     empirical ZP correction may help.  ``sigma_c_magzero`` vs
     ``sigma_c_smash_r`` shows whether MAGZERO or the SMASH calibration
     gives better intra-group consistency.

   Inter-filter:
     ``residual_pct`` is the minimum percentage of stellar flux that will
     remain as a residual after CleanStars continuum subtraction.
     ``sigma_residual`` is the scatter after removing the color-term trend;
     this is the irreducible floor set by stellar color diversity.
     ``std_dzp_magzero`` measures image-to-image stability of the
     inter-filter flux ratio (the inter-filter analogue of ``std_dzp`` in
     the INTRA table).

   History
   -------

   260701 ksl  Written
   260702 ksl  Split ZP stats by (catalog, ref_col); add -ref flag
   260703 ksl  Image-group approach from DeMCELS_images.txt; -summary flag; column renames



Functions
---------

.. autoapisummary::

   CheckPhot.check_phot
   CheckPhot.steer
   CheckPhot.summarize


Module Contents
---------------

.. py:function:: check_phot(field, mef_dir='DECam_MEF', tabphot_dir='TabPhot', config_file=None, np_proc=8, redo=False, run_mefphot=True, checkdir=None, extra_refs=None)

   Run the full photometric consistency check for *field*.

   Parameters
   ----------
   field : str
       Field name, e.g. ``'LMC_c42'``.
   mef_dir : str
       Directory containing ``{field}/mef/*.fits.fz`` (default ``DECam_MEF``).
   tabphot_dir : str
       Directory where MefPhot writes its output (default ``TabPhot``).
   config_file : str or None
       Path to DeMCELS_images.txt.  Searched in current directory then
       $KRED/config/ if not given.
   np_proc : int
       MefPhot parallel workers (default 8).
   redo : bool
       Reprocess existing TabPhot files (default False).
   run_mefphot : bool
       If False, skip MefPhot and assume TabPhot/ is already populated.
   checkdir : str or None
       Directory for all CheckPhot outputs (intermediate and summary).
       Default: ``CheckPhot``.
   extra_refs : list of str or None
       Additional catalog column names to compare against in CalcZeroPoint.

   Returns
   -------
   intra : astropy.table.Table
       One row per Image group with intra-filter consistency statistics.
   inter : astropy.table.Table
       One row per Image-group pair with inter-filter consistency statistics.


.. py:function:: steer(argv)

.. py:function:: summarize(fits_file)

   Read a ``*_phot_check.fits`` file and print the comparison summary.

   Can be called as::

       CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits


