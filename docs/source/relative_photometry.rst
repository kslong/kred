====================
Relative Photometry
====================

Relative photometry ensures that all images in a field are on the same
flux scale, so that the same star has the same measured brightness regardless
of which exposure it was taken from.  This matters in two contexts:

* **Mosaicking**: SWarp co-adds overlapping CCD images into tile mosaics.
  If two exposures are not on the same flux scale, the mosaic will have
  visible brightness steps at the seams between them.

* **Continuum subtraction**: ``CleanStars`` builds a continuum image from
  broadband (r) and off-line (N708) data and subtracts it from each
  emission-line image.  The stellar flux must be the same in both images
  so that stars cancel cleanly.  A 5% error in the relative scale between
  filters leaves 5% stellar residuals in the subtracted image.


The Standard Approach: MAGZERO
--------------------------------

The NOIRLAB community pipeline assigns a photometric zero point (``MAGZERO``)
to each MEF file.  ``MefPrep`` uses this to rescale every CCD image:

.. math::

   \text{DN}_{\rm prep} = \text{DN}_{\rm raw} \times 10^{\,0.4\,(28 - \mathrm{MAGZERO})}

After this step, 1 DN in every output image corresponds to a source of
magnitude 28, regardless of the original exposure time, filter, or night.
All images can then be directly compared or co-added.

In practice the NOIRLAB ``MAGZERO`` values are consistent across exposures
to ~10–15 mmag, which is sufficient for most science.  The tools below
allow you to verify this for your specific dataset and, if necessary,
improve it.


Checking Consistency
---------------------

The diagnostic workflow uses ``MefPhot``, ``CalcZeroPoint``, and ``PhotEval``.
All three work on the raw MEF files — it is not necessary to re-run photometry
on the ``DECam_PREP/`` output.  The reason is that the magnitude you would
measure in a PREP image,

.. math::

   m_{\rm prep} = m_{\rm raw} + \mathrm{MAGZERO} - 28

is identical to the corrected magnitude ``magc`` that ``PhotEval`` already
computes from the ``TabPhot/`` catalogs.

.. code-block:: text

    MEF files
        │
        ▼
    MefPhot.py -cat smash          (or -cat gaia for emission-line filters)
        │
        ▼  TabPhot/*.smash.fits
        │
        ├── CalcZeroPoint.py -filter r
        │       ──► zeropoints.fits          (per-file ZP and delta_zp)
        │           Summary/{field}_mef.tab  (ZP columns appended)
        │
        └── PhotEval.py -filter r [-zp_table zeropoints.fits]
                ──► phot_eval_r.fits         (per-star scatter statistics)


Step 1: Forced photometry on MEF files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MefPhot`` measures stellar fluxes at catalog positions in each CCD
extension of every MEF file::

    # r-band: use SMASH (DECam photometric system, deep, good star-galaxy separation)
    MefPhot.py -cat smash -r 6 DECam_MEF/LMC_c42/*.fits

    # Emission-line filters or fields outside LMC/SMC: use Gaia
    MefPhot.py -r 6 DECam_MEF/LMC_c42/*.fits

Output tables are written to ``TabPhot/``:

* ``TabPhot/<root>.smash.fits`` — SMASH-matched photometry
* ``TabPhot/<root>.gaia.fits``  — Gaia-matched photometry

Each table carries ``phot_mag`` (instrumental magnitude, ZP = 28) and
``MAGZERO`` (from the MEF header).

By default ``MefPhot`` restricts photometry to catalog stars in the range
14 < R < 22.  Stars brighter than 14 are likely saturated; stars fainter
than 22 have negligible weight in the calibration fits and would only slow
processing.  The limits are stored in the output header (``MAGBRITE``,
``MAGFAINT``) and can be overridden with ``-mag_bright`` and ``-mag_faint``.

This step needs to be run only once.  The same ``TabPhot/`` files feed both
the consistency checks below and the absolute calibration tools in
:doc:`absolute_photometry`.


Step 2: Per-frame zero-point summary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CalcZeroPoint`` reads the ``TabPhot/`` files and computes a sigma-clipped
median zero point for each file independently of the header ``MAGZERO``::

    CalcZeroPoint.py -filter r      # r-band SMASH files
    CalcZeroPoint.py -filter N662   # Ha, Gaia files
    CalcZeroPoint.py -filter N673   # [SII], Gaia files

Key output columns in ``zeropoints.fits``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Meaning
   * - ``zp_calc``
     - Sigma-clipped median ZP measured from reference stars
   * - ``MAGZERO``
     - Pipeline value from the MEF header
   * - ``delta_zp``
     - ``zp_calc − MAGZERO``: deviation of the header value from measured
   * - ``zp_std``
     - Scatter of per-star ZP estimates within the file (reflects PSF
       and crowding effects, not frame-to-frame calibration)
   * - ``n_stars``
     - Number of reference stars used after sigma-clipping

``CalcZeroPoint`` also updates ``Summary/{field}_mef.tab`` automatically,
adding columns ``ZP_{cat}_{band}``, ``ZP_std_{cat}_{band}``, and
``ZP_n_{cat}_{band}`` (e.g. ``ZP_smash_r``, ``ZP_gaia_g``).

**Interpreting** ``delta_zp``:

A *mean* offset that is the same for all files in a filter indicates a
systematic difference between the reference catalog and the pipeline
photometric system (a color term).  This does not affect how consistently
frames are scaled relative to each other.

A *scatter* in ``delta_zp`` across files — quantified by ``std(delta_zp)``
— indicates that some frames have genuinely different ``MAGZERO`` values.
This is what hurts consistency and is the case where the empirical correction
below is useful.


Step 3: Multi-frame scatter check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``PhotEval`` groups all detections of the same star across overlapping frames
and computes per-star scatter statistics::

    PhotEval.py -filter r -min_n 3 -o phot_eval_r.fits

The printed summary gives an immediate answer:

.. code-block:: text

    --- Summary (N sources) ---
      magc_std    : median=0.051   (header MAGZERO correction)
      mag_err     : median=0.056   (photon noise per measurement)
      chi2_nu_c   : median=0.996   (header MAGZERO)

The key metric is ``chi2_nu_c``: the reduced chi-squared of
``phot_mag + MAGZERO − 28`` relative to the photon-noise expectation.

* **≈ 1** — the ``MAGZERO`` values are mutually consistent at the
  photon-noise level.  No further calibration is needed.
* **Significantly > 1** — there is excess frame-to-frame scatter beyond
  photon noise.  The empirical ZP correction (:ref:`empirical-zp`) may help.

The 90th-percentile ``chi2_nu_c`` characterises the tail of inconsistent
sources (variable stars, blends, or chip-edge artefacts) even when the
bulk of the data is well-behaved.

To compare the header ``MAGZERO`` directly against the empirical ZP from
``CalcZeroPoint``, pass the ``zeropoints.fits`` table::

    PhotEval.py -filter r -zp_table zeropoints.fits -o phot_eval_r_compare.fits

The summary will show both corrections side by side:

.. code-block:: text

    chi2_nu_c   : median=0.996  (header MAGZERO)
    chi2_nu_emp : median=1.024  (empirical ZP)

If ``chi2_nu_emp`` is not lower than ``chi2_nu_c``, the header values are
already as consistent as possible and no further action is needed.


.. _inter-filter:

Inter-Filter Consistency
--------------------------

Intra-filter consistency (checked above with ``PhotEval``) ensures that all
r-band exposures are mutually consistent.  A separate question is whether
the r-band images and the emission-line images (N662, N673) are on the same
flux scale, which is what ``CleanStars`` requires.

``CleanStars`` performs a straight subtraction:

.. math::

   \text{ha\_clean} = \text{N662} - r_{\rm pure}

A star cancels only if it has the same DN in both images.  After ``MefPrep``
has placed all images on the mag-28 scale, this is equivalent to asking
whether the star has the same MAGZERO-corrected magnitude
(:math:`m_c = m_{\rm phot} + \mathrm{MAGZERO} - 28`)
in both filters.

**Diagnostic workflow**::

    # Run MefPhot with Gaia for ALL filters (no -cat smash)
    MefPhot.py -r 6 DECam_MEF/LMC_c42/*.fits

    # Compare r vs Hα
    PhotAnal.py -filter1 r -filter2 N662 TabPhot/*.gaia.fits

    # Compare r vs [SII]
    PhotAnal.py -filter1 r -filter2 N673 TabPhot/*.gaia.fits

    # Compare off-line vs Hα (N708 subtraction path)
    PhotAnal.py -filter1 N708 -filter2 N662 TabPhot/*.gaia.fits

**Interpreting the output**:

The printed summary reports two scatter values:

.. code-block:: text

    --- Inter-filter consistency: r vs N662 ---
      N stars matched          : 3842
      Mean Δmag (r−N662)       : +0.021 mag
      Scatter (total)          : σ = 0.048 mag  → 4.5% typical stellar residual
      Color term b (G−R)       : +0.028 mag/mag
      Scatter (after color)    : σ = 0.031 mag  → 2.9% typical stellar residual  [floor]

* **Mean Δmag**: a systematic offset between the two zero points.  This
  shifts all stars by the same factor and is the dominant error when it is
  non-zero.  If ``|mean Δmag| > 0.05``, the MAGZERO for one filter may
  need adjustment.

* **Scatter (total)**: includes both the color-term spread and any
  calibration noise.

* **Scatter (after color)**: the irreducible floor set by stellar color
  diversity.  Stars of different temperatures have different flux ratios
  between the two filters; this scatter cannot be eliminated by an
  overall scale adjustment.  It sets the minimum stellar residual that
  will remain after ``CleanStars`` subtraction.

The **color term** :math:`b` describes how the filter-to-filter offset
depends on stellar color (Gaia G−R).  A large :math:`b` means that hot
blue stars and cool red stars subtract differently even with the best
mean scale factor.

**Diagnostic plots**: ``PhotAnal`` saves two panels:

1. Δmag vs Gaia G−R color with the best-fit color term (left)
2. Histogram of residuals after removing the color fit (right)

**Output files**:

* ``filter_compare_<f1>_<f2>.fits`` — per-star table with columns
  ``Source_name``, ``magc_1``, ``magc_2``, ``delta_mag``,
  ``G_R`` (Gaia G−R color), ``delta_mag_corrected`` (residual after
  color-term removal).
* ``filter_compare_<f1>_<f2>.png`` — diagnostic plot.

.. note::

   ``PhotAnal`` statistics are unweighted: every star in the magnitude range
   contributes equally to the scatter and color-term fit regardless of
   photometric S/N.  This means the reported ``sigma_residual`` is a slight
   overestimate of the true calibration floor, particularly near the faint
   magnitude limit.  For a cleaner result restrict the Gaia files to a
   brighter limit (e.g. ``-mag_faint 20``) when running ``MefPhot``.


Running the complete check: CheckPhot
--------------------------------------

``CheckPhot.py`` runs the full intra-filter and inter-filter workflow in a
single command and writes two concise summary tables::

    CheckPhot.py LMC_c42
    CheckPhot.py -np 16 LMC_c42           # more parallel workers for MefPhot
    CheckPhot.py -no_mefphot LMC_c42      # skip MefPhot if TabPhot/ already exists

It runs ``MefPhot`` (SMASH + Gaia), ``CalcZeroPoint``, ``PhotEval``, and
``PhotAnal`` in sequence, then prints a two-table summary and writes:

* ``Summary/{field}_phot_check_intra.fits``
* ``Summary/{field}_phot_check_inter.fits``


Interpreting the intra-filter table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    Filter N_files N_stars_zp mean_delta_zp std_delta_zp N_sources chi2_nu_c_med chi2_nu_c_90
    ------ ------- ---------- ------------- ------------ --------- ------------- ------------
         r      60   14199208       -0.0233       0.0220    522529        0.9961       4.3871
      N662      34    7753892        0.0044       0.0580    533454        1.6982      19.5661
      N673      60   13613626        0.0063       0.0503    536902        1.1741      21.2342
      N708      72   14500032        0.0640       0.0616    488251        1.3733      11.7555

**std_delta_zp** (MAD scatter of per-file MAGZERO deviations, robust to
outlier exposures):

* < 0.03 mag — excellent; header ``MAGZERO`` values are mutually consistent.
* 0.03–0.06 mag — moderate; narrowband filters often fall here due to
  atmospheric variability in the filter bandpass.
* > 0.06 mag — investigate; consider the empirical ZP correction below.

**mean_delta_zp** (systematic offset of header MAGZERO from measured ZP):

* A small mean (< 0.02 mag) with low scatter is fine; it reflects a color
  term between the catalog and the DECam system.
* A large mean (> 0.05 mag) for one filter, as seen for N708 above (+0.064),
  means the pipeline MAGZERO is systematically wrong for that filter.  This
  will leave a scale error after ``MefPrep`` that affects CleanStars.
  Check whether the same offset is present in both emission-line filters
  used for subtraction — if the offset is equal it cancels in the
  inter-filter comparison.

**chi2_nu_c_med** (frame-to-frame consistency relative to photon noise):

* ≈ 1 — perfect (r-band in the example); MAGZERO values are consistent at
  the photon-noise level.
* 1.1–1.7 — moderate excess scatter; typical for narrowband filters where
  atmospheric airglow or extended nebulosity affects the background.
* chi2_nu_c_90 >> 1 — the high-chi2 tail in N662 and N673 consists largely
  of genuine emission-line objects (Be stars, symbiotic stars, etc.) varying
  intrinsically in H\ |alpha| / [SII].  This is astrophysics, not a
  calibration problem.

.. |alpha| unicode:: U+03B1


Interpreting the inter-filter table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    Filter1 Filter2 N_stars mean_delta_mag sigma_total sigma_residual residual_pct color_term_b
    ------- ------- ------- -------------- ----------- -------------- ------------ ------------
          r    N662  788480        -0.0279      0.0781         0.0612       5.8025       0.1063
          r    N673  808748        -0.0284      0.0666         0.0510       4.8122       0.0870
       N708    N662  786383        -0.0043      0.0761         0.0610       5.7791      -0.0908
       N708    N673  803152         0.0054      0.0712         0.0537       5.0664      -0.1100

**mean_delta_mag** (systematic inter-filter offset):

* |mean| < 0.03 mag — acceptable; a small multiplicative correction is
  applied implicitly by the mean subtraction in CleanStars.
* |mean| > 0.05 mag — significant; consider adjusting MAGZERO for one of
  the two filters.

**residual_pct** (minimum stellar residual after best mean + color
correction):

* This is the irreducible floor for CleanStars stellar subtraction.  Values
  of 5–6% are typical for DECam broadband vs narrowband comparisons because
  of the stellar color diversity in the LMC/SMC.
* Expressed as a fraction of stellar DN: a 6% residual means the brightest
  stars will leave ±6% flux rings after continuum subtraction.

**color_term_b** (slope of Δmag vs Gaia G−R):

* A positive :math:`b` for r vs N662/N673 means blue stars have a smaller
  r−narrowband offset than red stars; hot OB stars will be
  over-subtracted while cool giants are under-subtracted.
* The sign flips to negative for N708 vs narrowband because N708 lies
  red-ward of the emission lines.

**What the numbers above imply for this field**: the r-band calibration
is excellent.  The N708 systematic offset of +0.064 mag warrants attention,
but since N708 vs N662 shows a near-zero inter-filter mean (−0.004 mag),
the offset is shared by both filters and largely cancels in the
CleanStars subtraction.  The ~5–6% residual floor is intrinsic to the
stellar population and cannot be improved by better photometric calibration.


.. _empirical-zp:

Empirical Zero-Point Correction (optional)
-------------------------------------------

If ``std(delta_zp)`` is larger than ~0.03 mag, or if ``chi2_nu_c`` is
significantly above 1, you can replace the header ``MAGZERO`` with a value
measured directly from reference stars.  This requires a second pass of
``MefPrep``.

``CalcZeroPoint`` writes the empirical ZPs to ``Summary/{field}_mef.tab``
automatically (Steps 1–2 above must already have been run).  Re-run
``MefPrep`` with ``-zp`` pointing at the appropriate column::

    MefPrep.py -np 8 -zp ZP_smash_r LMC_c42      # r-band, SMASH
    MefPrep.py -np 8 -zp ZP_gaia_g  LMC_c42      # emission-line, Gaia

``MefPrep`` falls back to header ``MAGZERO`` for any row where the column
value is absent or outside the range 20–35, with a warning.  Every output
file records the ZP actually used:

* ``ZP_USE`` — zero point applied for flux scaling
* ``ZP_SRC`` — source: the column name (e.g. ``ZP_smash_r``) or ``MAGZERO``

After the second pass, continue the pipeline from ``SetupTile`` as normal.

.. note::

   Use ``ZP_smash_r`` only for r-band images of LMC/SMC fields; SMASH does
   not cover other regions or emission-line filters.  Use ``ZP_gaia_g`` for
   N662, N673, or fields outside the Magellanic Cloud footprint.  Both
   columns can coexist in the same ``_mef.tab`` file.


Output Files
------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - File
     - Created by
     - Contents
   * - ``TabPhot/<root>.smash.fits``
     - MefPhot -cat smash
     - Per-star photometry at SMASH positions; ``phot_mag``, ``MAGZERO``
   * - ``TabPhot/<root>.gaia.fits``
     - MefPhot
     - Per-star photometry at Gaia positions; same structure
   * - ``zeropoints.fits``
     - CalcZeroPoint
     - Per-file: ``zp_calc``, ``MAGZERO``, ``delta_zp``, ``zp_std``,
       ``n_stars``, ``Root``, ``Field``, ``Filter``, ``Exptime``
   * - ``Summary/{field}_mef.tab`` (updated)
     - CalcZeroPoint
     - Appended columns: ``ZP_{cat}_{band}``, ``ZP_std_{cat}_{band}``,
       ``ZP_n_{cat}_{band}``
   * - ``phot_eval_<filter>.fits``
     - PhotEval
     - Per-star: ``n_detect``, ``magc_std``, ``chi2_nu_c``,
       ``mag_err_mean``, and (with ``-zp_table``) ``magc_emp_std``,
       ``chi2_nu_emp``
