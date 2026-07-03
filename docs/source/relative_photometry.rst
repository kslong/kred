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
single command::

    CheckPhot.py LMC_c42
    CheckPhot.py -np 16 LMC_c42           # more parallel workers for MefPhot
    CheckPhot.py -no_mefphot LMC_c42      # skip MefPhot if TabPhot/ already exists
    CheckPhot.py -o mydir LMC_c42         # write all outputs to mydir/ instead
    CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits   # reprint summary only

It runs ``MefPhot`` (SMASH + Gaia), ``CalcZeroPoint``, ``PhotEval``, and
``PhotAnal`` in sequence.  All output files are collected in a single
``CheckPhot/`` directory (override with ``-o DIR``) to avoid cluttering
the working directory.  All filenames include the field name so multiple
fields share the directory without collision.  The final summary is a
single FITS file with two named table extensions::

    CheckPhot/LMC_c42_phot_check.fits   # ext INTRA + INTER


Image groups
^^^^^^^^^^^^^

``CheckPhot`` reads ``config/DeMCELS_images.txt`` (or the file named by
``-config``) to determine which exposures belong together.  Each named
*Image group* identifies a set of exposures that will be co-added into one
mosaic tile:

.. code-block:: text

    # Image   FILTER  EXPTIME
    N662       N662    400
    N662       N662    600
    N662       N662    800
    N662_s     N662     30
    N662_s     N662     60
    r          r        30
    r          r        60
    r_s        r         5

Groups whose label ends in ``_s`` are *short-exposure* (shallow) groups;
all others are *long-exposure*.  The INTRA table has one row per Image
group that has matching files in ``TabPhot/``.  The INTER table compares
every pair of groups within the same tier (long or short) that have
different filters — these are auto-generated, so no filter list needs to
be hard-coded.


Intermediate products in ``CheckPhot/``:

* ``{field}_zeropoints_{label}.fits``         — CalcZeroPoint per-file ZP table
* ``{field}_phot_eval_{label}.fits``          — PhotEval per-source chi² table
* ``{field}_filter_compare_{l1}_{l2}.fits``   — PhotAnal matched-star table

To read the summary tables in Python::

    from astropy.table import Table
    intra = Table.read('CheckPhot/LMC_c42_phot_check.fits', hdu='INTRA')
    inter = Table.read('CheckPhot/LMC_c42_phot_check.fits', hdu='INTER')

The ``-summary`` flag re-reads an existing FITS file and reprints the
formatted comparison without rerunning any photometry::

    CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits

The same output can be generated from Python::

    import CheckPhot
    CheckPhot.summarize('CheckPhot/LMC_c42_phot_check.fits')

or, if you already have the tables in memory::

    from astropy.table import Table
    intra = Table.read('CheckPhot/LMC_c42_phot_check.fits', hdu='INTRA')
    inter = Table.read('CheckPhot/LMC_c42_phot_check.fits', hdu='INTER')
    CheckPhot._print_summary(intra, inter, 'LMC_c42')


Interpreting the intra-filter table (INTRA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The INTRA extension has one row per active Image group.  Zero-point
columns are named ``{stat}_{catalog}_{ref}`` where ``{catalog}`` is
``smash`` or ``gaia`` and ``{ref}`` is the reference band (``r``, ``g``).
The per-source consistency columns are at the right.

.. code-block:: text

  Image      Filt     N_src  chi2_emp   MagZero    GaiaG   SmashR  best
  ---------- ----- --------  --------  -------- -------- --------  -------
  N662       N662    525648     1.913    0.0123   0.0441   0.0131  MagZero/SmashR
  N662_s     N662    225356     0.813    0.0362   0.0450   0.0341  MagZero/SmashR
  N673       N673    527994     1.112    0.0113   0.0161   0.0104  MagZero/SmashR
  N673_s     N673    270392     0.762    0.0146   0.0188   0.0129  MagZero/SmashR
  N708       N708    473965     1.329    0.0106   0.0508   0.0148  MagZero
  N708_s     N708    280882     0.921    0.0634   0.0259   0.0146  SmashR
  r          r       511339     0.926    0.0104   0.0131   0.0098  MagZero/SmashR
  r_s        r       336795     0.853    0.0149   0.0420   0.0129  MagZero/SmashR

The three ``sigma_c_*`` columns show the per-star scatter of
:math:`m_c = m_{\rm phot} + \mathrm{ZP_{ref}} - 28` across all frames
of the same Image group, for three different choices of zero-point
reference.  They answer the question: *if we scale all frames of this
group to ZP = 28 using reference X, how consistently does the same star
come out?*

**sigma_c_magzero** — uses the header ``MAGZERO`` (pipeline default).

**sigma_c_gaia_g** — uses a ZP measured directly from Gaia G magnitudes
for each frame.  For broadband (r) filters this is competitive; for
narrowband filters the broad Gaia G passband introduces a color-term
scatter that makes it worse than ``MAGZERO``.

**sigma_c_smash_r** — uses a ZP measured from SMASH r magnitudes.  This
is competitive with ``MAGZERO`` for most groups and wins clearly when
``MAGZERO`` is inconsistent (e.g. N708_s, where sigma_c_magzero = 0.063
vs sigma_c_smash_r = 0.015).

**std_dzp_{ref}** (MAD scatter of per-file MAGZERO deviations, robust to
outlier exposures):

* < 0.03 mag — excellent; header ``MAGZERO`` values are mutually consistent.
* 0.03–0.06 mag — moderate; narrowband filters often fall here due to
  atmospheric variability in the filter bandpass.
* > 0.06 mag — investigate; consider the empirical ZP correction below.

**mean_dzp_{ref}** (systematic offset of header MAGZERO from measured ZP):

* A small mean (< 0.02 mag) with low scatter is fine; it reflects a color
  term between the catalog and the DECam system.
* A large mean (> 0.05 mag) means the pipeline MAGZERO is systematically
  wrong for that group.  Check whether the same offset is present in both
  filters used for subtraction — if the offset is equal it cancels in the
  inter-filter comparison.

**chi2_nu_c_emp_med** (frame-to-frame consistency relative to photon noise,
empirical noise model):

* ≈ 1 — frames are mutually consistent at the photon-noise level.
* 1.1–1.7 — moderate excess scatter; typical for narrowband filters where
  atmospheric airglow or extended nebulosity affects the background.
* ``chi2_nu_c_90 >> 1`` — the high-chi2 tail in N662 and N673 consists
  largely of genuine emission-line objects (Be stars, symbiotic stars, etc.)
  varying intrinsically in H\ |alpha| / [SII].  This is astrophysics, not a
  calibration problem.

.. |alpha| unicode:: U+03B1


Interpreting the inter-filter table (INTER)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The INTER extension has one row per Image-group pair.  Pairs are formed
automatically within each tier: all combinations of long-exposure groups
with different filters, and separately all combinations of short-exposure
groups with different filters.  The columns ``Image1``, ``Image2``,
``Filter1``, ``Filter2`` identify each pair.

.. code-block:: text

  Pair                    N_stars     ---- sigma_c ----       resid_MZ  resid%    ---- std_dzp ----
                                    MagZero    GaiaG   SmashR                    MagZero    GaiaG   SmashR
  ---------------------- --------  -------- -------- --------  --------  ------  -------- -------- --------
  r × N662                 780502    0.0724   0.0726   0.0723    0.0568     5.4%    0.0257   0.0159   0.0370
  r × N673                 803780    0.0617   0.0617   0.0617    0.0471     4.4%    0.0166   0.0148   0.0233
  r × N708                 800114    0.1126   0.1128   0.1125    0.0793     7.6%    0.0215   0.0124   0.0293
  N662 × N673              780732    0.0323   0.0328   0.0322    0.0325     3.0%    0.0219   0.0130   0.0360
  N662 × N708              782266    0.0726   0.0732   0.0721    0.0571     5.4%    0.0263   0.0122   0.0403
  N673 × N708              800184    0.0628   0.0630   0.0627    0.0465     4.4%    0.0207   0.0138   0.0268

**sigma_c_magzero / sigma_c_gaia_g / sigma_c_smash_r** (scatter of
:math:`m_{c,F1} - m_{c,F2}` per matched star, for three ZP references):

These three columns are the inter-filter analogue of the same-named
columns in the INTRA table.  The meaning is: *after correcting both
filters to ZP = 28 using reference X, how much does the same star's
flux ratio scatter from star to star?*  For the long-exposure LMC_c42
pairs above, all three references give essentially identical values —
the dominant scatter is stellar color diversity, not calibration noise.

**sigma_residual / residual_pct** (continuum subtraction floor):

The scatter that remains after fitting and subtracting the best-fit
color term (Gaia G−R).  This is the irreducible floor for ``CleanStars``
stellar subtraction.  Values of 4–8% are typical for DECam broadband vs
narrowband comparisons because of the stellar color diversity in the
LMC/SMC.  A 6% residual means the brightest stars leave ±6% flux rings
after continuum subtraction.  This floor cannot be improved by better
photometric calibration — it is set by the physics of the stellar
population.

**mean_delta_mag** (systematic inter-filter offset):

* |Δmag| < 0.03 mag — acceptable; a small multiplicative correction is
  applied implicitly by the mean subtraction in CleanStars.
* |Δmag| > 0.05 mag — significant; consider adjusting MAGZERO for one of
  the two filters.

.. |Δmag| replace:: \|mean_delta_mag\|

**color_term_b** (slope of Δmag vs Gaia G−R):

* A positive :math:`b` for r vs N662/N673 means blue stars have a smaller
  r−narrowband offset than red stars; hot OB stars will be
  over-subtracted while cool giants are under-subtracted.
* The sign flips to negative for N708 vs narrowband because N708 lies
  red-ward of the emission lines.

**std_dzp_magzero / std_dzp_gaia_g / std_dzp_smash_r** (image-to-image
variation in the inter-filter mean offset, for three ZP references):

Each (Image1 file, Image2 file) pair contributes one mean Δmag; ``std_dzp``
is the MAD scatter of those per-pair means.  This is the inter-filter
analogue of ``std_dzp`` in the INTRA table: it measures how consistently
*any two specific images* from the two groups are scaled relative to each
other.

* Values ≲ 0.02 mag mean image-to-image calibration is stable.
* Values > 0.05 mag suggest individual image pairs are poorly matched and
  the continuum subtraction will vary from one set of images to another.

Note that for the long-exposure pairs above, Gaia G gives the smallest
``std_dzp`` (~0.012–0.016 mag vs ~0.017–0.026 for MAGZERO), though the
difference is small compared with ``sigma_residual``.


Printed summary
^^^^^^^^^^^^^^^^^

At the end of every ``CheckPhot`` run (and when called with ``-summary``),
a formatted comparison table is printed that places the three ZP references
side by side and marks which is best for each Image group.  The INTRA
section shows ``sigma_c`` for each ZP reference and the empirical chi2.
The INTER section shows ``sigma_c`` and ``resid_MZ`` (the continuum
subtraction floor using MAGZERO), with pairs sorted so that broadband ×
narrowband pairs appear first, r listed first within each pair, and
short-exposure pairs grouped separately below the long-exposure pairs.

**Choosing a zero-point reference**: for most Image groups and filter pairs,
``MAGZERO`` and ``SMASH r`` give essentially identical ``sigma_c``, and
``GAIA G`` is noticeably worse for narrowband filters.  The main exception
is any Image group where the MAGZERO values are internally inconsistent
(``sigma_c_magzero >> sigma_c_smash_r``); in that case ``SMASH r`` is
the better choice for MefPrep rescaling.  The ``best`` column in the
printed INTRA summary flags this automatically.

The continuum subtraction floor (``residual_pct``) is the same regardless
of which ZP reference is used, because it is set by stellar color diversity
rather than calibration accuracy.


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
   * - ``CheckPhot/{field}_phot_check.fits``
     - CheckPhot
     - Two-extension FITS: ``INTRA`` (one row per Image group) and
       ``INTER`` (one row per Image-group pair).  Contains
       ``sigma_c_magzero``, ``sigma_c_gaia_g``, ``sigma_c_smash_r``,
       ``std_dzp_*``, ``sigma_residual``, ``residual_pct``
