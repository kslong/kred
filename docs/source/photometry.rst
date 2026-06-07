==============================
Photometry and Calibration
==============================

This page describes the tools available for stellar photometry and photometric
calibration.  There are two main goals: (1) measuring instrumental fluxes of
stars using a reference catalog — either Gaia DR3 or, for Magellanic Cloud
fields, the SMASH DR2 catalog — and (2) placing those measurements on an
absolute calibration scale, either in magnitudes (by comparison to broadband
photometry from the reference catalog) or in physical flux units (by
comparison to Gaia XP spectra).  A secondary goal is to check whether the
``MAGZERO`` value carried in each MEF file header is consistent with an
independently derived zero point.


Overview
--------

The photometry system comprises the following modules:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Module
     - Purpose
   * - :doc:`GaiaCat <api/GaiaCat/index>`
     - Retrieve Gaia DR3 source positions and photometry for a field
   * - :doc:`Smash <api/Smash/index>`
     - Retrieve SMASH DR2 catalog for Magellanic Cloud fields (alternative to Gaia)
   * - :doc:`MefPhot <api/MefPhot/index>`
     - Forced aperture photometry on MEF images at Gaia or SMASH catalog
       positions (``-cat gaia|smash``); primary source of ``TabPhot/`` files
   * - :doc:`CalcZeroPoint <api/CalcZeroPoint/index>`
     - Quick per-file zero-point estimate (sigma-clipped median) from
       ``TabPhot/`` files; produces a multi-file summary table with
       ``delta_zp`` relative to the MEF header ``MAGZERO``
   * - :doc:`ZeroCalc <api/ZeroCalc/index>`
     - Weighted regression fit (with optional color term) of the zero point
       against Gaia or SMASH broadband magnitudes; preferred when a color
       correction is scientifically important
   * - :doc:`PhotEval <api/PhotEval/index>`
     - Multi-frame photometric consistency check: stacks all ``TabPhot/``
       catalogs for a given filter, groups by SMASH ``Source_name``, and
       computes per-source scatter statistics (both raw and ZP-corrected)
   * - :doc:`PhotCompare <api/PhotCompare/index>`
     - Aperture photometry on tile/swarped images with catalog cross-match and
       comparison plots (``-cat gaia|smash``)
   * - :doc:`ZeroPoint <api/ZeroPoint/index>`
     - Derive a physical flux zero point using Gaia XP spectra
   * - :doc:`PhotAnal <api/PhotAnal/index>`
     - Multi-filter consistency check across a smaller set of files
   * - :doc:`StarFind <api/StarFind/index>`
     - Automatically detect stars for PSF construction
   * - :doc:`PsfBuild <api/PsfBuild/index>`
     - Select PSF stars and build Gaussian / Moffat / empirical PSF models
   * - :doc:`PsfPhot <api/PsfPhot/index>`
     - PSF-fitting photometry for crowded fields


Calibration Concepts
--------------------

The DECam community pipeline assigns a photometric zero point to each
exposure, stored as the ``MAGZERO`` keyword in the primary FITS header of
each MEF file.  This value is used downstream (e.g., in ``MefPrep``) to
place images on a common flux scale.  The routines described here allow you
to independently derive a zero point from reference catalog stars (Gaia, or
SMASH for Magellanic Cloud fields) and compare it to the pipeline-supplied
value.

Three calibration / evaluation paths are available:

**Quick zero-point summary** (``CalcZeroPoint``)
    For each ``TabPhot/`` file, computes ``zp_i = 28 + (ref_mag - phot_mag)``
    for every star, then takes the sigma-clipped median.  Produces a single
    summary table (``zeropoints.fits``) with one row per file, including
    ``delta_zp = zp_calc − MAGZERO``.  Use this for a rapid sanity check of
    the pipeline-supplied ``MAGZERO`` values across an entire observing run.

**Magnitude calibration by regression** (``ZeroCalc``)
    Fits instrumental magnitudes (measured at a reference zero point of 28)
    against Gaia or SMASH broadband magnitudes.  Two modes are available:

    *Simple fit* (default):

    .. math::

        m_{\rm ref} = m_{\rm inst} + c_0

    *Color-corrected fit* (``-color`` flag):

    .. math::

        m_{\rm ref} = m_{\rm inst} + c_0 + c_1 \,\times\, \mathrm{color}

    where the color predictor is independent of the target band (Gaia R: G−R;
    Gaia G: B−R; SMASH R: G−R; SMASH G: U−R).  All fits are weighted by the
    per-source photometric uncertainty so that bright, well-measured stars
    dominate the solution.  The derived zero point is :math:`28 + c_0` and
    can be compared directly to ``MAGZERO``.

**Physical flux calibration** (``ZeroPoint``)
    Uses Gaia XP spectra to predict the expected flux at the Hα (6563 Å)
    and [SII] (6720 Å) wavelengths for each reference star, then compares
    these spectroscopic fluxes to the measured instrumental counts to derive
    the factor converting 1 DN to physical flux units (erg s⁻¹ cm⁻² Å⁻¹).
    This is the appropriate calibration for the emission-line images produced
    by this pipeline.


Standard Workflow
-----------------

.. _photometry-workflow:

The typical workflow is illustrated below.  The ``GaiaCat`` module is called automatically by ``MefPhot`` and
``PhotCompare`` when using the default Gaia catalog; the ``Smash`` module is
called automatically by either script when ``-cat smash`` is used.  You do
not need to invoke either catalog module separately.

.. code-block:: text

    MEF file(s)
        │
        ▼
    MefPhot.py [-cat gaia|smash]
        │
        ├──► TabPhot/<name>.gaia.fits
        │    TabPhot/<name>.smash.fits
        │
        ├── CalcZeroPoint.py [-filter F]   ──► zeropoints.fits
        │       (quick sigma-clipped median ZP per file; delta_zp vs MAGZERO)
        │
        ├── ZeroCalc.py [-G|-R] [-color] [-smash]
        │       ──► MagZero.<band>.<cat>[.color].<date>.txt
        │           FigZero/<band>_<name>.<cat>[.color].png
        │       (regression ZP fit with optional color term, per file)
        │
        └── PhotEval.py [-filter F] [-exp T]
                ──► phot_eval_<filter>[_<exp>s].fits
                (per-source scatter across all frames; raw + ZP-corrected)

    MEF or tile image(s)
        │
        ▼
    PhotCompare.py [-cat gaia|smash]  ──► TabPhot/xmatch_<name>.fits
                                           Figs_phot/<name>.png
        │
        ▼
    ZeroPoint.py         ──► PhotMaster.txt
                               (DN→flux conversion for Hα and [SII])

    Multiple TabPhot files
        │
        ▼
    PhotAnal.py          ──► summary table + residual histogram


Step-by-Step Instructions
--------------------------

Step 1: Run forced photometry on MEF files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MefPhot`` reads each CCD extension, queries Gaia for sources in the field
of view (caching results in ``GAIA/``), and measures fluxes at each Gaia
position using circular aperture photometry with sigma-clipped annulus
background subtraction.  The output FITS table contains one row per source
per CCD and carries the MEF primary header in extension 0, making it easy
to compare derived quantities with the original header values.

::

    MefPhot.py -r 6 -b 8 12 DECam_MEF/c4d_241122_*.fits          # Gaia (default)
    MefPhot.py -cat smash -r 6 -b 8 12 DECam_MEF/c4d_241122_*.fits  # SMASH

Key parameters:

* ``-r RADIUS`` – aperture radius in pixels (default 6); also sets the
  background annulus to RADIUS+1 and RADIUS+4
* ``-b INNER OUTER`` – override background annulus radii explicitly
* ``-np N`` – number of parallel processes (default 8)
* ``-cat gaia|smash`` – reference catalog (default ``gaia``)

Output tables are written to ``TabPhot/`` with a catalog suffix in the
filename to prevent overwriting when both catalogs are used::

    TabPhot/<root>.gaia.fits
    TabPhot/<root>.smash.fits

Each table includes ``phot_mag`` (instrumental magnitude assuming ZP = 28),
``MAGZERO`` (from the MEF header), the reference catalog magnitudes ``G``
and ``R``, and a ``Catalog`` column (``'Gaia'`` or ``'SMASH'``).

Step 2: Derive the magnitude zero point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ZeroCalc`` reads one or more ``TabPhot/`` files and fits the model
described above.  The default is a simple one-parameter fit; add ``-color``
for a two-parameter color-corrected fit.  All fits are weighted by the
per-source photometric uncertainty.

::

    ZeroCalc.py -R TabPhot/*.gaia.fits                  # simple, Gaia R
    ZeroCalc.py -R -color TabPhot/*.gaia.fits           # color-corrected, Gaia R
    ZeroCalc.py -G -color TabPhot/*.gaia.fits           # color-corrected, Gaia G (uses B-R)
    ZeroCalc.py -R -smash TabPhot/*.smash.fits          # simple, SMASH R
    ZeroCalc.py -R -color -smash TabPhot/*.smash.fits   # color-corrected, SMASH R

Output filenames encode the band, catalog, and fit mode:

* ``MagZero.<band>.gaia.txt`` / ``MagZero.<band>.gaia.color.txt`` – simple
  and color-corrected Gaia runs; one row per input file with columns
  ``Filter``, ``Exptime``, ``Root``, ``MagZero`` (= 28 + c_0), ``c_0``,
  ``c_1``, ``rms``, ``HdrZero`` (pipeline MAGZERO), ``Catalog``, ``Filename``
* ``MagZero.<band>.smash.txt`` / ``MagZero.<band>.smash.color.txt`` – same
  for SMASH runs
* ``FigZero/<band>_<root>.gaia[.color].png`` /
  ``FigZero/<band>_<root>.smash[.color].png`` – per-file diagnostic plots

Step 2b: Quick zero-point summary with CalcZeroPoint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CalcZeroPoint`` provides a faster alternative to ``ZeroCalc`` for obtaining
a per-file zero-point estimate.  Instead of a regression fit it takes the
sigma-clipped median of the per-star zero-point estimates and reports the
scatter and number of stars used.  Running with no arguments processes all
``*.smash.fits`` and ``*.gaia.fits`` files in ``TabPhot/``::

    CalcZeroPoint.py                     # all filters
    CalcZeroPoint.py -filter r           # r-band only
    CalcZeroPoint.py -filter N662 -snr 20

The output table ``zeropoints.fits`` has one row per input file.  The
``delta_zp`` column (= ``zp_calc − MAGZERO``) directly answers the question
"how far is the pipeline zero point from what we measure?".  A consistent
offset across many files in the same filter suggests a systematic difference
between the reference catalog's photometric system and the pipeline, rather
than file-by-file variation.

.. note::

   **CalcZeroPoint vs ZeroCalc**: use ``CalcZeroPoint`` for a rapid overview
   and to spot outlier files; use ``ZeroCalc`` when you need a color-corrected
   fit or want per-file diagnostic plots.  Both read the same ``TabPhot/``
   files and can be run independently.

Step 2c: Multi-frame consistency check with PhotEval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``PhotEval`` addresses a different question: not what the zero point of each
frame is, but whether sources measured in multiple overlapping frames give
consistent fluxes.  It stacks all ``TabPhot/`` catalogs for a chosen filter
and exposure time, groups detections by SMASH ``Source_name``, and computes
per-source statistics::

    PhotEval.py -filter r -exp 30 -snr 10 -min_n 5
    PhotEval.py -filter N662 -min_n 3 -o ha_eval.fits

The output table (e.g. ``phot_eval_r_30s.fits``) contains one row per SMASH
source with the following key columns:

* ``mag_std`` / ``magc_std`` – scatter of raw and ZP-corrected magnitudes
  across frames.  ``magc_std`` (using ``MAGZERO`` from the MEF header to
  correct each detection) is the primary repeatability metric.
* ``mag_err_mean`` – expected photon-noise magnitude error (2.5/ln10 × ErrNet/Net).
* ``chi2_nu`` / ``chi2_nu_c`` – reduced chi-squared of raw and corrected
  magnitudes.  Values near 1 indicate photon-noise-limited scatter; values
  significantly larger than 1 indicate excess frame-to-frame variation.

A ``chi2_nu_c`` median near 1 after applying the ``MAGZERO`` correction means
the header zero points are accurately capturing the frame-to-frame throughput
variation.  A persistent excess indicates either an additional systematic (PSF
variations, flat-fielding) or astrophysical variability in the source sample.

Step 3: Compare to the MEF header MAGZERO
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The derived zero point is :math:`28 + c_0` (from the model
``m_ref = m_inst + c_0 + c_1(G-R)`` with ``m_inst = 28 - 2.5\,\log_{10}(\text{flux})``).
The ``HdrZero`` column in the ``ZeroCalc`` summary table already carries the
pipeline ``MAGZERO`` read from the MEF header, so a comparison is
straightforward:

.. code-block:: python

    from astropy.table import Table
    import numpy as np

    mz = Table.read('MagZero.G.gaia.txt', format='ascii.fixed_width_two_line')
    mz['ZP_derived'] = 28.0 + mz['c_0']
    mz['delta'] = mz['ZP_derived'] - mz['HdrZero']

    print(mz['Root', 'Filter', 'HdrZero', 'ZP_derived', 'delta', 'rms'])
    print(f"\nMean offset (derived - pipeline): {np.mean(mz['delta']):+.3f}")
    print(f"Std  offset                      : {np.std(mz['delta']):.3f}")

If many files have been processed, the summary table can be cross-matched
against ``MefSum`` output (which carries ``MAGZERO`` per file) to produce a
statistical comparison across an entire observing campaign.

Step 4 (optional): Multi-filter consistency check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``PhotAnal`` reads several ``TabPhot/`` files from different filters and
reports the median offset and scatter of instrumental magnitudes relative
to the Gaia R band.  This is useful for identifying nights or CCDs with
anomalous zero points.

::

    PhotAnal.py TabPhot/c4d_241122_*.fits

Step 5 (optional): Physical flux calibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Hα and [SII] emission-line images the physical flux corresponding
to 1 DN is needed.  This requires Gaia XP spectra.

First, run ``PhotCompare`` to produce cross-matched tables:

::

    PhotCompare.py -forced DECam_SWARP/LMC_c42/T01/*.fits

Then run ``ZeroPoint`` on the resulting ``xmatch_*.fits`` tables:

::

    ZeroPoint.py TabPhot/xmatch_*.fits

Output is accumulated in ``PhotMaster.txt``, which contains one row per
image giving the flux (in erg s⁻¹ cm⁻² Å⁻¹) corresponding to 1 DN at
the Hα and [SII] wavelengths.


PSF Photometry (crowded fields)
--------------------------------

For fields where stellar crowding makes aperture photometry unreliable, a
PSF-fitting path is available.  The PSF-fitted fluxes can also be fed into
``ZeroCalc`` for zero-point determination.

1. **Detect sources** with ``StarFind``::

       StarFind.py image.fits

   Outputs ``TabPhot/<root>_all_stars.fits`` with detected positions and
   aperture photometry quality metrics.

2. **Build a PSF model** with ``PsfBuild``::

       PsfBuild.py image.fits TabPhot/<root>_all_stars.fits

   ``PsfBuild`` selects stars using a weighted quality score (SNR 35%,
   FWHM 25%, eccentricity 20%, background contamination 10%, concentration
   10%) and constructs Gaussian, Moffat, and empirical PSF models.

3. **Fit the PSF** with ``PsfPhot``::

       PsfPhot.py image.fits <root>_psf_model.fits TabPhot/<root>_all_stars.fits

   Output is a photometry table compatible with ``ZeroCalc``.


.. note::

   ``MefPhot`` (forced at Gaia positions) and ``PsfPhot`` (PSF fitting at
   detected positions) serve different scientific goals.  ``MefPhot`` is
   preferred for zero-point determination because the source list is fixed
   and unambiguous.  ``PsfPhot`` is preferred when completeness or accurate
   fluxes in crowded regions are needed.


Using SMASH as the Reference Catalog
-------------------------------------

For Magellanic Cloud fields the SMASH DR2 catalog (Survey of the MAgellanic
Stellar History, `Nidever et al. 2017 <https://doi.org/10.3847/1538-3881/aa6af6>`_)
is an attractive alternative to Gaia because:

* Its photometry is in the **DECam ugriz system** — the same system as the
  images being calibrated.  For r-band images the color term :math:`c_1`
  should be very close to zero.
* It is **deeper** (r ~ 22) than Gaia (G ~ 20), giving more reference
  stars per CCD in dense LMC/SMC fields.
* It uses an explicit stellar probability parameter (``prob``) for
  star-galaxy separation, which helps in the crowded Magellanic Cloud
  background.

SMASH is **not all-sky**.  For fields outside the LMC/SMC footprint, or for
any filter other than the standard ugriz set, Gaia remains the only option.
``ZeroPoint.py`` (spectroscopic flux calibration using Gaia XP spectra) has
no SMASH equivalent and is unaffected by this choice.

Retrieving a SMASH catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Smash.py`` now provides a ``get_smash(ra, dec, size)`` function with the
same interface as ``GaiaCat.get_gaia``.  It wraps ``Smash.do_one`` and caches
the result in the ``Smash/`` subdirectory, so repeated calls for the same
field are fast.  You can also retrieve a catalog manually for inspection::

    Smash.py -rad 0.5 -plot 81.9 -69.0    # LMC centre, 0.5 deg radius

Running MefPhot with SMASH
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass ``-cat smash`` to ``MefPhot``::

    MefPhot.py -cat smash -r 6 DECam_MEF/c4d_*.fits

The output tables are identical in structure to Gaia-based runs.  A
``Catalog`` column (value ``'SMASH'``) is added so downstream tools and the
output FITS header record which reference was used.

Running PhotCompare with SMASH
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass ``-cat smash`` to ``PhotCompare`` to perform photometry and generate
comparison figures using SMASH catalog positions::

    PhotCompare.py -cat smash -forced DECam_SWARP/LMC_c42/T01/*.fits

The cross-matched output tables written to ``TabPhot/`` are otherwise
identical in structure to Gaia-based runs, and can be passed directly to
``ZeroPoint.py`` for flux calibration.

Deriving zero points from SMASH
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass ``-smash`` to ``ZeroCalc``::

    ZeroCalc.py -R -smash TabPhot/*.fits   # fit to SMASH R
    ZeroCalc.py -G -smash TabPhot/*.fits   # fit to SMASH G (DECam g)

The output summary file is named ``MagZero.<band>.smash.txt`` (rather than
``MagZero.<band>.gaia.txt``) so Gaia and SMASH runs do not overwrite each
other.  The ``Catalog`` column in the summary table records the source.

Comparing Gaia and SMASH zero points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running both pipelines on the same images provides an internal consistency
check.  Because SMASH R ≈ DECam r, a well-calibrated image should give a
smaller :math:`c_1` with SMASH than with Gaia and similar :math:`c_0`:

.. code-block:: python

    from astropy.table import Table
    import numpy as np

    gaia  = Table.read('MagZero.R.gaia.txt',  format='ascii.fixed_width_two_line')
    smash = Table.read('MagZero.R.smash.txt', format='ascii.fixed_width_two_line')

    # Match on Root filename
    from astropy.table import join
    both = join(gaia, smash, keys='Root', table_names=['gaia', 'smash'])

    delta_zp = both['MagZero_gaia'] - both['MagZero_smash']
    print(f"Mean ZP difference (Gaia - SMASH): {np.mean(delta_zp):+.3f}")
    print(f"Std  ZP difference               : {np.std(delta_zp):.3f}")

A large systematic offset between the two would indicate either a problem
with the SMASH coverage in that field or a significant systematic in the
community pipeline MAGZERO for that filter.


Output File Summary
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File
     - Created by
     - Contents
   * - ``TabPhot/<name>.gaia.fits``
     - MefPhot (default)
     - Per-source photometry using Gaia positions; MEF header in
       extension 0; extraction parameters in extension 1 header
   * - ``TabPhot/<name>.smash.fits``
     - MefPhot -cat smash
     - Per-source photometry using SMASH positions; same structure
   * - ``TabPhot/xmatch_<name>.fits``
     - PhotCompare, StarFind, PsfPhot
     - Cross-matched photometry tables (no catalog suffix)
   * - ``zeropoints.fits``
     - CalcZeroPoint
     - Per-file zero-point summary (sigma-clipped median): zp_calc, zp_wmean,
       zp_std, zp_err, zp_mad, n_stars, n_total, delta_zp
   * - ``phot_eval_<filter>[_<exp>s].fits``
     - PhotEval
     - Per-source scatter statistics across all frames: n_detect, mag_mean,
       mag_std, magc_mean, magc_std, mag_err_mean, chi2_nu, chi2_nu_c
   * - ``MagZero.<band>.gaia.txt``
     - ZeroCalc (default)
     - Per-file zero-point fit results against Gaia (c_0, c_1, rms, Catalog)
   * - ``MagZero.<band>.smash.txt``
     - ZeroCalc -smash
     - Per-file zero-point fit results against SMASH (c_0, c_1, rms, Catalog)
   * - ``PhotMaster.txt``
     - ZeroPoint
     - Per-image physical flux zero points at Hα and [SII]
   * - ``Figs_phot/``
     - PhotCompare, ZeroCalc
     - Diagnostic comparison plots
   * - ``GAIA/``
     - GaiaCat
     - Cached Gaia source catalogs (auto-populated)
