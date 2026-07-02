====================
Absolute Photometry
====================

Absolute photometry converts instrumental counts (DN) into physical flux
units (erg s⁻¹ cm⁻² Å⁻¹).  This is needed to compare emission-line
measurements with other observations or theoretical models, and to compute
line luminosities.

Two levels of calibration are available:

* **Magnitude calibration** — places instrumental magnitudes on the Gaia or
  SMASH photometric system.  Useful for checking the pipeline ``MAGZERO``
  against an independent reference and for deriving color terms.

* **Physical flux calibration** — converts DN directly to erg s⁻¹ cm⁻² Å⁻¹
  using Gaia XP spectra as spectroscopic standard references.  This is the
  appropriate calibration for the Hα and [SII] emission-line mosaics.


Step 1: Forced photometry on MEF files
---------------------------------------

All calibration paths start from ``MefPhot``, which measures stellar fluxes
at catalog positions in every CCD extension of each MEF file::

    MefPhot.py -r 6 -b 8 12 DECam_MEF/LMC_c42/*.fits          # Gaia (default)
    MefPhot.py -cat smash -r 6 -b 8 12 DECam_MEF/LMC_c42/*.fits  # SMASH

Key parameters:

* ``-r RADIUS`` — aperture radius in pixels (default 6)
* ``-b INNER OUTER`` — background annulus inner and outer radii (default 8, 12)
* ``-np N`` — parallel processes (default 8)
* ``-cat gaia|smash`` — reference catalog (default ``gaia``)

Output tables are written to ``TabPhot/``::

    TabPhot/<root>.gaia.fits    (Gaia-matched)
    TabPhot/<root>.smash.fits   (SMASH-matched, with -cat smash)

Each table carries ``phot_mag`` (instrumental magnitude, ZP = 28), the
reference catalog magnitudes (``G``, ``R`` for Gaia; ``R``, ``G`` for SMASH),
and file-level metadata in the extension 1 header (``FILTER``, ``EXPTIME``,
``MAGZERO``, ``FIELD``, ``ROOT``, ``DETNAME``).

.. note::

   If you have already run ``MefPhot`` for the relative photometry checks
   (:doc:`relative_photometry`), the same ``TabPhot/`` files are used here.
   You do not need to re-run ``MefPhot``.


Step 2: Magnitude calibration
-------------------------------

``ZeroCalc`` reads one or more ``TabPhot/`` files and fits a zero-point model
relating instrumental magnitudes to reference catalog magnitudes.

*Simple fit* (default):

.. math::

   m_{\rm ref} = m_{\rm inst} + c_0

*Color-corrected fit* (``-color``):

.. math::

   m_{\rm ref} = m_{\rm inst} + c_0 + c_1 \times \mathrm{color}

where the color predictor is chosen to be independent of the target band:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20

   * - Target
     - Catalog
     - Color
   * - R
     - Gaia
     - G − R
   * - G
     - Gaia
     - B − R
   * - R
     - SMASH
     - G − R
   * - G
     - SMASH
     - U − R

The derived zero point is :math:`28 + c_0`.  All fits are weighted by the
per-source photometric uncertainty.

::

    ZeroCalc.py -R TabPhot/*.gaia.fits                  # simple, Gaia R
    ZeroCalc.py -R -color TabPhot/*.gaia.fits           # color-corrected, Gaia R
    ZeroCalc.py -G -color TabPhot/*.gaia.fits           # Gaia G
    ZeroCalc.py -R -smash TabPhot/*.smash.fits          # SMASH R
    ZeroCalc.py -R -color -smash TabPhot/*.smash.fits   # color-corrected, SMASH R

Output filenames encode the band, catalog, and fit mode:

* ``MagZero.<band>.gaia.txt`` / ``MagZero.<band>.gaia.color.txt``
* ``MagZero.<band>.smash.txt`` / ``MagZero.<band>.smash.color.txt``
* ``FigZero/<band>_<root>.<cat>[.color].png`` — per-file diagnostic plots

Each summary table has one row per input file with columns ``Filter``,
``Exptime``, ``Root``, ``MagZero`` (= 28 + c_0), ``c_0``, ``c_1``, ``rms``,
and ``HdrZero`` (the pipeline ``MAGZERO``).

Comparing the derived ``MagZero`` to ``HdrZero`` across many files gives a
direct statistical measure of how well the NOIRLAB pipeline calibration
tracks the true zero point for your dataset:

.. code-block:: python

    from astropy.table import Table
    import numpy as np

    mz = Table.read('MagZero.R.gaia.color.txt', format='ascii.fixed_width_two_line')
    delta = mz['MagZero'] - mz['HdrZero']
    print(f"Mean offset (derived - pipeline): {np.mean(delta):+.3f}")
    print(f"Std  offset                      : {np.std(delta):.3f}")


Step 3: Physical flux calibration
-----------------------------------

For the Hα (N662) and [SII] (N673) emission-line mosaics, the conversion
from DN to physical flux (erg s⁻¹ cm⁻² Å⁻¹) uses Gaia XP spectra as
spectroscopic standards.

**Why XP spectra?**  Most reference catalogs provide broad-band magnitudes
(G, R, …) that must be transformed to the narrow bandpass of each
emission-line filter using synthetic photometry, which introduces color-term
uncertainties.  Gaia XP spectra cover 330–1050 nm at low resolution and
allow a direct synthetic-photometry prediction for any filter.

Sub-step 3a: Cross-match photometry on tile images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run ``PhotCompare`` on the swarped tile images to produce cross-matched tables
linking each detected star to its Gaia source::

    PhotCompare.py -forced DECam_SWARP/LMC_c42/T01/*.fits
    PhotCompare.py -forced DECam_SWARP/LMC_c42/T02/*.fits

Output is written to ``TabPhot/xmatch_<name>.fits`` (one file per tile image)
and diagnostic figures to ``Figs_phot/``.

Sub-step 3b: Derive the DN → flux conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run ``ZeroPoint`` on the cross-matched tables::

    ZeroPoint.py TabPhot/xmatch_*.fits

``ZeroPoint`` retrieves the Gaia XP spectrum for each matched star, integrates
it through the filter bandpass to predict the expected count rate, and
compares it to the measured DN to determine the conversion factor.  Results
are accumulated in ``PhotMaster.txt``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Contents
   * - ``Filename``
     - Tile image filename
   * - ``Ha_flux``
     - DN-to-flux factor at Hα (erg s⁻¹ cm⁻² Å⁻¹ per DN)
     -
   * - ``SII_flux``
     - DN-to-flux factor at [SII] (erg s⁻¹ cm⁻² Å⁻¹ per DN)


Using SMASH as the Reference Catalog
--------------------------------------

For Magellanic Cloud fields the SMASH DR2 catalog (`Nidever et al. 2017
<https://doi.org/10.3847/1538-3881/aa6af6>`_) offers advantages over Gaia
for broadband calibration:

* Its photometry is in the **DECam ugriz system**, so the color term
  :math:`c_1` is near zero for r-band images.
* It is **deeper** (r ~ 22) than Gaia (G ~ 20), giving more reference stars
  per CCD in the crowded LMC/SMC fields.
* It includes a stellar probability parameter (``prob``) for star-galaxy
  separation.

SMASH is **not all-sky** and has no coverage for the emission-line filters
(N662, N673).  ``ZeroPoint`` (spectroscopic flux calibration) has no SMASH
equivalent and always uses Gaia XP spectra.

PSF Photometry (crowded fields)
---------------------------------

In regions where stellar crowding makes aperture photometry unreliable, a
PSF-fitting path is available.

1. **Detect sources** with ``StarFind``::

       StarFind.py image.fits

2. **Build a PSF model** with ``PsfBuild``::

       PsfBuild.py image.fits TabPhot/<root>_all_stars.fits

3. **Fit the PSF** with ``PsfPhot``::

       PsfPhot.py image.fits <root>_psf_model.fits TabPhot/<root>_all_stars.fits

   Output is a photometry table with the same structure as ``MefPhot`` output
   and can be passed directly to ``ZeroCalc``.

``MefPhot`` (forced at fixed catalog positions) is preferred for zero-point
determination because the source list is unambiguous.  ``PsfPhot`` is
preferred when completeness or accurate fluxes in crowded regions matter.


Output Files
------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - File
     - Created by
     - Contents
   * - ``TabPhot/<name>.gaia.fits``
     - MefPhot
     - Per-star photometry at Gaia positions
   * - ``TabPhot/<name>.smash.fits``
     - MefPhot -cat smash
     - Per-star photometry at SMASH positions
   * - ``MagZero.<band>.gaia.txt``
     - ZeroCalc
     - Per-file ZP fit vs Gaia: ``c_0``, ``c_1``, ``rms``, ``HdrZero``
   * - ``MagZero.<band>.smash.txt``
     - ZeroCalc -smash
     - Per-file ZP fit vs SMASH
   * - ``FigZero/``
     - ZeroCalc
     - Per-file diagnostic plots
   * - ``TabPhot/xmatch_<name>.fits``
     - PhotCompare
     - Stars cross-matched between tile images and Gaia
   * - ``Figs_phot/``
     - PhotCompare
     - Comparison plots (DECam vs catalog magnitudes)
   * - ``PhotMaster.txt``
     - ZeroPoint
     - Per-tile DN → erg s⁻¹ cm⁻² Å⁻¹ at Hα and [SII]
   * - ``GAIA/``
     - GaiaCat
     - Cached Gaia source catalogs
