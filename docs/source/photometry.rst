==========
Photometry
==========

Although the scientific goal of the KRED pipeline is to map extended
emission-line nebulae, stellar photometry plays a central role at two
distinct stages of the reduction.

**Relative photometry** — ensuring that all images are on the same flux
scale — is a prerequisite both for combining overlapping exposures into
mosaics and for subtracting the stellar continuum from emission-line images.
This is handled automatically by ``MefPrep`` using the ``MAGZERO`` keyword
carried in each MEF header, but tools are provided to verify and, if
necessary, improve that calibration.

**Absolute photometry** — converting instrumental counts into physical flux
units (erg s⁻¹ cm⁻² Å⁻¹) — is needed for science-ready emission-line
images.  This requires comparison to reference stars with known spectra
(Gaia XP) and is handled by a separate set of tools.

.. toctree::
   :maxdepth: 1

   relative_photometry
   absolute_photometry


Where Photometry Enters the Pipeline
--------------------------------------

Relative photometry enters at two points:

1. **Mosaicking** (``MefPrep`` → ``Swarp``).  Each DECam exposure of the
   same field is taken on a different night and possibly under different
   atmospheric conditions.  The ``MAGZERO`` keyword encodes the throughput
   correction for each exposure.  ``MefPrep`` applies this to rescale every
   CCD image so that 1 DN = mag 28, placing all images on a common flux
   scale before SWarp combines them.  If ``MAGZERO`` is inconsistent across
   exposures, the mosaic will have visible seams between images taken on
   different nights.

2. **Continuum subtraction** (``CleanStars``).  Emission-line images are
   cleaned of stellar continuum by constructing a continuum image from
   broadband (r) and a narrowband off-line (N708) filter, then subtracting
   it from the emission-line image.  This subtraction relies on the stellar
   fluxes in the two images being on the same flux scale.  An error in the
   relative photometric calibration between the r-band and emission-line
   images will leave stellar residuals in the continuum-subtracted output.

Absolute photometry enters once the mosaics are ready:

3. **Flux calibration** (``ZeroPoint``).  Stars in the swarped tile images
   are cross-matched with Gaia XP spectra, which provide synthetic
   photometry at the exact wavelength of each filter.  The ratio of
   predicted to measured counts gives the conversion from DN to physical
   flux units.


Tool Overview
-------------

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Module
     - Used for
     - Purpose
   * - :doc:`MefPhot <api/MefPhot/index>`
     - Both
     - Forced aperture photometry on MEF images at Gaia or SMASH positions;
       primary source of ``TabPhot/`` files used by all calibration tools
   * - :doc:`CalcZeroPoint <api/CalcZeroPoint/index>`
     - Relative
     - Quick per-file ZP from sigma-clipped median; ``delta_zp`` shows
       how far each header ``MAGZERO`` is from the measured value
   * - :doc:`PhotEval <api/PhotEval/index>`
     - Relative
     - Multi-frame scatter check; ``chi2_nu_c`` measures whether
       overlapping frames are on the same flux scale
   * - :doc:`ZeroCalc <api/ZeroCalc/index>`
     - Relative / Absolute
     - Regression fit of ZP against Gaia or SMASH broadband magnitudes,
       with optional color term
   * - :doc:`PhotCompare <api/PhotCompare/index>`
     - Absolute
     - Aperture photometry on tile images with catalog cross-match
       and comparison plots
   * - :doc:`ZeroPoint <api/ZeroPoint/index>`
     - Absolute
     - Physical flux calibration (DN → erg s⁻¹ cm⁻² Å⁻¹) using
       Gaia XP spectra
   * - :doc:`GaiaCat <api/GaiaCat/index>`
     - Both
     - Retrieve and cache Gaia DR3 source positions and photometry
   * - :doc:`Smash <api/Smash/index>`
     - Relative
     - Retrieve SMASH DR2 catalog for Magellanic Cloud fields
   * - :doc:`PhotAnal <api/PhotAnal/index>`
     - Relative
     - Inter-filter consistency: compares the same star's DN across two
       filters (e.g. r and N662); reports scatter and color term that
       set the floor on CleanStars continuum subtraction quality
   * - :doc:`CheckPhot <api/CheckPhot/index>`
     - Relative
     - One-command wrapper: runs MefPhot (SMASH + Gaia), CalcZeroPoint,
       PhotEval, and PhotAnal in sequence; groups exposures by Image label
       from ``DeMCELS_images.txt``; writes ``CheckPhot/{field}_phot_check.fits``
       with INTRA and INTER table extensions and prints a formatted
       comparison summary
   * - :doc:`StarFind <api/StarFind/index>`
     - Absolute
     - Detect stars for PSF construction
   * - :doc:`PsfBuild <api/PsfBuild/index>`
     - Absolute
     - Build PSF models from selected stars
   * - :doc:`PsfPhot <api/PsfPhot/index>`
     - Absolute
     - PSF-fitting photometry for crowded fields
