ZeroPoint
=========

.. py:module:: ZeroPoint

.. autoapi-nested-parse::

   ZeroPoint - Compute flux zero points using Gaia XP spectra

   Space Telescope Science Institute

   Synopsis
   --------

   Estimate the flux corresponding to 1 DN (data number) in the H-alpha
   (6563 Å) and [SII] (6720 Å) filters by cross-referencing aperture
   photometry from PhotCompare with Gaia DR3 XP continuous spectra.

   Command Line Usage
   ------------------

   ::

       ZeroPoint.py [-h] xmatch_file [xmatch_file ...]

   **Required Arguments:**

   xmatch_file
       One or more cross-match files (``TabPhot/xmatch_*.txt``) produced
       by ``PhotCompare.py``.

   **Optional Arguments:**

   -h
       Display help message and exit

   Description
   -----------

   For each xmatch file this routine:

   1. Selects stars in a temperature window (8000–15000 K by default) and
      magnitude range (14 < R < 18).  If the table contains an
      ``xp_spec_exists`` column (written by the Gaia catalog builder), only
      stars with confirmed XP spectra are kept, avoiding wasted retrieval
      attempts.

   2. Retrieves Gaia DR3 XP continuous spectra for those stars via
      ``get_gaia_spec()``, which first checks a local per-star cache
      (``GaiaSpec/``) and falls back to the Gaia archive on a cache miss.

   3. Interpolates the calibrated spectrum at 6563 Å (H-alpha) and 6720 Å
      ([SII]) to obtain the expected physical flux in erg s-1 cm-2 Å-1.

   4. Divides each physical flux by the measured net counts (``Net`` column)
      to derive the conversion factor.

   Primary Routines
   ----------------

   do_one, gaia_choose, get_gaia_flux

   Output
   ------

   * ``PhotMaster.txt`` — running summary table of median, mean, and
     standard-deviation conversion factors for each xmatch file, for both
     H-alpha and [SII].  New results are appended (or replaced if the same
     filename already appears in the existing table).

   * ``Figs_phot/<xmatch_root>.png`` — six-panel diagnostic plot: flux
     distributions and flux-vs-net-counts scatter for both filters.

   Notes
   -----

   * Input xmatch files are produced by ``PhotCompare.py`` and live in
     ``TabPhot/``.

   * Spectrum retrieval requires the ``gaiaxpy`` package and valid Gaia
     archive credentials in ``~/.gaia_credentials``.  Stars without XP
     spectra (G > 19 or otherwise absent) are silently skipped.

   * If the ``xp_spec_exists`` column is absent from the xmatch file, the
     routine still works but will attempt retrieval for every star in the
     temperature/magnitude window, which is slower.

   * This routine was used to estimate sensitivity in Points+24.

   Version History
   ---------------

   240505 ksl
       Coding begun

   260429 ksl
       Filter on ``xp_spec_exists`` column to skip stars without confirmed
       spectra.  Unpack ``(wave, flux)`` tuple from ``get_gaia_spec`` into
       an astropy Table before interpolation.  Suppress gaiaxpy boilerplate
       output during batch retrieval.



Functions
---------

.. autoapisummary::

   ZeroPoint.do_one
   ZeroPoint.do_plot
   ZeroPoint.gaia_choose
   ZeroPoint.get_flux
   ZeroPoint.get_gaia_flux
   ZeroPoint.steer
   ZeroPoint.tab_remove_bad


Module Contents
---------------

.. py:function:: do_one(xmatch_file, Rmin=14, Rmax=18, tmin=8000, tmax=15000, nmax=1000)

.. py:function:: do_plot(ztab, outroot='foo')

.. py:function:: gaia_choose(filename='TabPhot/xmatch_Gaia.LMC_c32_T07_LMC_c32_T07.N673.t800_phot.txt', Rmin=14, Rmax=18, tmin=8000, tmax=15000, nmax=1000)

.. py:function:: get_flux(xtab, wavelength=6563)

   Get the flux of a star at a particular wavelength



.. py:function:: get_gaia_flux(xtab, key='Ha', wave=6563)

.. py:function:: steer(argv)

.. py:function:: tab_remove_bad(qtab, key='teff')

   Return a table with masked values of a column removed


