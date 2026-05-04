ZeroCalc
========

.. py:module:: ZeroCalc

.. autoapi-nested-parse::

   Calculate the magnitude zero point for an image

   Space Telescope Science Institute

   Synopsis
   --------

   Given one or more forced-photometry tables (output of MefPhot), fit a
   linear model to determine the magnitude zero point and optionally a color
   term, placing instrumental magnitudes on the Gaia or SMASH photometric
   scale.  The derived zero point can be compared directly to the MAGZERO
   keyword carried in the MEF file header.

   Command Line Usage
   ------------------

   ::

       ZeroCalc.py [-h] [-R] [-G] [-color] [-smash] [-out ROOT] file1.fits ...

       -h        Print this help and exit
       -R        Fit to the reference catalog R band (default)
       -G        Fit to the reference catalog G band
       -color    Add an independent color term to the fit
       -smash    Input tables were produced with SMASH as the reference catalog
       -out ROOT Output table root name (default: MagZero)

   Description
   -----------

   Two fit modes are available:

   **Simple fit** (default)::

       m_ref = m_inst + c_0

   **Color-corrected fit** (``-color``)::

       m_ref = m_inst + c_0 + c_1 * color

   where ``m_inst = 28 - 2.5 * log10(flux)`` and the color predictor is
   chosen to be independent of the target band:

   +--------+---------+-------+
   | Target | Catalog | Color |
   +========+=========+=======+
   | R      | Gaia    | G−R   |
   +--------+---------+-------+
   | G      | Gaia    | B−R   |
   +--------+---------+-------+
   | R      | SMASH   | G−R   |
   +--------+---------+-------+
   | G      | SMASH   | U−R   |
   +--------+---------+-------+

   The fitted ``c_0`` gives the zero-point correction: ``ZP_derived = 28 + c_0``.

   All fits are weighted by the per-source photometric uncertainty
   ``sigma_mag = 1.0857 * ErrNet / Net``, floored at 0.001 mag, so
   bright well-measured stars dominate the solution.  Sources with
   non-positive Net flux are excluded.

   Output filenames encode the band, catalog, and fit mode::

       MagZero.<band>.gaia.txt          simple fit, Gaia (default)
       MagZero.<band>.gaia.color.txt    color-corrected fit, Gaia
       MagZero.<band>.smash.txt         simple fit, SMASH
       MagZero.<band>.smash.color.txt   color-corrected fit, SMASH

   Each row in the summary table contains: Filter, Exptime, Root, MagZero
   (= 28 + c_0), c_0, c_1, rms, HdrZero (pipeline MAGZERO from the MEF
   header), Catalog, and Filename.

   Diagnostic plots are written to ``FigZero/`` with matching suffixes::

       FigZero/<band>_<root>.gaia.png
       FigZero/<band>_<root>.gaia.color.png
       FigZero/<band>_<root>.smash.png
       FigZero/<band>_<root>.smash.color.png

   Primary Routines
   ----------------

   do_one
       Process a single photometry table and return fit results.

   do_many
       Process multiple tables and accumulate results into a summary file.

   Notes
   -----

   This routine was developed to assess the consistency of MAGZERO as
   delivered by the DECam community pipeline.  Comparing the HdrZero column
   (pipeline MAGZERO) with the MagZero column (28 + c_0) across many
   exposures reveals systematic trends with filter, time, or CCD.  Running
   with both Gaia and SMASH provides an additional cross-check because the
   SMASH color term should be close to zero for r-band data.

   The simple fit (no color term) is the appropriate default for SMASH R-band
   calibration where the color term is expected to be near zero.  The
   color-corrected fit is useful for tightening scatter when a significant
   color term is present or for cross-checking Gaia G-band calibration.

   Version History
   ---------------

   251128 ksl
       Coding begun

   251228 ksl
       Updated to allow fitting to the Gaia G band

   260414 ksl
       Added SMASH support (-smash flag).
       Output filenames now include catalog suffix (.gaia / .smash).
       Plot axis labels and titles reflect the reference catalog used.

   260504 ksl
       Add simple fit mode (default) and optional color-corrected fit (-color).
       Color predictor is now always independent of the target band, eliminating
       the algebraic degeneracy where -G and -R gave identical c_0.
       Fits are now weighted by per-source photometric uncertainty.
       Output filenames include .color suffix when -color is used.


Functions
---------

.. autoapisummary::

   ZeroCalc.do_fig
   ZeroCalc.do_many
   ZeroCalc.do_one
   ZeroCalc.fit_magnitude_model
   ZeroCalc.get_filter_from_filename
   ZeroCalc.steer


Module Contents
---------------

.. py:function:: do_fig(xtab, band='R', outroot='', catalog='gaia', use_color=False)

   Plot diagnostic figures for a zero-point fit.

   The top two panels show instrumental vs. reference magnitudes before and
   after the fit.  The bottom two panels show residuals vs. reference
   magnitude, colour-coded by the colour index used in the fit.  Output is
   written to ``FigZero/`` with a filename that encodes the band, catalog,
   and fit mode.


.. py:function:: do_many(filenames, band='R', outroot='MagZero', catalog='gaia', use_color=False)

   Process multiple photometry tables and write a summary file.

   Parameters
   ----------
   filenames : list of str
       Paths to photometry FITS tables (output of MefPhot).
   band : str
       Reference band: ``'R'`` (default) or ``'G'``.
   outroot : str
       Root name for the output summary table. Default ``'MagZero'``.
   catalog : str
       Reference catalog: ``'gaia'`` (default) or ``'smash'``.
   use_color : bool
       If True, include a color term in the fit. Default False.


.. py:function:: do_one(filename, option='R', catalog='gaia', use_color=False)

   Process a single photometry table and return fit results.

   Parameters
   ----------
   filename : str
       Path to photometry FITS table (output of MefPhot).
   option : str
       Reference band to fit against: ``'R'`` (default) or ``'G'``.
   catalog : str
       Reference catalog: ``'gaia'`` (default) or ``'smash'``.
   use_color : bool
       If True, include an independent color term in the fit. Default False.
       Color predictor is chosen to be independent of the target band:
       Gaia R → G−R, Gaia G → B−R, SMASH R → G−R, SMASH G → U−R.

   Returns
   -------
   tuple
       ``(MagZero, c_0, c_1, rms, HdrZero, Filter, Exptime)``


.. py:function:: fit_magnitude_model(data, use_color=False)

   Fit catalogued magnitudes to a simple or color-corrected model.

   Simple (``use_color=False``)::

       m_ref = m_inst + c_0

   Color-corrected (``use_color=True``)::

       m_ref = m_inst + c_0 + c_1 * color

   Sources are weighted by ``sigma_mag = 1.0857 * ErrNet / Net``, floored
   at 0.001 mag.  Sources with non-positive Net flux are excluded.

   Parameters
   ----------
   data : astropy.table.Table
       Table with columns ``phot_mag``, ``target_mag``, ``Net``, ``ErrNet``,
       and (if use_color) ``target_color``.
   use_color : bool
       If True, fit a two-parameter model with a color term. Default False.

   Returns
   -------
   dict
       Keys: ``c_0``, ``c_0_err``, ``c_1``, ``c_1_err``, ``rms``, ``table``.
       ``c_1`` and ``c_1_err`` are 0.0 when ``use_color=False``.


.. py:function:: steer(argv)

   Command-line interface for ZeroCalc.

   usage: ``ZeroCalc.py [-h] [-R] [-G] [-color] [-smash] [-out ROOT] file1.fits ...``

