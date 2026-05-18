ZeroCalc
========

.. py:module:: ZeroCalc

.. autoapi-nested-parse::

   Calculate the magnitude zero point for an image

   Space Telescope Science Institute

   Synopsis
   --------

   Given one or more forced-photometry tables (output of MefPhot), fit a
   linear model to determine the magnitude zero point and color term that
   place instrumental magnitudes on the Gaia or SMASH photometric scale.
   The derived zero point can be compared directly to the MAGZERO keyword
   carried in the MEF file header.

   Command Line Usage
   ------------------

   ::

       ZeroCalc.py [-h] [-R] [-G] [-color] [-smash] [-fig] [-np N] [-out ROOT] file1.fits file2.fits ...

       -h        Print this help and exit
       -R        Fit to the reference catalog R band (default)
       -G        Fit to the reference catalog G band
       -color    Add a color term to the fit (default: simple zero-point fit only)
       -smash    Input tables were produced with SMASH as the reference catalog
       -fig      Write diagnostic plots to FigZero/ (default: no plots)
       -np N     Number of parallel worker processes (default: 1)
       -out ROOT Output table root name (default: MagZero)

   Description
   -----------

   Two fit modes are available:

   **Simple fit** (default)::

       m_ref = m_inst + c_0

   **Color-corrected fit** (``-color``)::

       m_ref = m_inst + c_0 + c_1 * color

   where the color predictor is chosen to be independent of the target band:

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

   In all cases ``m_inst = 28 - 2.5 * log10(flux)`` and the fitted ``c_0``
   gives the zero-point correction: ``ZP_derived = 28 + c_0``.

   Output filenames encode the band, catalog, fit mode, and run date::

       MagZero.<band>.gaia.<date>.txt          simple fit, Gaia (default)
       MagZero.<band>.gaia.color.<date>.txt    color-corrected fit, Gaia
       MagZero.<band>.smash.<date>.txt         simple fit, SMASH
       MagZero.<band>.smash.color.<date>.txt   color-corrected fit, SMASH

   where ``<date>`` is a six-digit YYMMDD string (e.g. ``260518``).  Runs on
   the same day append to the same file; runs on different days produce
   separate files.

   Each row in the summary table contains: Filter, Exptime, Root, MagZero
   (= 28 + c_0), c_0, c_1, rms, HdrZero (pipeline MAGZERO from the MEF
   header), Catalog, and Filename.

   Diagnostic plots are written to ``FigZero/`` only when ``-fig`` is given::

       FigZero/<band>_<root>.gaia.png
       FigZero/<band>_<root>.gaia.color.png
       FigZero/<band>_<root>.smash.png
       FigZero/<band>_<root>.smash.color.png

   Figure generation is the dominant per-file cost; omitting ``-fig`` when
   processing large batches gives a significant speedup.  The remaining
   bottleneck is FITS read I/O; ``fitsio`` with selective column reads would
   be the next step if further speedup is needed.

   Primary Routines
   ----------------

   do_one
       Process a single photometry table and return fit results.

   do_many
       Process multiple tables in parallel and accumulate results into a
       summary file.

   Notes
   -----

   This routine was developed to assess the consistency of MAGZERO as
   delivered by the DECam community pipeline.  Comparing the HdrZero column
   (pipeline MAGZERO) with the MagZero column (28 + c_0) across many
   exposures reveals systematic trends with filter, time, or CCD.  Running
   with both Gaia and SMASH provides an additional cross-check because the
   SMASH color term should be close to zero for r-band data.

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
       Color predictor is now always independent of the target band:
       Gaia R: G-R, Gaia G: B-R, SMASH R: G-R, SMASH G: U-R.
       Output filenames include .color suffix when -color is used.

   260518 ksl
       Performance improvements and usability changes.
       Add -np N flag for parallel processing via multiprocessing.Pool.
       Add rasterized=True to all scatter calls in do_fig (no quality loss
       for PNG output, faster rendering).
       Set matplotlib Agg backend explicitly so worker processes need no display.
       Make figure generation opt-in with -fig flag (previously always-on);
       figure generation was found to be the main per-file bottleneck when
       processing large batches. Without -fig, throughput is limited by FITS
       read I/O rather than CPU; fitsio with selective column reads would be
       the next step if further speedup is needed.
       Progress reporting now prints index/total and filename per file.
       Output filenames now include a date suffix (e.g. MagZero.R.gaia.260518.txt)
       to avoid overwriting previous runs.
       Fixed latent bug in do_many where a failed file would cause a
       column-length mismatch when building the output Table.



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

   Plot results.  The top two panels plot the
   magnitudes as measured by aperstats, assuming
   a zeropoint of 28

   The bottom two panels plot the fitted
   fluxes


.. py:function:: do_many(filenames, band='G', outroot='MagZero', catalog='gaia', use_color=False, n_processes=1, use_fig=False)

   Process multiple photometry tables and accumulate results into a summary
   file.  When ``n_processes`` > 1 the files are distributed across a
   ``multiprocessing.Pool``.

.. py:function:: do_one(filename='TabPhot/c4d_241122_023910_ooi_N673_v1.fits', option='R', catalog='gaia', use_color=False, use_fig=False)

   Process a single photometry table and return fit results.

   Parameters
   ----------
   filename : str
       Path to photometry FITS table (output of MefPhot).
   option : str
       Reference band to fit against: 'R' (default) or 'G'.
   catalog : str
       Reference catalog: 'gaia' (default) or 'smash'.
   use_color : bool
       If True, include an independent color term in the fit. Default False.
       Color predictor is chosen independent of the target band:
       Gaia R: G-R, Gaia G: B-R, SMASH R: G-R, SMASH G: U-R.
   use_fig : bool
       If True, write a diagnostic plot to FigZero/. Default False.


.. py:function:: fit_magnitude_model(data, use_color=False)

   Fit catalogued magnitudes to either a simple or color-corrected model.

   Simple (use_color=False):   m_ref = m_inst + c_0
   Color-corrected (use_color=True): m_ref = m_inst + c_0 + c_1 * color

   Sources are weighted by their photometric uncertainty:
       sigma_mag = 1.0857 * ErrNet / Net
   floored at 0.001 mag so a handful of very bright stars do not dominate.
   Sources with non-positive Net flux are excluded.

   Parameters
   ----------
   data : astropy.table.Table
       Table with columns 'phot_mag', 'target_mag', 'Net', 'ErrNet', and
       (if use_color) 'target_color'.
   use_color : bool
       If True, include a color term in the fit. Default False.

   Returns
   -------
   dict with keys c_0, c_0_err, c_1, c_1_err, rms, table.
   c_1 and c_1_err are 0.0 when use_color=False.


.. py:function:: get_filter_from_filename(filename)

.. py:function:: steer(argv)

   usage: ZeroCalc.py [-h] [-R] [-G] [-color] [-smash] [-fig] [-np N] [-out ROOT] file1.fits file2.fits ...
