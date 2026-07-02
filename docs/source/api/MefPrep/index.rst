MefPrep
=======

.. py:module:: MefPrep

.. autoapi-nested-parse::

   Prepare files for combining with Swarp

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       MefPrep.py [-h] [-all] [-finish] [-back_none] [-back_min] [-np N]
                  [-zp COLUMN] FIELD [FIELD ...]

   Rescales individual DECam CCD images to a common flux scale (1 DN = mag 28)
   and optionally subtracts an initial background estimate.  Output files are
   written to ``DECam_PREP/{field}/data/``.

   Prerequisites: MEF files downloaded into the standard directory structure and
   ``MefSum.py`` run to create the ``Summary/{field}_mef.tab`` inventory tables.

   Optional Arguments
   ------------------

   -h
       Print this help and exit.

   -all
       Process all fields found under ``DECam_MEF/`` (interactive confirmation
       required).

   -finish
       Skip output files that already exist (do not reprocess).

   -back_none
       Do not subtract any background estimate.

   -back_min
       Subtract the minimum CCD background across all detectors in each exposure
       (instead of the default mode estimate).

   -np N
       Number of parallel worker processes (default: 1).

   -zp COLUMN
       Use empirical zero points from the named column in
       ``Summary/{field}_mef.tab`` (e.g. ``-zp ZP_smash_r``) instead of the
       ``MAGZERO`` header keyword.  The value must satisfy ``20 < ZP < 35``; if
       the column is absent or the value is invalid, the script falls back to the
       header ``MAGZERO`` with a warning.  Generates ``ZP_USE`` and
       ``ZP_SRC`` keywords in every output header.

   Output Header Keywords
   ----------------------

   * ``ZP_USE`` -- zero point actually applied for flux scaling
   * ``ZP_SRC`` -- source of that zero point (column name or ``MAGZERO``)
   * ``XFACTOR``     -- multiplicative scale factor applied to pixel data
   * ``SUB_BACK``    -- ``TRUE`` / ``FALSE`` depending on whether background was subtracted
   * ``GAIN``, ``SATURAT``, ``RDNOISE`` -- rescaled detector properties

   Examples
   --------

   Standard single-field run with background subtraction (mode)::

       MefPrep.py -np 8 LMC_c42

   Re-run using empirical SMASH r-band zero points::

       MefPrep.py -np 8 -zp ZP_smash_r LMC_c42

   Skip already-processed files, no background subtraction::

       MefPrep.py -finish -back_none LMC_c42

   Notes
   -----

   Unlike ``MefSum.py``, individual files are processed in separate threads so
   all requested cores are typically utilised.  If too many files are open
   simultaneously, reduce ``-np``.

   Version History
   ---------------

   230513 ksl
       Coding begun

   230621 ksl
       Revised so that the default is to subtract background

   260608 ksl
       Add -zp COLUMN to use empirical zero points from Summary/_mef.tab instead
       of the header MAGZERO.  ZP_USE and ZP_SRC written to output headers.
       Fixed pre-existing bug where computed back value was not passed to prep_one_det.



Attributes
----------

.. autoapisummary::

   MefPrep.CCDDIR
   MefPrep.MEFDIR
   MefPrep.SUMDIR


Functions
---------

.. autoapisummary::

   MefPrep.get_no_jobs
   MefPrep.prep_one_det
   MefPrep.prep_one_field
   MefPrep.prep_one_mef
   MefPrep.steer
   MefPrep.xprep_one_field


Module Contents
---------------

.. py:data:: CCDDIR
   :value: 'DECam_CCD/'


.. py:data:: MEFDIR
   :value: 'DECam_MEF/'


.. py:data:: SUMDIR
   :value: 'Summary/'


.. py:function:: get_no_jobs(jobs)

   Check how many jobs are running


.. py:function:: prep_one_det(filename='DECam_MEF/LMC_c42/mef/c4d_211111_024404_ooi_N673_v1.fits.fz', ext='S2', prepdir='DECam_PREP/LMC_c42/T07/', back=0, redo=True, magzero=None, zp_col='')

   Prepare one file for combining with swarp.  This version scales by MAGZERO
   where 1DN will be a flux corresponding to mag 28
   History:

   230504 - Removed .fz from output file names.  Note that this does not mean that the
   data is not compressed.   Added extra keywords to primary header of output file
   that may be useful for swarp, including the rescaled saturation level

   230611 - Addapted from prep_one_file.  This version allows for several backgroud
   options, based on information assembled by MEFSum

   :


.. py:function:: prep_one_field(field='LMC_c42', back_type='none', redo=False, outdir='', zp_col='')

   Prep a field with a single processor

   This routine identifies the mefs to be processed when one
   is in single processor mode


.. py:function:: prep_one_mef(field='LMC_c42', root='c4d_190109_061931_ooi_N662_v1', back_type='min', redo=False, outdir='', zp_col='')

   Split the mef files into their individual extenstions and prepare them for the downstream parts of the processing

   This version stores all of the actual output files in a single directory for each field. Setup_Tiles (will be modified
   for the new file structure


   Note that background subtraction here, if is carriout out, subtractis the same value from all of 
   the CCDs in a single image.  

   What is done is slightly confusing.  MEFSum.py has estimated a different background in each of
   the CCDs in the exposure and stored this in the det file.  Here we take all of those measurements
   for the mef image and either take the median of those values or the minimum.




.. py:function:: steer(argv)

   This is just a steering routine


.. py:function:: xprep_one_field(field='LMC_c42', back_type='none', redo=False, outdir='', nproc=4, zp_col='')

   Process the mef files in a field usingwith using more than one core.

   This routine creates a bundle of files to process 
   at once. This is necessary because multiprocessing
   opens many files (note that in linux everything is
   a file, so this does not refer to fits files), and
   this can exceed the file limit.  The sympton of this
   is a message that says too many files are open.  If
   one sees this message, one may need to reduce n
   below.



