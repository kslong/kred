MefSum
======

.. py:module:: MefSum

.. autoapi-nested-parse::

   Sumarize mef files that constitute the raw data from

   Space Telescope Science Institute

   Synopsis
   --------

   Sumarize mef files that constitute the raw data from
   DECam observations of the LMC and SMC

   This routines assume a standard setup with the a
   local link to the raw DECam mef files (or part of
   them with the same format as on Box

   Command Line Usage
   ------------------

   ::

       usage: MefSum.py [-h] [-all] [-np 3]  [-det] [-mef] field1 field2 ...

       where
       -h prints the documentation
       -np 3  causes the processing to be carried out with a a given no of threads
       -all causes all files in the MEF directories to be processed. This should
       only be used with caution since it will take considerable time, and so
       the user is asked to confirm this option.
       -det just runs the individual ccd portion of the process (this is diagnostic)
       -mef just runs the overall mef portions (this is diagnostic)
       -sigma_clipped (Causes some statitics to be calculaed using sigma clipping).
       This is a diagnostic mode and is not normally used.  The estimates
       are now based on the mode

       and fields comprise one or more fields e.g LMC_c42  to be processed

   Description
   -----------

   The routine reads the MEF files and summarizes information from
       the headers in two tables, a mef.tab file and a det.tab file
       The routine also calculates some statstics associated with the
       counts in each detector, that are used to carry out an
       intial background subtraction.

   Primary Routines
   ----------------

   Notes:

       At present this does not check that a field has been processed

       Here multiple processors are used only to handle individual files
       so it makes sense to us multiple threads even if only one
       field is being processed.,

       Note that this does not read the  file produced by MefCheck, which is
       really u

   Notes
   -----

   At present this does not check that a field has been processed

       Here multiple processors are used only to handle individual files
       so it makes sense to us multiple threads even if only one
       field is being processed.,

       Note that this does not read the  file produced by MefCheck, which is
       really u

   Version History
   ---------------

   230513 ksl
       Coding begun



Attributes
----------

.. autoapisummary::

   MefSum.MEFDIR


Functions
---------

.. autoapisummary::

   MefSum.check_exposures
   MefSum.do_one
   MefSum.do_one_det_overview
   MefSum.get_det_overview
   MefSum.get_keyword
   MefSum.get_mef_overview
   MefSum.halfsamplemode
   MefSum.steer


Module Contents
---------------

.. py:data:: MEFDIR
   :value: 'DECam_MEF/'


.. py:function:: check_exposures(field, instrument_file='DeMCELS_images.txt')

   Check that the exposures are contained in the instrument file


.. py:function:: do_one(field='LMC_c45', nproc=1, do_mef=True, do_det=True, calc_sigma_clipped=False)

   Proecess a single field.                                              


.. py:function:: do_one_det_overview(filename, keys, calc_sigma_clipped=False)

   Get stastictics and keyword values for the various ccd images in a single mef image    


.. py:function:: get_det_overview(field='LMC_c45', nproc=1, calc_sigma_clipped=False)

   Get information about the individual CCDs, including 
   some statistics associated with the data that is 
   useful for estimating background

   Normally this produces the standard median, mean and standard deviation
   for the image, but asa diagnosic mode one can produced sigma clipped
   statistics.   The mode is also produced, and usually this is what
   is used for initial sky subtraction.


.. py:function:: get_keyword(key, ext)

   where key is just the name of a keyword
   and ext points to a specific extension
   of an already opened fits file, e.g x[0]

   Note that one does not include .header in
   the call

   History:

   230520 - Fix to handle the fact that there
   are some files for which the Band is not filled
   in.



.. py:function:: get_mef_overview(field='LMC_c45')

   Get keyword  information from the 0th extension of the mef files
   in a single field


.. py:function:: halfsamplemode(inputData, axis=None)

   Compute mode of an array using the half-sample algorithm


.. py:function:: steer(argv)

   This is just a steering routine, 


