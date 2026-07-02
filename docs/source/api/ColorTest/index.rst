ColorTest
=========

.. py:module:: ColorTest

.. autoapi-nested-parse::

   See how well images have been subtracted

   Space Telescope Science Institute

   Synopsis
   --------

   See how well images have been subtracted

   Command Line Usage
   ------------------

   ::

       Usage ColorTest.py -h -dir dirname file1..

       where:
       -h prints out this help and exits
       -dir searches the directory dirname and
       all subdirectories for continuum subtracted
       emission line files, processes these images
       file1 ... if -xdir is not specified carries out
       the analysis on specific files

   Description
   -----------

   The routine access the GAIA database to find stars in the
       field.  It then carries out forced photometry on emission
       line image from which the stars were subtracted and on
       the continuum image that was used for subtraction.  It
       produces a plot that shows the amount of over-subtaction
       or under-subtraction

   Primary Routines
   ----------------

   doit

   Notes
   -----

   The routine looks for files with specific names to decide
       what images are are those that can be analyzed.  If new
       filters are involved some changes are requried.

   Version History
   ---------------

   251108 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   ColorTest.color_compare
   ColorTest.do_dir
   ColorTest.doit
   ColorTest.find_cont
   ColorTest.get_base_name
   ColorTest.plot_both
   ColorTest.steer


Module Contents
---------------

.. py:function:: color_compare(cont_image, subtracted_image, forced=False)

.. py:function:: do_dir(xdir='DECam_SUB2/LMC_c37/T07', nrow_max=-1)

   Sort out how to call doit based on the directories


.. py:function:: doit(continuum_file='DECam_SWARP2/LMC_c35/T06/LMC_c35_T06.r.fits', subtracted_file='DECam_SUB2/LMC_c35/T06/LMC_c35_T06.ha_sub_r.fits')

.. py:function:: find_cont(filenames)

   Locate the companion files for subtracted files


.. py:function:: get_base_name(filename='DECam_SWARP2/LMC_c35/T06/LMC_c35_T06.r.fits')

.. py:function:: plot_both(gaia, final, title='')

.. py:function:: steer(argv)

   Run the script given choices from the command line

   Usage ColorTest.py -h -xdir whatever file1..



