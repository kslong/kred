CleanStars
==========

.. py:module:: CleanStars

.. autoapi-nested-parse::

   Produce pure emision (and continuum) images via

   Space Telescope Science Institute

   Synopsis
   --------

   Produce pure emision (and continuum) images via
   simple subtraction

   Command Line Usage
   ------------------

   ::

       usage: CleanStars.py [-all] [-bsub] field [tiles]

       where
       -all will cause CleanStars to be run on all 16 tiles in a field
       and if -all is not given, then one or more tiles should be listed
       -bsub will search for inputs in the DECam_SWARP2 directories which
       have "better" background matching, while if it is absence
       one will use the data in the DECam_SWARP directories, which use
       the background form the overal fields

   Description
   -----------

   This routine produces simple images that
       have been subtracted to remove the coinuum
       (or from the r-band continuum the emissionlines)

       It must be run from the top level kred directory

       The routine is currently hardwired to use the
       longest exposures of a particular filter type.


       The routine is NOT sophisticated.

       For creating the r-band results, we first subtract the two emission
       line images from the r-band images; this should produce an image
       with the emission lines removed.  We then subtract this from the emission
       line images.  There is some scaling involved.

       For creating the N708 results, we do a simple subtraction.

       (Note that before we actually do the subtractions, we calculated
       a background level n the continuum images.   This is intended
       to leave whatever backgrouground level in the emission line images
       unaffected.)

       The files that are produced are nominally the following (for LMC_c42_T07):
       LMC_c42_T07.ha_sub_r.fits     - pure ha image using r-band for continuum
       LMC_c42_T07.s2_sub_r.fits     - pure s2 image using r-rand for continuum
       LMC_c42_T07.r_sub.fits        - pure r-band image after emission lines are subtracted
       LMC_c42_T07.ha_s2_sub.fits    - the emission line portion of the r-band image

       LMC_c42_T07.ha_sub_N708.fits  - pure ha image based on subtrcting the n708 image
       LMC_c42_T07.s2_sub_N708.fits  - pure ha image based on subtrcting the n708 image
       LMC_c42_T07.n708_sub.fits     - n708 image (after subtracting a biased median)

       The outputs are written to a subdirectory of DECam_SUB

   Primary Routines
   ----------------

   doit
       make_rband_subtractions
       make_n708_subtractions
       make_n540_subtractions

   Notes
   -----

   History:

   230717 ksl Coding begun

   Version History
   ---------------

   230717 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   CleanStars.doit
   CleanStars.fits_deep_copy
   CleanStars.make_n540_subtractions
   CleanStars.make_n708_subtractions
   CleanStars.make_rband_subtractions
   CleanStars.steer


Module Contents
---------------

.. py:function:: doit(xdir='data', outdir='data')

   Routine to produce star subtracted images


.. py:function:: fits_deep_copy(orig)

   Make a deep copy of a fits image

   Astropy does not have a way, apparently to make 
   a deep copy of an image; instead it makes shallow
   copies, which means that often one does not know
   what one is dealing with. 

   The basic problem here is that if one does not
   make deep copies, one will not necessary know
   what image you are working with, one that is
   pristine or one that has been modified

   This little routine
   avoids this issue.




.. py:function:: make_n540_subtractions(o3='data/LMC_c42_T07.N501.t800.fits', n540='data/LMC_c42_T07.N540.t300.fits', outroot='test3')

   For subtraction of N540 from N501 and assume no emission line 
   contamination, and since everything is scaled to the same level
   we just do a straight subtraction.


.. py:function:: make_n708_subtractions(ha='data/LMC_c42_T07.N662.t800.fits', s2='data/LMC_c42_T07.N673.t800.fits', n708='data/LMC_c42_T07.N708.t400.fits', outroot='test2')

   For narrow band subtraction we assume there is no emission line
   contamination, and since everything is scaled to the same level
   we just do a straight subtractin


.. py:function:: make_rband_subtractions(ha='data/LMC_c42_T07.N662.t800.fits', s2='data/LMC_c42_T07.N673.t800.fits', r='data/LMC_c42_T07.r.t060.fits', outroot='test')

.. py:function:: steer(argv)

   This is just a steering routine for running swarp on one or more
   tiles from the command line



