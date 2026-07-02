MefCheck
========

.. py:module:: MefCheck

.. autoapi-nested-parse::

   Check whether one or more fits files have headers

   Space Telescope Science Institute

   Synopsis
   --------

   Check whether one or more fits files have headers
   that conform to the fits standard, by looking at
   the warnings that astropy generates and for certain
   keywrods in the header.

   Command Line Usage
   ------------------

   ::

       usage: MefCheck.py [-np 8] [-all]  -dup_only Field ...

   Description
   -----------

   where

           -np 8 imples to run multiple threads, where each
               thread processes one field
           -all implies to do all fields in the DECam_MEF
               directory
           -dup_only just check for duplicate files in the MEF
               directory structure

   Primary Routines
   ----------------

   doit

   Notes
   -----

   Normally, this will not generate any errors.  The
       routine generates a table mef_qual.tab in the
       Summary directory for each of the fields.

   Version History
   ---------------

   230706 ksl
       Coding begun



Attributes
----------

.. autoapisummary::

   MefCheck.MEFDIR


Functions
---------

.. autoapisummary::

   MefCheck.check_for_duplicates
   MefCheck.check_mefs
   MefCheck.count_warnings
   MefCheck.do_one
   MefCheck.do_one_field
   MefCheck.steer
   MefCheck.test


Module Contents
---------------

.. py:data:: MEFDIR
   :value: 'DECam_MEF/'


.. py:function:: check_for_duplicates()

   Check to see if the directory structure is clean in
   the sense that there is only one version of the same
   file in the directory structure


.. py:function:: check_mefs(field='LMC_c45')

.. py:function:: count_warnings(filename)

.. py:function:: do_one(filename)

.. py:function:: do_one_field(field='LMC_c45')

   Simply get all the information one wants about the imaages associated
   with a single field


.. py:function:: steer(argv)

   This is just a steering routine, 
   and at present is not very useful,


.. py:function:: test(filename='c4d_190108_063145_ooi_N662_v1.fits.fz', outputfile='out.txt')

