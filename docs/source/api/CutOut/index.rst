CutOut
======

.. py:module:: CutOut

.. autoapi-nested-parse::

   Create cutout images of objects contained in a "masterfile"

   Space Telescope Science Institute

   Synopsis
   --------

   Create cutout images of objects contained in a "masterfile"
   in the images contained in a specific directory or set of directories

   Command Line Usage
   ------------------

   ::

       usage: CutOut.py [-r] [-size 10] image_directory object_table

   Description
   -----------

   where

           image_directory is the directory to be searched for fits images
           object_table is a list of objects with their positions

       and
           -h prints the documentatio and quits

           -r is an optional parameter that says not
           only to search the directory names above, but also
           all of the subdirectories

           -size 10 is an optional parameter that gives the size
           of the cutouts in arc minutes

   Primary Routines
   ----------------

   doit

   Notes
   -----

   The object_table should be an astropy table that contains (at least)
       the following columns:

       Source_name   - the name of the object (one word, with no spaces)
       RA
       Dec

   Version History
   ---------------

   230930 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   CutOut.extract_region
   CutOut.make_cutouts
   CutOut.steer


Module Contents
---------------

.. py:function:: extract_region(source_name, ra, dec, size_arcmin, input_fits, outdir='test', default_value=0, frac_off=0.01)

   Extract a region of a given size, but do not write an image if a fraction of an
   image has not data exceeds frac_off




.. py:function:: make_cutouts(object_table='lmc_snr.txt', size_arcmin=10, directory='T07', subdirectories=True, xoutdir='')

   where 
       object table is an ascii table containing the names and positions of 
       a specific object

       size_in_arcmin  is the size of the region to be extracted

       director is the directory to search

       subdirectories indicates whether just to search the top level directory, or to
       search all subdirectories as well



.. py:function:: steer(argv)

   Just a steering routine


