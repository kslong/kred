BackSub
=======

.. py:module:: BackSub

.. autoapi-nested-parse::

   Use calculated backgrounds and subtract these

   Space Telescope Science Institute

   Synopsis
   --------

   Use calculated backgrounds and subtract these
   from the existing images

   Command Line Usage
   ------------------

   ::

       usage: BackSub.py [-h] [-all] -np 8  field T01 T02 ...

       where

       * -h gives this documentation
       - all causes all tiles in the field to be processed
       - np 8  causes 8 procesors to be used

   Description
   -----------

   The routine reads the files in the DECam_Prep/field/tile
       directory and writes the subtacted in imaages to
       DECam_Prep2/field/tile

   Primary Routines
   ----------------

   Notes:
       The input tables  for this have names like  Summary/Field_Tile_bbb.txt

   Notes
   -----

   The input tables  for this have names like  Summary/Field_Tile_bbb.txt

   Version History
   ---------------

   230619 ksl
       Coding begun



Attributes
----------

.. autoapisummary::

   BackSub.CWD
   BackSub.DATADIR
   BackSub.PREPDIR


Functions
---------

.. autoapisummary::

   BackSub.do_one_file
   BackSub.steer
   BackSub.subtract_one_tile


Module Contents
---------------

.. py:data:: CWD

.. py:data:: DATADIR

.. py:data:: PREPDIR

.. py:function:: do_one_file(filein, fileout, b)

   This simply subtracts a background from an image

   250201 - For reasons that are unclear the current version of
   astropy gave an error associated with PCOUNT at this
   stage.  The value for an image should be 0, but was
   a large number.  So I have forced a fix to this.


.. py:function:: steer(argv)

   This is just a steering routine


.. py:function:: subtract_one_tile(field='LMC_c42', tile='T07', np=1)

