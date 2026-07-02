FindOverlaps
============

.. py:module:: FindOverlaps

.. autoapi-nested-parse::

   Calculate a which images obtained with a filter overlap

   Space Telescope Science Institute

   Synopsis
   --------

   Calculate a which images obtained with a filter overlap
   with which others

   Command Line Usage
   ------------------

   ::

       usage: FindOverlaps.py [-all] Field Tile1  Tile2 .....

       where -all implies that the ovelaps for all of the
       tiles in afiled will be dound

   Description
   -----------

   Find the images that ovelap that have the same
       filter.

       The routine calculates overlaps of the same filter,
       regardless of the exposure times.

       The overlaps are found based on the RA, DEC corrners
       of the immages.

       The routine outputs a single file that indicates
       which images overlap subject to thise conditions

       The routine DOES NOT calculate fluxes in the
       overlap regions, just that the ovelaps exist

   Primary Routines
   ----------------

   do_one_tile

   Notes
   -----

   It is important to recognize that care has
       been taken to avoid double listing, that is
       to say in separate lines of the output file
       that

       file_a matches file_b
       file_b matches file_a

       There is no significance to the order of
       file_a and file_b.

       This routine depends on having run either SetupTile or SetupSpecial, but does
       not require SwarpSetup, or Swarp.

   Version History
   ---------------

   230524 ksl
       Coding begun

   250128 ksl
       Modified the code so that both exposure times are recorded. Choosing

       what exposure times to match background to occurs later.



Functions
---------

.. autoapisummary::

   FindOverlaps.do_one_filter
   FindOverlaps.do_one_tile
   FindOverlaps.get_files
   FindOverlaps.get_overlap
   FindOverlaps.steer


Module Contents
---------------

.. py:function:: do_one_filter(field='LMC_c45', tile='T07', xfilter='Ha')

   Find the overlaps for one filter in one tile


.. py:function:: do_one_tile(field='LMC_c45', tile='T07')

   Find the overlaps for all of the filters in one tile


.. py:function:: get_files(field, tile, xfilter='Ha')

   Idenfify all of the files associated
   with a specific filter in a tile

   Note: 
   This reads the Field_Tile.txt files in the
   Summary directory


.. py:function:: get_overlap(row, xxtab)

   Get all of the images which have the
   same exposure time as the given row
   and also overlap with a given row 
   in the table

   Note eliminating images that have
   different exposure times implies
   that we do not intend to match backgrounds
   from images that have differnt times


.. py:function:: steer(argv)

   This is just a steering routine for running swarp on one or more
   tiles from the command line



