MCLocate
========

.. py:module:: MCLocate

.. autoapi-nested-parse::

   Read an astropy table containg RA's and Dec's

   Space Telescope Science Institute

   Synopsis
   --------

   Read an astropy table containg RA's and Dec's
   and return the field and tile located cose
   to each object

   Command Line Usage
   ------------------

   ::

       usage: MCLocate.py [-h] [-out whatever] filename

       where -h prints this documentation and quits
       and  -out whatever prints the output to a file named whatever.
       (if this is missing then the output file will be X+filename

   Description
   -----------

   The routine is hardwired to look at MC_tiles.txt in
       the kred/config directory.

       The input table needs to have fields with columns
       named 'RA', and 'Dec' with the positions of objects
       in degrees.

       The routine finds the field and tile closest to
       each position and adds that (and the separation to the output
       table, which will be printed to the scrren,
       and written to a file.

   Primary Routines
   ----------------

   doit

   Notes
   -----

   The routine does not check that that the object is actually
       in a tile; it just reports the distance in arcmin to the
       closest tile.

   Version History
   ---------------

   240422 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   MCLocate.doit
   MCLocate.steer


Module Contents
---------------

.. py:function:: doit(filename='aida.txt', output='')

.. py:function:: steer(argv)

   Just a steering routine


