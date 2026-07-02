MakeXcut
========

.. py:module:: MakeXcut

.. autoapi-nested-parse::

   This is a diagnositc routine intended to allow one to

   Space Telescope Science Institute

   Synopsis
   --------

   This is a diagnositc routine intended to allow one to
   make plots of rows and columns in regions of images
   that show breaks in the apparent flux

   Command Line Usage
   ------------------

   ::

       usage: MakeXcut.py  -filt N662 -exptime 400 -size 200 -ymin 30 -ymax 100 xdir ra dec

       where

       -xfilt designates a fileter to use
       -exptime designates and expousre time to chose
       -size indicates the size in pixels of the xcuts
       -ymin and -ymax set the min and maxim y leve

       and
       xdir is the directory in which to look for fits files
       ra and dec are the ra and dec of a position to inspect

   Description
   -----------

   Primary routines:

       doit

   Primary Routines
   ----------------

   doit

   Notes
   -----

   History:

   231018 ksl Coding begun

   Version History
   ---------------

   231018 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   MakeXcut.convert_coordinates
   MakeXcut.do_xcut_plot
   MakeXcut.doit
   MakeXcut.find_overlaps
   MakeXcut.radec2pix
   MakeXcut.round2integer
   MakeXcut.steer


Module Contents
---------------

.. py:function:: convert_coordinates(ra, dec, is_degrees=True)

   Convert coordiantes to degrees if
   necessary


.. py:function:: do_xcut_plot(ztab, xfilt='N662', exptime=400.0, size=100, ymin=30, ymax=100)

.. py:function:: doit(ra='01:02:38.6', dec='-71:59:08.23', xdir='.', xfilt='N662', exptime=400, size=200, ymin=30, ymax=100)

   Run the script


.. py:function:: find_overlaps(xdir='.', ra='01:02:38.6', dec='-71:59:08.23', size=1000)

.. py:function:: radec2pix(filename='c4d_221221_004943_ooi_N673_req_N24.fits', ra='01:02:38.6', dec='-71:59:08.23')

   Open a file and see if the postion is inside the image.  If so incdicate the location and
   flux nearby.  Also collect some other information about the object


.. py:function:: round2integer(arr)

.. py:function:: steer(argv)

   MakeXcut.py  -filt N662 -exptime 400 -size 200 -ymin 30 -ymax 100 xdir ra dec


