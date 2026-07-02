xMakeXcut
=========

.. py:module:: xMakeXcut

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

   xMakeXcut.do_xcut_plot
   xMakeXcut.doit
   xMakeXcut.get_pixel_scale
   xMakeXcut.make_xcut_tab
   xMakeXcut.radec2deg
   xMakeXcut.radec2pix
   xMakeXcut.round2integer
   xMakeXcut.steer


Module Contents
---------------

.. py:function:: do_xcut_plot(ztab_name, ra=10, dec=-70, xfilt='N662', exptime=400.0, size=100, ymin=30, ymax=100)

.. py:function:: doit(ra='01:02:38.6', dec='-71:59:08.23', xdir='.', xfilt='N662', exptime=400, size=200, ymin=30, ymax=100)

   Run the script


.. py:function:: get_pixel_scale(filename)

   Return the pixel scale from and image 
   that is presumed to be in extension 1


.. py:function:: make_xcut_tab(xdir)

.. py:function:: radec2deg(ra, dec)

   Convert coordiantes to degrees if
   necessary


.. py:function:: radec2pix(x, ra=15.0, dec=-71)

   Open a file and see if the position is inside the image.  If so incdicate the location and
   flux nearby.  Also collect some other information about the object

   If the position is not in the image, then -99, -99  is returned for the pixel position


.. py:function:: round2integer(arr)

.. py:function:: steer(argv)

   MakeXcut.py  -filt N662 -exptime 400 -size 200 -ymin 30 -ymax 100 xdir ra dec


