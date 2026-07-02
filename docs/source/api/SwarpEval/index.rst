SwarpEval
=========

.. py:module:: SwarpEval

.. autoapi-nested-parse::

   Create plots of the images that are in the DECam_Swarp  after running Swarp.py, CleanStars.py

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   directory.



Attributes
----------

.. autoapisummary::

   SwarpEval.exec_dir


Functions
---------

.. autoapisummary::

   SwarpEval.display_fits_image
   SwarpEval.get_images
   SwarpEval.get_stats
   SwarpEval.make_plots
   SwarpEval.steer


Module Contents
---------------

.. py:function:: display_fits_image(image_file, scale='linear', invert=False, vmin=None, vmax=None, outfile='')

.. py:data:: exec_dir

.. py:function:: get_images(xdir='DECam_SWARP2', field='LMC_c01', tile='T01')

.. py:function:: get_stats(xfiles)

   Characterize the stats of the images, assuing
   that 0 represents a pixel should be masked.


.. py:function:: make_plots(xdirectory='DECam_SWARP2', field='LMC_c01', tile='T01')

.. py:function:: steer(argv)

   This is just a steering routine for running swarp on one or more
   tiles from the command line



