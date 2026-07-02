BackCalc
========

.. py:module:: BackCalc

.. autoapi-nested-parse::

   Calculate backgrounds for images for which

   Space Telescope Science Institute

   Synopsis
   --------

   Calculate backgrounds for images for which
   the statistics have been gathered

   Command Line Usage
   ------------------

   ::

       usage: BackCalc.py [-all]   Field T01 ...

   Description
   -----------

   Uset a variety of techingiques to try to find
       the best backgrounds to subtract from the
       fileds to match backgrounds

       The main output is a file field_tile_bbb.txt
       that contains the calculated values for
       the offsets.

       Files (xxx) are also generated that show how
       good the the actual fits is.

   Primary Routines
   ----------------

   Notes:

   Notes
   -----

   History:

   230609 ksl Coding begun

   Version History
   ---------------

   230609 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   BackCalc.create_inputs
   BackCalc.diff_eval
   BackCalc.do_one_tile
   BackCalc.do_svd
   BackCalc.monte
   BackCalc.plot_fit_results
   BackCalc.steer
   BackCalc.svdskyfit
   BackCalc.xeval
   BackCalc.xpath


Module Contents
---------------

.. py:function:: create_inputs(infile='data/LMC_c45_T07_xxx.txt', ximage='N673')

   Create the inputs needed to fit different backgrounds for a single filter and
   exposure time.


.. py:function:: diff_eval(xdelta, xback, return_diff=False)

   Given a table that contains the measured differences between 
   the flux (xdelta) in different images, and another the table (xback) that 
   gives the current version the backgrounds for the images,
   calculated the 'effective chi*2'

   Note that this is the only statistic we have access to
   for a real case, since we do not know what the actual offsets
   are.

   The goodness of fit is returned 


.. py:function:: do_one_tile(field, tile)

   Generate a set of backgrounds to be used for matching images
   in a tile based on differences in the fluxes in images..


.. py:function:: do_svd(xcross, bfile, threshold=0.1, new=True, verbose=False)

   The driving routine for a SVD fit of the backgrounds


.. py:function:: monte(xdelta, files)

   This simply evaluates differences in a prediction
   for the background.  


.. py:function:: plot_fit_results(field, tile, summary_dir='', zlim=20)

.. py:function:: steer(argv)

   This is just a steering routine


.. py:function:: svdskyfit(xcross, threshold=0.1, new=True, verbose=True)

   Do a single value decomposition fit of the backgounds contained
   in the table xcross


.. py:function:: xeval(xall)

   xeval returns a metric that descibes now close we actually are to the
   original offsets.  Note that we do not have access to this information
   for real observations; it is only known here because we are using fake
   data.


.. py:function:: xpath(xdelta, files, nstart=2)

   Do a path search to find the offsets

   where xdelta is a table containing the differeces
   beween measurments of the "background" in overlapping
   regions, and files is a table listing the files,
   and the offsets for the files.  nstart is
   an index for the file that is used as the start 
   of the search.

   The routine returns the files table, with
   the offsets determined

   What is returned depends on the order of the
   files in the xdelta table



