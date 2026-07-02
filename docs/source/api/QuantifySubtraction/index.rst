QuantifySubtraction
===================

.. py:module:: QuantifySubtraction

.. autoapi-nested-parse::

   Quantify in various ways how good star subtraction

   Space Telescope Science Institute

   Synopsis
   --------

   Quantify in various ways how good star subtraction
   is between two images, in which star fluxes have been
   calculated using MefPhot

   Command Line Usage
   ------------------

   ::

       usage: QuantifySubtraction.py orig_file subtracted_file [more pairs...]

   Description
   -----------

   Primary routines:

       doit - Process a single pair of files
       do_many - Process multiple pairs and create comparison table

   Primary Routines
   ----------------

   doit - Process a single pair of files
       do_many - Process multiple pairs and create comparison table

   Notes
   -----

   History:

   251229 ksl Coding begun
   251229 ksl Refactored for modularity and result capture

   Version History
   ---------------

   251229 ksl
       Coding begun

   251229 ksl
       Refactored for modularity and result capture



Functions
---------

.. autoapisummary::

   QuantifySubtraction.calculate_metrics
   QuantifySubtraction.create_combined_figure
   QuantifySubtraction.do_many
   QuantifySubtraction.doit
   QuantifySubtraction.get_pairs
   QuantifySubtraction.steer


Module Contents
---------------

.. py:function:: calculate_metrics(orig, subtracted)

   Calculate all subtraction quality metrics

   Returns a dictionary with summary statistics


.. py:function:: create_combined_figure(orig, subtracted, outroot)

   Create a single combined figure with all analysis types


.. py:function:: do_many(file_pairs, output_table='SubQual/comparison_metrics.fits')

   Process multiple file pairs and create comparison table

   Parameters:
   -----------
   file_pairs : list of tuples
       List of (orig_phot, subtracted_phot) pairs
   output_table : str
       Path to save the comparison metrics table

   Returns:
   --------
   results : astropy.table.Table
       Table with metrics for all file pairs


.. py:function:: doit(orig_phot, subtracted_phot, outroot='', combined=True)

   Analyze subtraction quality for a single pair of files

   Parameters:
   -----------
   orig_phot : str
       Path to original photometry table
   subtracted_phot : str
       Path to subtracted photometry table
   outroot : str
       Root name for output files (derived from subtracted_phot if empty)
   combined : bool
       If True, create a single combined figure

   Returns:
   --------
   metrics : dict
       Dictionary of calculated metrics


.. py:function:: get_pairs(xdir='TabPhot', prefix='LMC')

.. py:function:: steer(argv)

   Steering routine for command line usage


