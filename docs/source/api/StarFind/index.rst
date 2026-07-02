StarFind
========

.. py:module:: StarFind

.. autoapi-nested-parse::

   StarFind - Star Identification for PSF Creation

   Space Telescope Science Institute

   Synopsis
   --------

   Identifies stars suitable for creating PSF functions from astronomical images.

   Command Line Usage
   ------------------

   ::

       Usage: StarFind.py -h -out root file1 file2 etc.

   Description
   -----------

   This routine identifies stars suitable for PSF (Point Spread Function) creation.

   Version History
   ---------------
   251210 ksl
       Coding begun

   251222 ksl
       Added to kred



Attributes
----------

.. autoapisummary::

   StarFind.XDIR


Functions
---------

.. autoapisummary::

   StarFind.do_forced_photometry
   StarFind.do_one
   StarFind.get_objects_from_image
   StarFind.read_table
   StarFind.steer


Module Contents
---------------

.. py:data:: XDIR
   :value: ''


.. py:function:: do_forced_photometry(filename='LMC_c48_T08.r.t060.fits', image_ext=0, object_file='objects.txt', nrows_max=-1, rstar=6, b_in=8, b_out=12, add_psf_metrics=False)

   Do forced photometry based on RA/Dec positions from object file.

   Parameters
   ----------
   filename : str, optional
       FITS file with one or more image extensions. Default is 'LMC_c48_T08.r.t060.fits'.
   image_ext : int, optional
       Extension to be analyzed. Default is 0.
   object_file : str, optional
       File containing source positions. Default is 'objects.txt'.
   nrows_max : int, optional
       Maximum number of objects for forced photometry. If -1, processes all. Default is -1.
   rstar : float, optional
       Aperture radius in pixels for source extraction. Default is 6.
   b_in : float, optional
       Inner radius of background annulus in pixels. Default is 8.
   b_out : float, optional
       Outer radius of background annulus in pixels. Default is 12.
   add_psf_metrics : bool, optional
       If True, adds columns useful for PSF star selection. Default is False.

   Returns
   -------
   phot_table : astropy.table.Table
       Photometry results table.

   Notes
   -----
   This function returns the table but does not write it to disk.



.. py:function:: do_one(filename, root='test')

   Find stars and perform photometry on an image.

   Parameters
   ----------
   filename : str
       Path to FITS image file.
   root : str, optional
       Output filename root. If empty, derived from filename.

   Notes
   -----
   Writes {root}_all_stars.fits containing all detected stars with photometry.
   Use PsfBuild to select PSF stars and build PSF models.


.. py:function:: get_objects_from_image(filename='LMC_c48_T08.r.t060.fits', outroot='')

.. py:function:: read_table(filename)

   This is a generic routine to try to read a table
   in fits or ascii format.  It is intended to accommodate 
   several different types of formats.


.. py:function:: steer(argv)

   This is generally just a steering routine

   Usage: StarFind -h -out root file1 file2 etc.


