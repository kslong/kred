ImageSum
========

.. py:module:: ImageSum

.. autoapi-nested-parse::

   ImageSum - Summarize FITS files in a directory

   Space Telescope Science Institute

   Synopsis
   --------

   Recursively search a directory for FITS files and produce a summary table
   with WCS-derived positions, image sizes, and header metadata. Useful for
   locating images that cover a given RA and Dec.

   Command Line Usage
   ------------------

   ::

       ImageSum.py [-h] [-out outname] dirname

   **Required Arguments:**

   dirname
       Directory to search recursively for FITS files.

   **Optional Arguments:**

   -h
       Print this documentation and exit.

   -out outname
       Set the output filename. Default: ``Image_Sum_{dirname}.txt``

   Description
   -----------

   Searches for all ``*.fits*`` files (including ``.fits.fz``) in the specified
   directory and subdirectories. For each file, extracts metadata from headers
   and calculates image center and size from WCS information.

   Handles multi-extension FITS (MEF) files by processing each image extension
   separately.

   Output Columns
   --------------

   Source_name   Object name from OBJECT header keyword
   Filter        Filter name from FILTER header keyword
   Exptime       Exposure time in seconds
   Image_type    Derived from filename (e.g., "ha", "sii", "r")
   RA            Right ascension of image center (degrees)
   Dec           Declination of image center (degrees)
   width         Image width in degrees (corrected for cos(dec))
   height        Image height in degrees
   mag           Photometric zeropoint (MAGZERO), or -999 if missing
   seeing        Seeing FWHM (SEEING), or -999 if missing
   filename      Full path to the FITS file
   ext           Extension number within the file

   Example
   -------

   ::

       ImageSum.py DECam_SUB2
       ImageSum.py -out all_tiles.txt DECam_SWARP

   Primary Routines
   ----------------

   table_create                        Main routine that builds the summary table
   list_image_extensions               Identify image extensions in MEF files
   get_image_center_and_size_from_header   Calculate WCS-derived image geometry

   History
   -------

   240813 ksl  Coding begun



Functions
---------

.. autoapisummary::

   ImageSum.get_image_center_and_size_from_header
   ImageSum.list_image_extensions
   ImageSum.print_header_image_info
   ImageSum.steer
   ImageSum.table_create


Module Contents
---------------

.. py:function:: get_image_center_and_size_from_header(header)

   Calculate the center coordinates (RA, Dec) and size in degrees 
   of a FITS image from its header containing WCS information.

   Parameters:
   -----------
   header : astropy.io.fits.Header
       FITS header containing WCS keywords and NAXIS information
       
   Returns:
   --------
   dict : Dictionary containing:
       - 'center_ra': RA of image center in degrees
       - 'center_dec': Dec of image center in degrees  
       - 'width_deg': Width of image in degrees
       - 'height_deg': Height of image in degrees
       - 'pixel_scale_ra_deg_per_pix': Pixel scale in RA direction (deg/pixel)
       - 'pixel_scale_dec_deg_per_pix': Pixel scale in Dec direction (deg/pixel)
       - 'image_shape': Tuple of (ny, nx) image dimensions
       - 'wcs': WCS object for further coordinate transformations


.. py:function:: list_image_extensions(filename)

.. py:function:: print_header_image_info(header, description='FITS Header')

   Convenience function to print image information from header in a readable format.

   Parameters:
   -----------
   header : astropy.io.fits.Header
       FITS header containing WCS and dimension information
   description : str
       Optional description for the output


.. py:function:: steer(argv)

.. py:function:: table_create(xdir='DECam_SUB2', outname='')

