ImageBuilder
============

.. py:module:: ImageBuilder

.. autoapi-nested-parse::

   Construct an artificial image of a field given

   Space Telescope Science Institute

   Synopsis
   --------

   Construct an artificial image of a field given
   an original image, a psf, and a list of sources
   with fluxes (calculated using aperstats) from
   the image. A residual image is also constructed.

   PSF Image Generator Module

   Generate artificial images from PSF models and star catalogs.
   Completely independent from PSF building - just reads PSF FITS files and star tables.

   Command Line Usage
   ------------------

   ::

       ImageBuilder.py [-h] [-out root] [-radec] [-flux COL] [-mag] [-zp ZP]
                       original_image psf_file star_table

   Options
   -------

   -h
       Print this help message and exit

   -out root
       Output filename root. Default: derived from original_image

   -radec
       Compute pixel positions from RA/Dec columns in star_table using
       the WCS of the reference image. By default, expects xcenter/ycenter
       columns in the star table.

   -flux COL
       Column name for flux values in star_table. If not specified,
       tries 'Net' then 'flux_fit'.

   -mag
       Convert the flux column from magnitudes to flux values using
       the formula: flux = 10^((zeropoint - mag) / 2.5)

   -zp ZP
       Magnitude zeropoint for mag-to-flux conversion. This is the
       magnitude that corresponds to flux = 1. Default: 28.0

   Description
   -----------

   Primary routine is doit(), which calls generate_artificial_image() and
   generate_residual_image().

   Notes
   -----

   History:

   251216 ksl Coding begun
   251223 ksl Added to kred
   250116 ksl Added -radec, -flux, -mag, -zp options



Classes
-------

.. autoapisummary::

   ImageBuilder.PSFImageGenerator


Functions
---------

.. autoapisummary::

   ImageBuilder.doit
   ImageBuilder.generate_artificial_image
   ImageBuilder.generate_residual_image
   ImageBuilder.steer


Module Contents
---------------

.. py:class:: PSFImageGenerator(psf_file, reference_image=None, star_table=None)

   Generate artificial images using PSF models and star catalogs.

   Parameters
   ----------
   psf_file : str
       Path to PSF FITS file (from PSF builder or any compatible PSF)
   reference_image : str, optional
       Path to reference FITS image to copy WCS and dimensions from
   star_table : astropy.Table, optional
       Table with star positions and fluxes (xcenter, ycenter, Net columns)


   .. py:method:: add_star_to_image(image, x, y, flux, subpixel=True)

      Add a single star to the image at position (x, y) with given flux.

      Parameters
      ----------
      image : ndarray
          Image array to add star to (modified in place)
      x : float
          X position (0-indexed)
      y : float
          Y position (0-indexed)
      flux : float
          Total flux of the star
      subpixel : bool, optional
          Use subpixel positioning (default: True)



   .. py:method:: generate_image(image_shape=None, background=0.0, subpixel=True, star_table=None, flux_column='Net')

      Generate artificial image with all stars from catalog.

      Parameters
      ----------
      image_shape : tuple, optional
          Shape of output image (ny, nx). If None, uses reference image shape.
      background : float, optional
          Background level to add (default: 0.0)
      subpixel : bool, optional
          Use subpixel positioning (default: True)
      star_table : astropy.Table, optional
          Star table to use. If None, uses self.star_table
      flux_column : str, optional
          Column name for flux values (default: 'Net')

      Returns
      -------
      image : ndarray
          Generated image



   .. py:attribute:: image_shape
      :value: None



   .. py:method:: load_reference_image(reference_image)

      Load reference image to get WCS and dimensions.

      Parameters
      ----------
      reference_image : str
          Path to reference FITS image



   .. py:method:: load_star_table(star_table, use_radec=False, ra_col='RA', dec_col='Dec', flux_col=None)

      Load star table and create standardized columns (xcenter, ycenter, Net).

      After loading, the table will always have 'xcenter', 'ycenter', and 'Net'
      columns, regardless of the original column names. This allows subsequent
      transformations (e.g., magnitude to flux) to be applied to the 'Net' column.

      Parameters
      ----------
      star_table : str or astropy.Table
          Path to FITS table or Table object
      use_radec : bool, optional
          If True, compute xcenter and ycenter from RA/Dec columns using
          the WCS of the reference image. Default: False
      ra_col : str, optional
          Column name for RA if use_radec=True. Default: 'RA'
      dec_col : str, optional
          Column name for Dec if use_radec=True. Default: 'Dec'
      flux_col : str, optional
          Column name for flux values. The values will be copied to a 'Net'
          column. If None, tries 'Net' then 'flux_fit'. Default: None



   .. py:method:: mag_to_flux(zeropoint=28.0)

      Convert the Net column from magnitudes to flux values.

      Applies the transformation: flux = 10^((zeropoint - mag) / 2.5)

      With the default zeropoint of 28, a magnitude of 28 corresponds
      to a flux of 1.

      Parameters
      ----------
      zeropoint : float, optional
          The magnitude that corresponds to flux = 1. Default: 28.0

      Examples
      --------
      >>> generator.load_star_table('stars.fits', flux_col='rmag')
      >>> generator.mag_to_flux(zeropoint=25.0)  # mag 25 -> flux 1



   .. py:attribute:: psf_file


   .. py:attribute:: psf_half


   .. py:attribute:: psf_size


   .. py:attribute:: ref_header
      :value: None



   .. py:method:: save_image(image, output_file, copy_header=True)

      Save generated image to FITS file.

      Parameters
      ----------
      image : ndarray
          Image array to save
      output_file : str
          Output FITS filename
      copy_header : bool, optional
          Copy full header from reference image if available (default: True)



   .. py:method:: save_star_table(output_file=None)

      Save the processed star table to a FITS file.

      The table includes the standardized xcenter, ycenter, Net columns
      and any transformations applied (e.g., mag_to_flux).

      Parameters
      ----------
      output_file : str, optional
          Output filename. If None, derives from original filename by
          adding '.update' before '.fits' (e.g., 'stars.fits' -> 'stars.update.fits').
          If original was not a file, uses 'star_table.update.fits'.



   .. py:attribute:: star_table
      :value: None



   .. py:attribute:: wcs
      :value: None



.. py:function:: doit(psf_file='my_psf_gaussian.fits', star_table='forced.fits', reference_image='data/SMC_c01_T08.N673.fits', root='test', use_radec=False, flux_col=None, convert_mag=False, zeropoint=28.0)

.. py:function:: generate_artificial_image(psf_file, star_table, reference_image, output_file=None, background=0.0, subpixel=True, use_radec=False, ra_col='RA', dec_col='Dec', flux_col=None, convert_mag=False, zeropoint=28.0, save_table=True)

   Quick function to generate an artificial image.

   Parameters
   ----------
   psf_file : str
       Path to PSF FITS file
   star_table : str or astropy.Table
       Path to star catalog or Table object
   reference_image : str
       Path to reference FITS image (for WCS and dimensions)
   output_file : str, optional
       Output filename. If None, returns image without saving
   background : float, optional
       Background level (default: 0.0)
   subpixel : bool, optional
       Use subpixel positioning (default: True)
   use_radec : bool, optional
       Compute pixel positions from RA/Dec using WCS (default: False)
   ra_col : str, optional
       RA column name if use_radec=True (default: 'RA')
   dec_col : str, optional
       Dec column name if use_radec=True (default: 'Dec')
   flux_col : str, optional
       Column name for flux values. If None, tries 'Net' then 'flux_fit'.
   convert_mag : bool, optional
       Convert flux column from magnitudes to flux (default: False)
   zeropoint : float, optional
       Magnitude zeropoint for conversion (default: 28.0)
   save_table : bool, optional
       Save the processed star table (default: True)

   Returns
   -------
   image : ndarray
       Generated artificial image

   Example
   -------
   >>> image = generate_artificial_image('psf_gaussian.fits',
   ...                                    'stars.fits',
   ...                                    'science.fits',
   ...                                    output_file='artificial.fits')


.. py:function:: generate_residual_image(original_image, artificial_image, output_file=None)

   Generate residual image (original - artificial).

   Parameters
   ----------
   original_image : str or ndarray
       Original FITS image or array
   artificial_image : str or ndarray
       Artificial FITS image or array
   output_file : str, optional
       Output filename for residual

   Returns
   -------
   residual : ndarray
       Residual image


.. py:function:: steer(argv)

   usage ImageBuilder [-h] [-out root] [-radec] [-flux COL] [-mag] [-zp ZP]
                      original_image psf_file star_table


