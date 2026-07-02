PsfBuild
========

.. py:module:: PsfBuild

.. autoapi-nested-parse::

   Build PSF from star catalog

   Space Telescope Science Institute

   Synopsis
   --------

   Build PSF models from an image and star catalog. Selects optimal PSF stars
   from the input catalog and builds Gaussian, Moffat, and summed PSF models.

   Command Line Usage
   ------------------

   ::

       Usage: PsfBuild [-out root] image_file star_catalog

   where star_catalog is typically the all_stars.fits output from StarFind.

   Description
   -----------

   This module handles all PSF-related decisions:

   1. Selects optimal stars for PSF construction based on SNR, FWHM, eccentricity
   2. Writes the selected PSF stars to {prefix}_psf_stars.fits
   3. Builds Gaussian, Moffat, and summed PSF models
   4. Writes PSF models to FITS files

   Primary Routines
   ----------------

   select_psf_stars
       Select optimal stars for PSF construction from photometry table.

   do_one
       Build PSF from an image and star catalog.

   quick_psf_build
       Quick function to build all PSF models.

   Notes
   -----

   History:

   251210 ksl Coding begun
   251222 ksl Added to kred
   260113 ksl Moved select_psf_stars from StarFind; PsfBuild now handles all PSF decisions



Classes
-------

.. autoapisummary::

   PsfBuild.PSFBuilder


Functions
---------

.. autoapisummary::

   PsfBuild.do_one
   PsfBuild.quick_psf_build
   PsfBuild.select_psf_stars
   PsfBuild.steer


Module Contents
---------------

.. py:class:: PSFBuilder(fits_file, star_table, stamp_size=25, normalize='peak')

   Build PSF models from FITS images and star positions.

   Parameters
   ----------
   fits_file : str or HDUList
       Path to FITS file or astropy HDUList object
   star_table : astropy.Table
       Table with star positions (must contain 'xcenter' and 'ycenter' columns)
   stamp_size : int, optional
       Size of stamp to extract around each star (default: 25)
   normalize : str, optional
       Normalization method: 'peak' (default), 'sum', or 'none'


   .. py:method:: build_all_psfs()

      Build all three PSF models.

      Returns
      -------
      results : dict
          Dictionary containing all PSF models and parameters, or None if no stamps available



   .. py:method:: build_gaussian_psf()

      Build elliptical Gaussian PSF model.

      Returns
      -------
      params : dict
          Fitted parameters: amplitude, x0, y0, sigma_x, sigma_y, theta, fwhm_x, fwhm_y
      model : ndarray
          2D array of the fitted PSF model (normalized)



   .. py:method:: build_moffat_psf()

      Build elliptical Moffat PSF model.

      Returns
      -------
      params : dict
          Fitted parameters: amplitude, alpha, beta, ellipticity, theta, fwhm
      model : ndarray
          2D array of the fitted PSF model (normalized)



   .. py:method:: build_summed_psf(subtract_background=True)

      Build empirical PSF by summing star stamps.

      Parameters
      ----------
      subtract_background : bool, optional
          Subtract local background from each stamp (default: True)

      Returns
      -------
      psf : ndarray
          2D array of the summed PSF (normalized)



   .. py:method:: extract_stamps(filter_outliers=True, sigma_clip=3.0)

      Extract postage stamps around each star.

      Parameters
      ----------
      filter_outliers : bool, optional
          Remove stamps with unusual flux levels (default: True)
      sigma_clip : float, optional
          Sigma clipping threshold for outlier removal (default: 3.0)

      Returns
      -------
      stamps : list of dict
          List of stamps, each containing 'data', 'x', 'y', 'xcenter', 'ycenter'



   .. py:attribute:: normalize
      :value: 'peak'



   .. py:attribute:: psf_gaussian
      :value: None



   .. py:attribute:: psf_moffat
      :value: None



   .. py:attribute:: psf_summed
      :value: None



   .. py:method:: save_psfs(output_prefix='psf')

      Save PSF models to FITS files.

      Parameters
      ----------
      output_prefix : str, optional
          Prefix for output files (default: 'psf')



   .. py:attribute:: stamp_size
      :value: 25



   .. py:attribute:: stamps
      :value: None



   .. py:attribute:: star_table


.. py:function:: do_one(image_file, star_file, prefix='', max_stars=1000, weights=None)

   Build PSF from an image and star catalog.

   Parameters
   ----------
   image_file : str
       Path to FITS image file.
   star_file : str
       Path to star catalog (typically all_stars.fits from StarFind).
   prefix : str, optional
       Output filename prefix. If empty, derived from image_file.
   max_stars : int, optional
       Maximum number of stars to use for PSF. Default: 1000.
   weights : dict, optional
       Weights for quality score metrics. See select_psf_stars for details.


.. py:function:: quick_psf_build(fits_file, star_table, stamp_size=25, normalize='peak', output_prefix=None)

   Quick function to build all PSF models.

   Parameters
   ----------
   fits_file : str
       Path to FITS file
   star_table : astropy.Table
       Table with star positions (xcenter, ycenter columns)
   stamp_size : int, optional
       Size of PSF stamp (default: 25)
   normalize : str, optional
       Normalization method: 'peak' (default), 'sum', or 'none'
   output_prefix : str, optional
       If provided, save PSFs to files with this prefix

   Returns
   -------
   results : dict
       Dictionary containing all PSF models and parameters

   Example
   -------
   >>> from astropy.table import Table
   >>> star_table = Table.read('stars.fits')
   >>> results = quick_psf_build('image.fits', star_table, normalize='peak')
   >>> print(f"Gaussian FWHM: {results['gaussian']['params']['fwhm_x']:.2f}")


.. py:function:: select_psf_stars(phot_table, max_stars=100, weights=None, verbose=True)

   Select optimal stars for PSF construction using percentile-based quality scoring.

   Uses a weighted combination of quality metrics to rank stars, rather than
   hard cutoffs that may reject all stars in difficult fields.

   Parameters
   ----------
   phot_table : astropy.table.Table
       Output from do_forced_photometry with add_psf_metrics=True.
       Required columns: SNR, FWHM, Eccentricity, BkgContam, Concentration
   max_stars : int
       Maximum number of PSF stars to return (default: 100)
   weights : dict, optional
       Weights for each metric. Default weights emphasize SNR and FWHM consistency:
       {'snr': 0.35, 'fwhm': 0.25, 'ecc': 0.20, 'bkg': 0.10, 'conc': 0.10}
   verbose : bool
       Print diagnostic statistics (default: True)

   Returns
   -------
   psf_stars : astropy.table.Table
       Subset of highest quality stars, sorted by quality score


.. py:function:: steer(argv)

   This is generally just a steering routine

   Usage: PsfBuild [-out root] image_file star_catalog


