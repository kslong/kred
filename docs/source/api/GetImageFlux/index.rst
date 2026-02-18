GetImageFlux
============

.. py:module:: GetImageFlux

.. autoapi-nested-parse::

   GetImageFlux - Extended Source Photometry from Region Files

   Space Telescope Science Institute

   Synopsis
   --------

   Calculate fluxes in source and background regions from FITS images using
   elliptical and annular apertures defined in a region table file. Computes
   both raw and background-subtracted fluxes for extended sources.

   Command Line Usage
   ------------------

   ::

       GetImageFlux.py [-h] [-viz] [-auto_back] [-gap N] image1.fits [image2.fits ...] region_table.txt
       GetImageFlux.py [-viz] [-auto_back] [-gap N] -match match_file.txt region_table.txt

   **Required Arguments:**

   image1.fits, image2.fits, ...
       One or more FITS files to process. Each file should have a primary
       (flux) extension and optionally a MASK extension for bad pixel exclusion.

   region_table.txt
       A region file containing source region definitions. If ``-auto_back`` is
       used, this file only needs source regions. Otherwise, it must contain both
       source and background region rows (SourceBack='Source' and 'Back').

   **Optional Arguments:**

   -h
       Display this help message and exit.

   -viz
       Create visualization plots for each source region. Plots are saved
       to the ``Figs_Flux/`` directory (created if necessary). Shows both
       source ellipse (red) and background annulus (magenta).

   -auto_back
       Automatically generate background annulus regions from a source-only
       input file. The background annulus is calculated to have the same area
       as the source ellipse. A new region file ``<input>_with_back.txt`` is
       created and used for processing.

   -gap N
       Gap in arcseconds between the source outer edge and the background
       annulus inner radius. Default is 3.0 arcsec. Only used with ``-auto_back``.

   -match match_file.txt
       Use a match file (output from ``ImageMatch2Source.py``) to drive
       which source/image combinations are processed. The match file must
       have columns ``Source_name`` and ``filename``. When ``-match`` is used,
       no explicit image files are needed on the command line — they come
       from the match file's ``filename`` column. Only the matched sources
       are processed for each image, rather than all sources in the region
       table.

   -filter name
       Insert a filter identifier into the output filename. For example,
       ``-filter ha`` produces ``Flux_ha_<region>.txt`` instead of
       ``Flux_<region>.txt``. Useful to avoid overwriting results from
       different filters.

   Description
   -----------

   This module performs aperture photometry on extended sources defined in
   a region table. For each source, it calculates:

   * Source flux from the defined elliptical region
   * Background flux from the corresponding annular region
   * Net flux (source minus scaled background)
   * An estimate of the flux error due to background uncertainty
   * Statistical measures (mean, median, mode, std)

   The net flux is computed as::

       net_flux = source_flux - num_pixels_used * back_median

   where ``source_flux`` is the summed pixel values in the source aperture
   (containing both the real source signal and background), ``back_median``
   is the median pixel value in the background annulus, and
   ``num_pixels_used`` is the number of unmasked pixels in the source
   aperture. The median is used because it is robust to contaminating
   sources within the background region.

   The flux error is estimated by splitting the background annulus into
   radial sub-annuli (default 4), computing the median in each, and taking
   the standard deviation of those medians. This measures how the
   background level varies with distance from the source, which is the
   dominant source of uncertainty for extended-source photometry. The
   per-pixel uncertainty is then scaled to the source aperture::

       flux_err = num_pixels_used * std(sub_annulus_medians)

   The region table must contain the following columns:

   * Source_name : str - Identifier for the source
   * RA : float - Right Ascension in degrees
   * Dec : float - Declination in degrees
   * RegType : str - Region type ('ellipse' or 'annulus')
   * Major : float - Semi-major axis in arcseconds
   * Minor : float - Semi-minor axis in arcseconds
   * Theta : float - Position angle in degrees
   * SourceBack : str - 'Source' or 'Back' to identify region type

   Output
   ------

   A single consolidated output file is written containing three row types
   for each source/image combination, identified by the ``SourceBack`` column:

   * **Source** -- photometry of the source region
   * **Back** -- photometry of the background region
   * **Net** -- background-subtracted values

   Key output columns:

   * ``flux`` -- total flux (DN) in the aperture; for Net rows this is background-subtracted
   * ``flux_err`` -- for Back rows, the per-pixel background uncertainty (std of
     radial sub-annulus medians); for Net rows, the total flux error
     (num_pixels_used * back_flux_err); 0 for Source rows
   * ``mean``, ``median`` -- per-pixel statistics; for Net rows, source minus background
   * ``std`` -- standard deviation of pixel values in the aperture
   * ``surface_brightness_per_arcsec2`` -- flux per square arcsecond
   * ``num_pixels_used`` -- number of unmasked pixels in the aperture
   * ``frac_in_image`` -- fraction of the aperture within the image (1.0 = fully contained)
   * ``Image`` -- which FITS file the measurement came from

   If ``-viz`` is specified, visualization plots are saved to:

   * ``Figs_Flux/<image>_<source_name>.png`` - Aperture visualization for each source

   If ``-auto_back`` is specified, a new region file is created:

   * ``<input>_with_back.txt`` - Region file with auto-generated background regions

   Notes
   -----

   **Source-only input format (for use with -auto_back):**

   When using ``-auto_back``, the input file only needs source definitions::

       Source_name        RA        Dec RegType  Major  Minor Theta
       -------------- --------- ---------- ------- ------ ------ -----
             SNR_N11L 73.702779 -66.429249 ellipse  35.00  35.00     0
              SNR_N86 73.904167 -68.646389 ellipse 183.00 183.00     0

   The background annulus is automatically calculated with:

   * Inner radius = max(Major, Minor) + gap (default gap = 3 arcsec)
   * Outer radius chosen so annulus area = source ellipse area

   **Full region table format (standard input):**

   Without ``-auto_back``, the input must have both source and background rows::

       Source_name        RA        Dec RegType  Major  Minor Theta SourceBack
       -------------- --------- ---------- ------- ------ ------ ----- ----------
             SNR_N11L 73.702779 -66.429249 ellipse  35.00  35.00     0     Source
             SNR_N11L 73.702779 -66.429249 annulus 105.00  70.00     0       Back
              SNR_N86 73.904167 -68.646389 ellipse 183.00 183.00     0     Source
              SNR_N86 73.904167 -68.646389 annulus 401.00 218.00     0       Back

   For 'ellipse' regions, Major and Minor are semi-axes of the ellipse.
   For 'annulus' regions, Major is the outer radius and Minor is the inner radius.

   History
   -------

   250808 ksl
       Initial coding

   250122 ksl
       Added Sphinx documentation

   .. moduleauthor:: KSL



Functions
---------

.. autoapisummary::

   GetImageFlux.add_net_rows
   GetImageFlux.background_uncertainty
   GetImageFlux.calculate_background_annulus
   GetImageFlux.circular_photometry
   GetImageFlux.create_bad_pixel_mask
   GetImageFlux.do_all
   GetImageFlux.do_many
   GetImageFlux.do_one
   GetImageFlux.elliptical_annulus_photometry
   GetImageFlux.elliptical_photometry
   GetImageFlux.elliptical_region_photometry
   GetImageFlux.generate_region_table
   GetImageFlux.get_pixel_scale
   GetImageFlux.load_fits_image
   GetImageFlux.print_photometry_summary
   GetImageFlux.results2table
   GetImageFlux.steer
   GetImageFlux.visualize_aperture_region


Module Contents
---------------

.. py:function:: add_net_rows(xtab)

   Add background-subtracted Net rows to a photometry table.

   For each source that has both Source and Back entries (matched by
   Source_name and Image), compute a Net row with background-subtracted
   flux and statistics. Rows are ordered as Source, Back, Net for each
   source in each image.

   The net flux is computed as::

       net_flux = source_flux - num_pixels_used_source * back_median

   where ``source_flux`` is the total (summed) flux in the source
   aperture (which includes both the real source signal and the
   background contribution), ``back_median`` is the median pixel value
   in the background annulus, and ``num_pixels_used_source`` is the
   number of unmasked pixels that contributed to the source flux.
   The median is used rather than the mean because it is robust to
   contaminating sources within the background annulus.

   The flux error is propagated from the per-pixel background
   uncertainty measured by ``background_uncertainty()``::

       net_flux_err = num_pixels_used_source * back_flux_err

   where ``back_flux_err`` is the standard deviation of the median
   pixel values measured in radial sub-annuli of the background region.

   Parameters
   ----------
   xtab : astropy.table.Table
       Table containing Source and Back entries identified by the
       SourceBack column.

   Returns
   -------
   astropy.table.Table
       Table with rows grouped as Source, Back, Net per source/image.
       The ``flux_err`` column contains: 0 for Source rows, the
       per-pixel background uncertainty for Back rows, and the
       propagated total flux error for Net rows.


.. py:function:: background_uncertainty(pixel_coord, r_inner_pix, r_outer_pix, data, bad_pixel_mask, n_sub=4)

   Estimate background uncertainty by radially subsampling an annulus.

   Splits the background annulus into concentric sub-annuli and computes
   the median in each. The standard deviation of these medians measures
   how the background level varies with radius, which is the dominant
   source of uncertainty for extended source photometry.

   Parameters
   ----------
   pixel_coord : tuple
       (x, y) pixel coordinates of the annulus center.
   r_inner_pix : float
       Inner radius of the background annulus in pixels.
   r_outer_pix : float
       Outer radius of the background annulus in pixels.
   data : numpy.ndarray
       2D image data array.
   bad_pixel_mask : numpy.ndarray
       Boolean mask where True indicates bad pixels.
   n_sub : int, optional
       Number of radial sub-annuli. Default is 4.

   Returns
   -------
   float
       Standard deviation of the sub-annulus medians (per-pixel
       background uncertainty). Returns 0.0 if fewer than 2
       sub-annuli have valid data.


.. py:function:: calculate_background_annulus(a_arcsec, b_arcsec, gap=3.0)

   Calculate background annulus radii to match source ellipse area.

   Given a source ellipse, calculates the inner and outer radii of a
   circular background annulus such that:
   - Inner radius = max(a, b) + gap (starts outside source with buffer)
   - Outer radius chosen so annulus area equals source ellipse area

   Parameters
   ----------
   a_arcsec : float
       Semi-major axis of source ellipse in arcseconds.
   b_arcsec : float
       Semi-minor axis of source ellipse in arcseconds.
   gap : float, optional
       Gap between source outer edge and background inner radius in
       arcseconds. Default is 3.0.

   Returns
   -------
   r_inner : float
       Inner radius of background annulus in arcseconds.
   r_outer : float
       Outer radius of background annulus in arcseconds.

   Notes
   -----
   The area calculation:
   - Source ellipse area: A = π * a * b
   - Annulus area: A = π * (r_out² - r_in²)
   - Setting equal: r_out = sqrt(r_in² + a * b)

   Examples
   --------
   >>> r_in, r_out = calculate_background_annulus(35.0, 35.0, gap=3.0)
   >>> print(f"Inner: {r_in:.1f}, Outer: {r_out:.1f}")
   Inner: 38.0, Outer: 51.7


.. py:function:: circular_photometry(fits_file, ra, dec, radius_arcsec, radius_in_arcsec=None, include_zero_mask=True, image_data=None)

   Perform circular aperture photometry.

   Convenience wrapper using elliptical functions with equal semi-axes.

   Parameters
   ----------
   fits_file : str
       Path to the FITS file.
   ra : float
       Right Ascension in degrees.
   dec : float
       Declination in degrees.
   radius_arcsec : float
       Outer radius in arcseconds.
   radius_in_arcsec : float, optional
       Inner radius for annulus in arcseconds. If provided, creates
       a circular annulus.
   include_zero_mask : bool, optional
       Whether to mask zero-valued pixels. Default is True.
   image_data : dict, optional
       Pre-loaded image data from load_fits_image(). Default is None.

   Returns
   -------
   dict
       Photometry results dictionary. See elliptical_photometry().


.. py:function:: create_bad_pixel_mask(data, value_range=None, specific_values=None)

   Create a boolean mask for bad pixels in an image.

   Identifies bad pixels based on NaN values, value range limits,
   and specific flagged values.

   Parameters
   ----------
   data : numpy.ndarray
       2D (or N-D) array of image data.
   value_range : tuple of (float, float), optional
       Mask pixels with values outside this (min, max) range.
       If None, no range masking is applied.
   specific_values : list, set, or tuple, optional
       Mask pixels matching these exact values (e.g., {0, -999, 9999}).

   Returns
   -------
   numpy.ndarray
       Boolean mask where True indicates bad pixels.

   Examples
   --------
   >>> mask = create_bad_pixel_mask(image_data, value_range=(-100, 100))
   >>> masked_data = np.ma.array(image_data, mask=mask)


.. py:function:: do_all(image_file, region_table, create_visualization=False, source_names=None)

   Process all sources in a region table for one image.

   Parameters
   ----------
   image_file : str
       Path to the FITS image file.
   region_table : str
       Path to the region table file.
   create_visualization : bool, optional
       If True, create visualization plots for each source region.
       Plots are saved to ``Figs_Flux/`` directory. Default is False.
   source_names : list of str, optional
       If provided, only process sources whose Source_name is in this list.
       Used by ``-match`` mode to restrict processing to matched sources.
       Default is None (process all sources).

   Returns
   -------
   astropy.table.Table or None
       Photometry results table with Source and Back rows and an
       'Image' column, or None if no regions overlap the image.


.. py:function:: do_many(xtab, image_file, create_visualization=False)

   Process all regions in a table for a given image.

   Iterates through the region table and performs photometry for each
   source and background region. When visualization is enabled, creates
   plots showing both the source region and background annulus.

   Parameters
   ----------
   xtab : astropy.table.Table
       Region table with columns: Source_name, RA, Dec, RegType,
       Major, Minor, Theta, SourceBack.
   image_file : str
       Path to the FITS image file.
   create_visualization : bool, optional
       If True, create visualization plots for source regions showing
       both the source ellipse and background annulus. Default is False.

   Returns
   -------
   list of dict
       List of photometry result dictionaries, one per region.

   Raises
   ------
   ValueError
       If an unknown RegType is encountered.


.. py:function:: do_one(source_name='SNR_N49', region_file='lmc_spec_ann_reg.txt', image_file='../LMC/LMC.ha.csub.fits.gz')

   Process a single source from a region file.

   Extracts photometry for one source and its background region,
   useful for testing and debugging.

   Parameters
   ----------
   source_name : str, optional
       Name of the source to process. Default is 'SNR_N49'.
   region_file : str, optional
       Path to the region table file.
   image_file : str, optional
       Path to the FITS image file.

   Returns
   -------
   list of dict
       Two-element list containing [source_results, background_results].


.. py:function:: elliptical_annulus_photometry(fits_file, ra, dec, a_out_arcsec, b_out_arcsec, a_in_arcsec, b_in_arcsec, theta_deg=0, include_zero_mask=True, image_data=None)

   Perform photometry in an elliptical annulus.

   Convenience wrapper for elliptical_photometry() for annular apertures.

   Parameters
   ----------
   fits_file : str
       Path to the FITS file.
   ra : float
       Right Ascension in degrees.
   dec : float
       Declination in degrees.
   a_out_arcsec : float
       Semi-major axis of outer ellipse in arcseconds.
   b_out_arcsec : float
       Semi-minor axis of outer ellipse in arcseconds.
   a_in_arcsec : float
       Semi-major axis of inner ellipse in arcseconds.
   b_in_arcsec : float
       Semi-minor axis of inner ellipse in arcseconds.
   theta_deg : float, optional
       Position angle in degrees. Default is 0.
   include_zero_mask : bool, optional
       Whether to mask zero-valued pixels. Default is True.
   image_data : dict, optional
       Pre-loaded image data from load_fits_image(). Default is None.

   Returns
   -------
   dict
       Photometry results dictionary. See elliptical_photometry().


.. py:function:: elliptical_photometry(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0, a_in_arcsec=None, b_in_arcsec=None, include_zero_mask=True, image_data=None)

   Perform elliptical aperture photometry at a sky position.

   Unified function that handles both simple elliptical apertures and
   elliptical annuli for background estimation.

   Parameters
   ----------
   fits_file : str
       Path to the FITS file containing the image data.
   ra : float
       Right Ascension of the aperture center in degrees.
   dec : float
       Declination of the aperture center in degrees.
   a_arcsec : float
       Semi-major axis of outer ellipse in arcseconds.
   b_arcsec : float
       Semi-minor axis of outer ellipse in arcseconds.
   theta_deg : float, optional
       Position angle of ellipse in degrees, measured counter-clockwise
       from the positive x-axis. Default is 0.
   a_in_arcsec : float, optional
       Semi-major axis of inner ellipse in arcseconds. If provided along
       with b_in_arcsec, creates an elliptical annulus.
   b_in_arcsec : float, optional
       Semi-minor axis of inner ellipse in arcseconds.
   include_zero_mask : bool, optional
       Whether to mask pixels with zero values. Default is True.
   image_data : dict, optional
       Pre-loaded image data from load_fits_image(). If provided, the FITS
       file is not opened again. Default is None.

   Returns
   -------
   dict
       Dictionary containing photometry results with keys:

       * flux : Total flux in aperture
       * surface_brightness_per_pixel : Flux per pixel
       * surface_brightness_per_arcsec2 : Flux per square arcsecond
       * area_pixels, area_arcsec2 : Aperture areas
       * num_pixels_total, num_pixels_used, num_pixels_masked : Pixel counts
       * mean, median, mode, std, min, max : Statistical measures
       * aperture_type : 'elliptical_aperture' or 'elliptical_annulus'
       * Aperture geometry parameters

   Notes
   -----
   If the FITS file contains a 'MASK' extension, it will be used for bad
   pixel masking. Otherwise, a mask is created from NaN values and the
   specified value_range.

   Examples
   --------
   >>> results = elliptical_photometry('image.fits', 73.7, -66.4, 35.0, 35.0)
   >>> print(f"Flux: {results['flux']:.2e}")


.. py:function:: elliptical_region_photometry(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0, include_zero_mask=True, image_data=None)

   Perform photometry in a simple elliptical aperture.

   Convenience wrapper for elliptical_photometry() for simple apertures
   (no annulus).

   Parameters
   ----------
   fits_file : str
       Path to the FITS file.
   ra : float
       Right Ascension in degrees.
   dec : float
       Declination in degrees.
   a_arcsec : float
       Semi-major axis in arcseconds.
   b_arcsec : float
       Semi-minor axis in arcseconds.
   theta_deg : float, optional
       Position angle in degrees. Default is 0.
   include_zero_mask : bool, optional
       Whether to mask zero-valued pixels. Default is True.
   image_data : dict, optional
       Pre-loaded image data from load_fits_image(). Default is None.

   Returns
   -------
   dict
       Photometry results dictionary. See elliptical_photometry().


.. py:function:: generate_region_table(input_file, gap=3.0)

   Generate a complete region table with auto-generated background regions.

   Reads a source-only region file and generates corresponding background
   annulus regions for each source, with area matching the source region.

   Parameters
   ----------
   input_file : str
       Path to input region file containing source regions only.
       Must have columns: Source_name, RA, Dec, RegType, Major, Minor, Theta.
   gap : float, optional
       Gap between source outer edge and background inner radius in
       arcseconds. Default is 3.0.

   Returns
   -------
   output_file : str
       Path to the generated output file with both source and background
       regions. Named as ``<input_basename>_with_back.txt``.
   combined_table : astropy.table.Table
       Table containing both source and background region rows.

   Notes
   -----
   The output table includes all columns from the input plus:
   - SourceBack : 'Source' or 'Back' to identify region type

   For each source row, a corresponding background row is created with:
   - RegType = 'annulus'
   - Major = outer radius of background annulus
   - Minor = inner radius of background annulus
   - SourceBack = 'Back'

   The output file is written in ascii.fixed_width_two_line format.

   Examples
   --------
   >>> outfile, table = generate_region_table('sources.txt', gap=5.0)
   >>> print(f"Generated {outfile} with {len(table)} rows")


.. py:function:: get_pixel_scale(fits_header)

   Calculate the pixel scale from a FITS header.

   Extracts the pixel scale using either CDELT keywords or the CD matrix
   from the WCS information in the header.

   Parameters
   ----------
   fits_header : astropy.io.fits.Header
       FITS file header containing WCS information.

   Returns
   -------
   float
       Pixel scale in arcseconds per pixel.

   Notes
   -----
   The function first attempts to use CDELT1/CDELT2 keywords if present.
   If not available, it calculates the scale from the CD matrix via the WCS.
   The returned value is the mean of the two axes' scales.


.. py:function:: load_fits_image(fits_file)

   Load a FITS image and return all derived quantities needed for photometry.

   Reads the FITS file once and extracts the data, header, bad pixel mask,
   WCS, and pixel scale. This avoids redundant I/O when processing multiple
   regions on the same image.

   Parameters
   ----------
   fits_file : str
       Path to the FITS file containing the image data.

   Returns
   -------
   dict
       Dictionary with keys:

       * data : numpy.ndarray - 2D image data from the primary extension
       * header : astropy.io.fits.Header - Primary extension header
       * bad_pixel_mask : numpy.ndarray - Boolean mask (True = bad pixel)
       * wcs : astropy.wcs.WCS - World Coordinate System object
       * pixel_scale : float - Pixel scale in arcseconds per pixel


.. py:function:: print_photometry_summary(results)

   Print a formatted summary of photometry results.

   Parameters
   ----------
   results : dict
       Photometry results dictionary from elliptical_photometry().

   Returns
   -------
   None
       Prints summary to stdout.


.. py:function:: results2table(results_list)

   Convert a list of result dictionaries to an Astropy Table.

   Parameters
   ----------
   results_list : list of dict
       List of photometry result dictionaries.

   Returns
   -------
   astropy.table.Table
       Table with formatted columns. Float columns use '.2f' format,
       and None values are converted to -99.0. The ``flux_err`` column
       is placed immediately after ``flux``.


.. py:function:: steer(argv)

   Parse command line arguments and execute photometry.

   Parameters
   ----------
   argv : list of str
       Command line arguments (typically sys.argv).

   Returns
   -------
   None
       Results are written to output files.

   Notes
   -----
   Arguments are parsed as follows:

   * Arguments containing 'fits' are treated as image files
   * The first non-fits argument is treated as the region table
   * ``-h`` prints help and exits
   * ``-viz`` enables visualization output to Figs_Flux/
   * ``-auto_back`` auto-generates background regions from source-only input
   * ``-gap N`` sets the gap between source and background (default 3 arcsec)
   * ``-match match_file`` uses a match file to drive source/image pairs
   * ``-filter name`` inserts filter name into output filename


.. py:function:: visualize_aperture_region(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0, back_outer_arcsec=None, back_inner_arcsec=None, display_size_arcsec=None, output_filename=None, include_zero_mask=True, show_plot=True, source_name=None, image_data=None)

   Create a visualization showing source and background regions on image data.

   Generates a three-panel plot showing: (1) original data with source ellipse
   (red) and background annulus (magenta), (2) masked data, and (3) source region
   mask overlay.

   Parameters
   ----------
   fits_file : str
       Path to the FITS file.
   ra : float
       Right Ascension of aperture center in degrees.
   dec : float
       Declination of aperture center in degrees.
   a_arcsec : float
       Semi-major axis of source ellipse in arcseconds.
   b_arcsec : float
       Semi-minor axis of source ellipse in arcseconds.
   theta_deg : float, optional
       Position angle in degrees. Default is 0.
   back_outer_arcsec : float, optional
       Outer radius of background annulus in arcseconds.
   back_inner_arcsec : float, optional
       Inner radius of background annulus in arcseconds.
   display_size_arcsec : float, optional
       Size of cutout region in arcseconds. Default is 2.5x the background
       outer radius (if provided) or 4x the source major axis.
   output_filename : str, optional
       Output filename for plot. Auto-generated if None.
   include_zero_mask : bool, optional
       Whether to mask zero-valued pixels in display. Default is True.
   show_plot : bool, optional
       Whether to display plot interactively. Default is True.
   source_name : str, optional
       Name of the source for the plot title.
   image_data : dict, optional
       Pre-loaded image data from load_fits_image(). If provided, the FITS
       file is not opened again. Default is None.

   Returns
   -------
   str
       Path to the saved visualization image file.


