#!/usr/bin/env python
# coding: utf-8
"""GetImageFlux - Extended Source Photometry from Region Files

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

Description
-----------

This module performs aperture photometry on extended sources defined in
a region table. For each source, it calculates:

* Source flux from the defined elliptical region
* Background flux from the corresponding annular region
* Net flux (source minus scaled background)
* Statistical measures (mean, median, mode, std)

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

For each image, two output files are created:

* ``SB_<image>.<region>.txt`` - Full photometry results for all regions
* ``Net_<image>.<region>.txt`` - Net flux summary for each source

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
"""

import sys



import warnings
from astropy.wcs import WCS, FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning, message=".*datfix.*")


from astropy.io import ascii
from astropy.table import Table,join, vstack
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.aperture import ApertureStats




def get_pixel_scale(fits_header):
    """Calculate the pixel scale from a FITS header.

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
    """
    # Always use WCS pixel_scale_matrix, which correctly handles
    # both the CD matrix and CDELT+PC matrix conventions.
    wcs = WCS(fits_header)
    cd_matrix = wcs.pixel_scale_matrix
    cdelt1 = np.sqrt(cd_matrix[0, 0]**2 + cd_matrix[0, 1]**2)
    cdelt2 = np.sqrt(cd_matrix[1, 0]**2 + cd_matrix[1, 1]**2)

    pixel_scale = np.mean([cdelt1, cdelt2]) * 3600
    return pixel_scale


def create_bad_pixel_mask(data, value_range=None, specific_values=None):
    """Create a boolean mask for bad pixels in an image.

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
    """
    # Start with NaN mask
    mask = np.isnan(data)

    # Mask outside range
    if value_range is not None:
        min_val, max_val = value_range
        mask |= (data < min_val) | (data > max_val)

    # Mask specific values
    if specific_values is not None:
        for val in specific_values:
            mask |= (data == val)

    return mask


def calculate_background_annulus(a_arcsec, b_arcsec, gap=3.0):
    """Calculate background annulus radii to match source ellipse area.

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
    """
    # Inner radius: outside the source ellipse plus gap
    r_inner = max(a_arcsec, b_arcsec) + gap

    # Outer radius: annulus area = ellipse area
    # π * (r_out² - r_in²) = π * a * b
    # r_out = sqrt(r_in² + a * b)
    source_area_term = a_arcsec * b_arcsec
    r_outer = np.sqrt(r_inner**2 + source_area_term)

    return r_inner, r_outer


def generate_region_table(input_file, gap=3.0):
    """Generate a complete region table with auto-generated background regions.

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
    """
    import os

    # Read input table
    try:
        source_tab = ascii.read(input_file)
    except Exception as e:
        print(f'Error: could not read {input_file}: {e}')
        return None, None

    # Check for required columns
    required_cols = ['Source_name', 'RA', 'Dec', 'RegType', 'Major', 'Minor', 'Theta']
    missing = [col for col in required_cols if col not in source_tab.colnames]
    if missing:
        print(f'Error: input file missing required columns: {missing}')
        return None, None

    # Add SourceBack column if not present
    if 'SourceBack' not in source_tab.colnames:
        source_tab['SourceBack'] = 'Source'

    # Create list to hold all rows (source + background)
    all_rows = []

    for row in source_tab:
        # Add the source row
        source_dict = {col: row[col] for col in source_tab.colnames}
        source_dict['SourceBack'] = 'Source'
        all_rows.append(source_dict)

        # Generate background annulus
        a = row['Major']
        b = row['Minor']
        r_inner, r_outer = calculate_background_annulus(a, b, gap)

        # Create background row (copy source row and modify)
        back_dict = {col: row[col] for col in source_tab.colnames}
        back_dict['RegType'] = 'annulus'
        back_dict['Major'] = r_outer  # Outer radius
        back_dict['Minor'] = r_inner  # Inner radius
        back_dict['SourceBack'] = 'Back'
        all_rows.append(back_dict)

    # Create combined table
    combined_table = Table(rows=all_rows)

    # Format numeric columns
    for col in ['RA', 'Dec', 'Major', 'Minor', 'Theta']:
        if col in combined_table.colnames:
            if combined_table[col].dtype in [np.float64, np.float32]:
                combined_table[col].format = '.2f'

    # Generate output filename
    base = os.path.basename(input_file)
    # Remove extension
    if base.endswith('.txt'):
        base = base[:-4]
    elif base.endswith('.tab'):
        base = base[:-4]
    output_file = f'{base}_with_back.txt'

    # Write output file
    combined_table.write(output_file, format='ascii.fixed_width_two_line', overwrite=True)
    print(f'Generated region table: {output_file} with {len(combined_table)} rows')
    print(f'  ({len(source_tab)} sources + {len(source_tab)} background regions)')
    print(f'  Background gap: {gap:.1f} arcsec, area-matched annuli')

    return output_file, combined_table


def elliptical_photometry(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0,
                         a_in_arcsec=None, b_in_arcsec=None, include_zero_mask=True):
    """Perform elliptical aperture photometry at a sky position.

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
    """
    from photutils.aperture import EllipticalAperture, EllipticalAnnulus
    
    # Open the FITS file
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        try:
            bad_pixel_mask=hdul['MASK'].data
            bad_pixel_mask=bad_pixel_mask.astype('bool')
        except:
            bad_pixel_mask= create_bad_pixel_mask(data, value_range=[-100,100], specific_values=None)
        
    
    # Get the WCS (World Coordinate System) from the FITS header

    wcs = WCS(header)

    
    # Convert the RA and Dec to pixel coordinates
    sky_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    pixel_coord = wcs.world_to_pixel(sky_coord)

    # Get the pixel scale
    pixel_scale = get_pixel_scale(header)

    # Convert semi-axes from arcseconds to pixels
    a_pixels = a_arcsec / pixel_scale
    b_pixels = b_arcsec / pixel_scale
    theta_rad = np.deg2rad(theta_deg)

    # Check if the aperture could possibly overlap the image
    px, py = float(pixel_coord[0]), float(pixel_coord[1])
    max_radius = max(a_pixels, b_pixels)
    ny, nx = data.shape
    if (px + max_radius < 0 or px - max_radius >= nx or
        py + max_radius < 0 or py - max_radius >= ny):
        return None

    # Determine if this is an annulus or simple aperture
    is_annulus = (a_in_arcsec is not None) and (b_in_arcsec is not None)
    
    if is_annulus:
        # Convert inner semi-axes to pixels
        a_in_pixels = a_in_arcsec / pixel_scale
        b_in_pixels = b_in_arcsec / pixel_scale
        
        # Create elliptical annulus
        aperture = EllipticalAnnulus(pixel_coord, 
                                   a_in=a_in_pixels, a_out=a_pixels,
                                   b_in=b_in_pixels, b_out=b_pixels,
                                   theta=theta_rad)
        
        # Calculate areas
        area_outer = np.pi * a_pixels * b_pixels
        area_inner = np.pi * a_in_pixels * b_in_pixels
        area_total = area_outer - area_inner
        
        aperture_type = 'elliptical_annulus'
    else:
        # Create simple elliptical aperture
        aperture = EllipticalAperture(pixel_coord, a=a_pixels, b=b_pixels, theta=theta_rad)
        
        # Calculate area
        area_total = np.pi * a_pixels * b_pixels
        area_outer = area_total
        area_inner = 0
        
        aperture_type = 'elliptical_aperture'

    # Perform aperture photometry
    phot_table = aperture_photometry(data, aperture, mask=bad_pixel_mask)
    
    # Calculate aperture statistics
    aperture_stats = ApertureStats(data, aperture, mask=bad_pixel_mask)
    
    # Extract flux
    flux = phot_table['aperture_sum'][0]

    # Calculate pixel counts and fraction of aperture in image
    aperture_mask = aperture.to_mask()
    full_aperture_pixels = 0
    in_image_pixels = 0

    if hasattr(aperture_mask, '__len__'):  # Multiple masks for annulus
        total_mask_data = np.zeros_like(data, dtype=bool)
        for mask in aperture_mask:
            if mask is not None:
                full_aperture_pixels += np.sum(mask.data > 0)
                mask_array = mask.to_image(data.shape)
                if mask_array is not None:
                    total_mask_data |= (mask_array > 0)
        in_image_pixels = np.sum(total_mask_data)
        num_pixels_total = in_image_pixels
        masked_pixels_in_aperture = np.sum(bad_pixel_mask & total_mask_data)
    else:
        if aperture_mask is not None:
            full_aperture_pixels = np.sum(aperture_mask.data > 0)
            aperture_mask_array = aperture_mask.to_image(data.shape)
            if aperture_mask_array is not None:
                in_image_pixels = np.sum(aperture_mask_array > 0)
                masked_pixels_in_aperture = np.sum(bad_pixel_mask & (aperture_mask_array > 0))
            else:
                masked_pixels_in_aperture = 0
            num_pixels_total = in_image_pixels
        else:
            num_pixels_total = 0
            masked_pixels_in_aperture = 0

    frac_in_image = in_image_pixels / full_aperture_pixels if full_aperture_pixels > 0 else 0
    num_pixels_used = num_pixels_total - masked_pixels_in_aperture
    
    # Calculate surface brightness (flux per unit area)
    if area_total > 0:
        surface_brightness_per_pixel = flux / area_total
        surface_brightness_per_arcsec2 = flux / (area_total * pixel_scale**2)
    else:
        surface_brightness_per_pixel = 0
        surface_brightness_per_arcsec2 = 0
    
    # Compile comprehensive results
    results = {
        # Basic photometry
        'flux': flux,
        'surface_brightness_per_pixel': surface_brightness_per_pixel,
        'surface_brightness_per_arcsec2': surface_brightness_per_arcsec2,
        
        # Areas
        'area_pixels': area_total,
        'area_arcsec2': area_total * pixel_scale**2,
        'area_outer_pixels': area_outer,
        'area_inner_pixels': area_inner,
        
        # Pixel statistics
        'num_pixels_total': num_pixels_total,
        'num_pixels_used': num_pixels_used,
        'num_pixels_masked': masked_pixels_in_aperture,
        'fraction_pixels_used': num_pixels_used / num_pixels_total if num_pixels_total > 0 else 0,
        'frac_in_image': frac_in_image,
        
        # Statistical measures
        'mean': aperture_stats.mean,
        'median': aperture_stats.median,
        'mode': getattr(aperture_stats, 'mode', None),  # mode might not always be available
        'std': aperture_stats.std,
        'min': aperture_stats.min,
        'max': aperture_stats.max,
        
        # Aperture parameters
        'aperture_type': aperture_type,
        'a_pixels': a_pixels,
        'b_pixels': b_pixels,
        'a_arcsec': a_arcsec,
        'b_arcsec': b_arcsec,
        'theta_deg': theta_deg,
        'pixel_scale': pixel_scale,
        
        # Additional parameters for annulus
        'a_in_pixels': a_in_pixels if is_annulus else None,
        'b_in_pixels': b_in_pixels if is_annulus else None,
        'a_in_arcsec': a_in_arcsec if is_annulus else None,
        'b_in_arcsec': b_in_arcsec if is_annulus else None,
        'is_annulus': is_annulus
    }

    return results


def elliptical_region_photometry(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0,
                                include_zero_mask=True):
    """Perform photometry in a simple elliptical aperture.

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

    Returns
    -------
    dict
        Photometry results dictionary. See elliptical_photometry().
    """
    return elliptical_photometry(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg,
                               include_zero_mask=include_zero_mask)


def elliptical_annulus_photometry(fits_file, ra, dec, a_out_arcsec, b_out_arcsec,
                                 a_in_arcsec, b_in_arcsec, theta_deg=0, include_zero_mask=True):
    """Perform photometry in an elliptical annulus.

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

    Returns
    -------
    dict
        Photometry results dictionary. See elliptical_photometry().
    """
    return elliptical_photometry(fits_file, ra, dec, a_out_arcsec, b_out_arcsec, theta_deg,
                               a_in_arcsec, b_in_arcsec, include_zero_mask)


def circular_photometry(fits_file, ra, dec, radius_arcsec, radius_in_arcsec=None,
                       include_zero_mask=True):
    """Perform circular aperture photometry.

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

    Returns
    -------
    dict
        Photometry results dictionary. See elliptical_photometry().
    """
    if radius_in_arcsec is not None:
        return elliptical_photometry(fits_file, ra, dec, radius_arcsec, radius_arcsec, 0,
                                   radius_in_arcsec, radius_in_arcsec, include_zero_mask)
    else:
        return elliptical_photometry(fits_file, ra, dec, radius_arcsec, radius_arcsec, 0,
                                   include_zero_mask=include_zero_mask)


def visualize_aperture_region(fits_file, ra, dec, a_arcsec, b_arcsec, theta_deg=0,
                            back_outer_arcsec=None, back_inner_arcsec=None,
                            display_size_arcsec=None, output_filename=None,
                            include_zero_mask=True, show_plot=True, source_name=None):
    """Create a visualization showing source and background regions on image data.

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

    Returns
    -------
    str
        Path to the saved visualization image file.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from photutils.aperture import EllipticalAperture
    import numpy as np
    from datetime import datetime
    import os

    # Open the FITS file
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        # Create bad pixel mask
        if include_zero_mask:
            bad_pixel_mask = np.isnan(data) | (data == 0)
        else:
            bad_pixel_mask = np.isnan(data)

    # Get WCS and convert coordinates
    wcs = WCS(header)
    sky_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    pixel_coord = wcs.world_to_pixel(sky_coord)
    x_center, y_center = pixel_coord

    # Get pixel scale
    pixel_scale = get_pixel_scale(header)

    # Convert source ellipse to pixels
    a_pixels = a_arcsec / pixel_scale
    b_pixels = b_arcsec / pixel_scale
    theta_rad = np.deg2rad(theta_deg)

    # Check if background annulus is provided
    has_background = (back_outer_arcsec is not None) and (back_inner_arcsec is not None)

    if has_background:
        back_outer_pixels = back_outer_arcsec / pixel_scale
        back_inner_pixels = back_inner_arcsec / pixel_scale

    # Determine display size - make sure background annulus fits
    if display_size_arcsec is None:
        if has_background:
            display_size_arcsec = 2.5 * back_outer_arcsec
        else:
            display_size_arcsec = 4 * max(a_arcsec, b_arcsec)

    display_size_pixels = display_size_arcsec / pixel_scale
    half_size = int(display_size_pixels / 2)

    # Define cutout bounds
    x_min = max(0, int(x_center - half_size))
    x_max = min(data.shape[1], int(x_center + half_size))
    y_min = max(0, int(y_center - half_size))
    y_max = min(data.shape[0], int(y_center + half_size))

    # Extract cutout
    cutout_data = data[y_min:y_max, x_min:x_max]
    cutout_mask = bad_pixel_mask[y_min:y_max, x_min:x_max]

    # Adjust center coordinates for cutout
    x_center_cutout = x_center - x_min
    y_center_cutout = y_center - y_min

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Calculate display range for consistent scaling
    valid_data = cutout_data[~np.isnan(cutout_data) & (cutout_data != 0)]
    if len(valid_data) > 0:
        vmin, vmax = np.percentile(valid_data, [1, 99])
    else:
        vmin, vmax = np.nanmin(cutout_data), np.nanmax(cutout_data)

    # Plot 1: Original data with aperture overlay
    im1 = axes[0].imshow(cutout_data, origin='lower', vmin=vmin, vmax=vmax, cmap='viridis')
    axes[0].set_title('Original Data + Apertures')

    # Plot 2: Masked data
    masked_data = cutout_data.copy()
    masked_data[cutout_mask] = np.nan
    im2 = axes[1].imshow(masked_data, origin='lower', vmin=vmin, vmax=vmax, cmap='viridis')
    axes[1].set_title('Data (Masked Pixels Removed)')

    # Plot 3: Source aperture mask
    aperture = EllipticalAperture((x_center_cutout, y_center_cutout),
                                  a=a_pixels, b=b_pixels, theta=theta_rad)

    # Create aperture mask image
    aperture_mask = aperture.to_mask()
    if aperture_mask is not None:
        mask_image = aperture_mask.to_image(cutout_data.shape)
        if mask_image is None:
            mask_image = np.zeros_like(cutout_data)
    else:
        mask_image = np.zeros_like(cutout_data)

    im3 = axes[2].imshow(mask_image, origin='lower', cmap='Reds', alpha=0.7)
    axes[2].imshow(cutout_data, origin='lower', vmin=vmin, vmax=vmax, cmap='viridis', alpha=0.3)
    axes[2].set_title('Source Region (Red Overlay)')

    # Add aperture outlines to all plots
    for ax in axes:
        # Source ellipse in red
        ellipse_src = patches.Ellipse((x_center_cutout, y_center_cutout),
                                      width=2*a_pixels, height=2*b_pixels,
                                      angle=theta_deg, linewidth=2,
                                      edgecolor='red', facecolor='none',
                                      label='Source')
        ax.add_patch(ellipse_src)

        # Background annulus in magenta (if provided)
        if has_background:
            # Outer circle of background annulus
            circle_out = patches.Circle((x_center_cutout, y_center_cutout),
                                        radius=back_outer_pixels, linewidth=2,
                                        edgecolor='magenta', facecolor='none',
                                        linestyle='--', label='Background (outer)')
            ax.add_patch(circle_out)

            # Inner circle of background annulus
            circle_in = patches.Circle((x_center_cutout, y_center_cutout),
                                       radius=back_inner_pixels, linewidth=2,
                                       edgecolor='magenta', facecolor='none',
                                       linestyle='--', label='Background (inner)')
            ax.add_patch(circle_in)

        # Mark center
        ax.plot(x_center_cutout, y_center_cutout, 'w+', markersize=10, markeredgewidth=2)

    # Add colorbars
    plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)

    # Add coordinate labels
    for ax in axes:
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')

    # Add legend to first plot
    legend_elements = [patches.Patch(facecolor='none', edgecolor='red', linewidth=2, label='Source')]
    if has_background:
        legend_elements.append(patches.Patch(facecolor='none', edgecolor='magenta',
                                             linewidth=2, linestyle='--', label='Background'))
    axes[0].legend(handles=legend_elements, loc='upper right', fontsize=8)

    # Add text annotations
    info_text = f"RA: {ra:.6f}°, Dec: {dec:.6f}°\n"
    info_text += f"Center: ({x_center:.1f}, {y_center:.1f}) pixels\n"
    info_text += f"Source: a={a_arcsec:.1f}\", b={b_arcsec:.1f}\", θ={theta_deg:.1f}°"
    if has_background:
        info_text += f"\nBackground: r_out={back_outer_arcsec:.1f}\", r_in={back_inner_arcsec:.1f}\""

    # Title with source name if provided
    if source_name:
        fig.suptitle(f'{source_name}', fontsize=14)
    else:
        fig.suptitle('Source and Background Apertures', fontsize=14)

    fig.text(0.02, 0.02, info_text, fontsize=9, verticalalignment='bottom',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    # Generate filename if not provided
    if output_filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = os.path.splitext(os.path.basename(fits_file))[0]
        output_filename = f"{base_name}_aperture_visualization_{timestamp}.png"

    # Save the plot
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Aperture visualization saved to: {output_filename}")

    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return output_filename


def print_photometry_summary(results):
    """Print a formatted summary of photometry results.

    Parameters
    ----------
    results : dict
        Photometry results dictionary from elliptical_photometry().

    Returns
    -------
    None
        Prints summary to stdout.
    """
    print(f"\n=== Photometry Summary ===")
    print(f"Aperture type: {results['aperture_type']}")
    print(f"Flux: {results['flux']:.6e}")
    print(f"Mean: {results['mean']:.6f}")
    print(f"Median: {results['median']:.6f}")
    print(f"Std: {results['std']:.6f}")
    if results['mode'] is not None:
        print(f"Mode: {results['mode']:.6f}")
    
    print(f"\nArea: {results['area_pixels']:.2f} pixels ({results['area_arcsec2']:.4f} arcsec²)")
    print(f"Surface brightness: {results['surface_brightness_per_arcsec2']:.6e} per arcsec²")
    
    print(f"\nPixel usage: {results['num_pixels_used']}/{results['num_pixels_total']} " +
          f"({results['fraction_pixels_used']:.1%} usable)")
    print(f"Masked pixels: {results['num_pixels_masked']}")
    
    if results['is_annulus']:
        print(f"\nAnnulus parameters:")
        print(f"  Outer: a={results['a_arcsec']:.2f}\", b={results['b_arcsec']:.2f}\"")
        print(f"  Inner: a={results['a_in_arcsec']:.2f}\", b={results['b_in_arcsec']:.2f}\"")
    else:
        print(f"\nEllipse parameters:")
        print(f"  a={results['a_arcsec']:.2f}\", b={results['b_arcsec']:.2f}\"")
    
    print(f"  θ={results['theta_deg']:.1f}°")
    print(f"  Pixel scale: {results['pixel_scale']:.4f}\"/pixel")





def do_one(source_name='SNR_N49', region_file='lmc_spec_ann_reg.txt',
           image_file='../LMC/LMC.ha.csub.fits.gz'):
    """Process a single source from a region file.

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
    """
    xtab = ascii.read(region_file)
    one_source = xtab[xtab['Source_name'] == source_name]

    source_line = one_source[one_source['SourceBack'] == 'Source']
    ra = source_line['RA'][0]
    dec = source_line['Dec'][0]
    a = source_line['Major'][0]
    b = source_line['Minor'][0]
    theta = source_line['Theta'][0]

    sresults = elliptical_region_photometry(fits_file=image_file, ra=ra, dec=dec,
                                            a_arcsec=a, b_arcsec=b, theta_deg=theta)
    x = {'Source_name': source_name, 'SourceBack': 'Source'}
    sresults = x | sresults

    back_line = one_source[one_source['SourceBack'] == 'Back']
    ra = back_line['RA'][0]
    dec = back_line['Dec'][0]
    a = back_line['Major'][0]
    b = back_line['Minor'][0]
    theta = back_line['Theta'][0]

    bresults = circular_photometry(fits_file=image_file, ra=ra, dec=dec,
                                   radius_arcsec=a, radius_in_arcsec=b)
    bresults = {'Source_name': source_name, 'SourceBack': 'Back'} | bresults

    return [sresults, bresults]



def do_many(xtab, image_file, create_visualization=False):
    """Process all regions in a table for a given image.

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
    """
    import os

    xresults = []

    # Build a lookup dictionary for background regions by source name
    back_lookup = {}
    for row in xtab:
        if row['SourceBack'] == 'Back':
            back_lookup[row['Source_name']] = {
                'outer': row['Major'],
                'inner': row['Minor']
            }

    for one in xtab:
        ra = one['RA']
        dec = one['Dec']
        a = one['Major']
        b = one['Minor']
        theta = one['Theta']

        x = {'Source_name': one['Source_name'], 'SourceBack': one['SourceBack']}

        if one['RegType'] in ('ellipse', 'circle'):
            if one['RegType'] == 'circle':
                b = a
                theta = 0
            results = elliptical_region_photometry(fits_file=image_file, ra=ra, dec=dec,
                                                   a_arcsec=a, b_arcsec=b, theta_deg=theta)

            if results is None:
                continue

            # Create visualization for source regions (not background)
            if create_visualization and one['SourceBack'] == 'Source':
                iname = os.path.basename(image_file).replace('.fits', '').replace('.gz', '')
                viz_filename = 'Figs_Flux/%s_%s.png' % (iname, one['Source_name'])

                # Look up the corresponding background region
                back_info = back_lookup.get(one['Source_name'])
                if back_info:
                    back_outer = back_info['outer']
                    back_inner = back_info['inner']
                else:
                    back_outer = None
                    back_inner = None

                visualize_aperture_region(
                    fits_file=image_file, ra=ra, dec=dec,
                    a_arcsec=a, b_arcsec=b, theta_deg=theta,
                    back_outer_arcsec=back_outer, back_inner_arcsec=back_inner,
                    output_filename=viz_filename, show_plot=False,
                    source_name=one['Source_name']
                )

        elif one['RegType'] == 'annulus':
            results = circular_photometry(fits_file=image_file, ra=ra, dec=dec,
                                          radius_arcsec=a, radius_in_arcsec=b)
            if results is None:
                continue
        else:
            print('Error: Unknown RegType:', one['RegType'])
            raise ValueError(f"Unknown RegType: {one['RegType']}")

        results = x | results
        xresults.append(results)

    return xresults
            


def results2table(results_list):
    """Convert a list of result dictionaries to an Astropy Table.

    Parameters
    ----------
    results_list : list of dict
        List of photometry result dictionaries.

    Returns
    -------
    astropy.table.Table
        Table with formatted columns. Float columns use '.2f' format,
        and None values are converted to -99.0.
    """
    table = Table(rows=results_list)

    for col_name in table.colnames:
        if table[col_name].dtype == np.float64:
            table[col_name].format = '.2f'

    for col_name in table.colnames:
        col = table[col_name]

        if col.dtype == object:
            table[col_name] = [float(val) if val is not None else -99.0 for val in col]

        if table[col_name].dtype in [np.float64, np.float32] or np.issubdtype(table[col_name].dtype, np.floating):
            table[col_name].format = '.2f'

    return table



def add_net_rows(xtab):
    """Add background-subtracted Net rows to a photometry table.

    For each source that has both Source and Back entries (matched by
    Source_name and Image), compute a Net row with background-subtracted
    flux and statistics. Rows are ordered as Source, Back, Net for each
    source in each image.

    Parameters
    ----------
    xtab : astropy.table.Table
        Table containing Source and Back entries identified by the
        SourceBack column.

    Returns
    -------
    astropy.table.Table
        Table with rows grouped as Source, Back, Net per source/image.
        Net rows have SourceBack='Net' and contain background-subtracted
        values for flux, mean, and median. The net flux is computed as
        source_flux - num_pixels_used * back_median.
    """
    has_image = 'Image' in xtab.colnames

    # Build lookups keyed by (Source_name, Image)
    src_lookup = {}
    back_lookup = {}
    key_order = []

    for row in xtab:
        if has_image:
            key = (row['Source_name'], row['Image'])
        else:
            key = (row['Source_name'],)

        if key not in src_lookup and key not in back_lookup:
            key_order.append(key)

        if row['SourceBack'] == 'Source':
            src_lookup[key] = row
        elif row['SourceBack'] == 'Back':
            back_lookup[key] = row

    # Build output rows: Source, Back, Net for each key
    out_rows = []
    for key in key_order:
        src_row = src_lookup.get(key)
        back_row = back_lookup.get(key)

        if src_row is not None:
            out_rows.append(dict(src_row))
        if back_row is not None:
            out_rows.append(dict(back_row))

        # Compute Net row if both Source and Back exist
        if src_row is not None and back_row is not None:
            net = dict(src_row)
            net['SourceBack'] = 'Net'
            net['flux'] = src_row['flux'] - src_row['num_pixels_used'] * back_row['median']
            net['mean'] = src_row['mean'] - back_row['mean']
            net['median'] = src_row['median'] - back_row['median']
            if src_row['area_arcsec2'] > 0:
                net['surface_brightness_per_arcsec2'] = net['flux'] / src_row['area_arcsec2']
                net['surface_brightness_per_pixel'] = net['flux'] / src_row['area_pixels']
            else:
                net['surface_brightness_per_arcsec2'] = 0
                net['surface_brightness_per_pixel'] = 0
            out_rows.append(net)

    result = Table(rows=out_rows)
    for col_name in result.colnames:
        if col_name in xtab.colnames and hasattr(xtab[col_name], 'format'):
            result[col_name].format = xtab[col_name].format

    return result



def do_all(image_file, region_table, create_visualization=False):
    """Process all sources in a region table for one image.

    Parameters
    ----------
    image_file : str
        Path to the FITS image file.
    region_table : str
        Path to the region table file.
    create_visualization : bool, optional
        If True, create visualization plots for each source region.
        Plots are saved to ``Figs_Flux/`` directory. Default is False.

    Returns
    -------
    astropy.table.Table or None
        Photometry results table with Source and Back rows and an
        'Image' column, or None if no regions overlap the image.
    """
    import os

    try:
        xtab = ascii.read(region_table)
    except Exception:
        print('Error: could not read %s' % region_table)
        return None

    iname = image_file.split('/')[-1]
    iname = iname.replace('.fits', '')
    iname = iname.replace('.gz', '')

    # Create visualization directory if needed
    if create_visualization:
        os.makedirs('Figs_Flux', exist_ok=True)

    results = do_many(xtab, image_file, create_visualization=create_visualization)

    if len(results) == 0:
        print('  No regions overlap with %s, skipping' % image_file)
        return None

    xresults = results2table(results)
    xresults['Image'] = iname
    print('  %s: %d regions measured' % (iname, len(xresults)))
    return xresults


def steer(argv):
    """Parse command line arguments and execute photometry.

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
    """
    images = []
    reg_file = ''
    create_viz = False
    auto_back = False
    gap = 3.0

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-viz':
            create_viz = True
        elif argv[i] == '-auto_back':
            auto_back = True
        elif argv[i] == '-gap':
            i += 1
            try:
                gap = float(argv[i])
            except (ValueError, IndexError):
                print('Error: -gap requires a numeric value')
                return
        elif argv[i][0] == '-':
            print('Error: unknown option on command line:', argv[i])
            return
        elif argv[i].count('fits'):
            images.append(argv[i])
        elif reg_file == '':
            reg_file = argv[i]
        else:
            print('Cannot parse command line:', argv)
            return
        i += 1

    if len(images) == 0 or reg_file == '':
        print('Error: Not enough arguments')
        print('Usage: GetImageFlux.py [-viz] [-auto_back] [-gap N] image1.fits [image2.fits ...] region_table.txt')
        return

    # If auto_back is enabled, generate the region table with background regions
    if auto_back:
        print(f'Auto-generating background regions with gap={gap:.1f} arcsec')
        generated_file, _ = generate_region_table(reg_file, gap=gap)
        if generated_file is None:
            print('Error: Failed to generate region table')
            return
        reg_file = generated_file

    print('Processing %d images with region file %s' % (len(images), reg_file))
    if create_viz:
        print('Visualization enabled - output to Figs_Flux/')

    all_tables = []
    for one_image in images:
        print('Processing %s' % one_image)
        sb_table = do_all(one_image, reg_file, create_visualization=create_viz)
        if sb_table is not None:
            all_tables.append(sb_table)

    # Write single consolidated file with Source, Back, and Net rows
    if len(all_tables) > 0:
        from astropy.table import vstack

        rname = reg_file.split('/')[-1]
        rname = rname.replace('.txt', '').replace('.tab', '')

        all_results = vstack(all_tables)
        all_results = add_net_rows(all_results)
        outfile = 'Flux_%s.txt' % rname
        all_results.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)

        n_src = np.sum(all_results['SourceBack'] == 'Source')
        n_back = np.sum(all_results['SourceBack'] == 'Back')
        n_net = np.sum(all_results['SourceBack'] == 'Net')
        print('\nWrote %s with %d rows (%d Source, %d Back, %d Net) from %d images' %
              (outfile, len(all_results), n_src, n_back, n_net, len(all_tables)))
    else:
        print('\nNo images had overlapping regions - no output written')

    return







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv) 
    else:
        print (__doc__)
