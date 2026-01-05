#!/usr/bin/env python
# coding: utf-8

# # Create an artificial image given a psf


'''
                    Space Telescope Science Institute

Synopsis:  

Construct an artificial image of a field given 
an original image, a psf, and a list of sources
with fluxes (calculated using aperstats) from
the image.  A residual image is alse constructed


PSF Image Generator Module

Generate artificial images from PSF models and star catalogs.
Completely independent from PSF building - just reads PSF FITS files and star tables.




Command line usage (if any):

    ussage ImageBuilder [-h] [-out root] original_image psf_file star_table

Description:  

Primary routines:

    doit

Notes:
                                       
History:

251216 ksl Coding begun
251223 ksl Added to kred

'''

import sys
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.wcs import WCS
import warnings


class PSFImageGenerator:
    """
    Generate artificial images using PSF models and star catalogs.

    Parameters
    ----------
    psf_file : str
        Path to PSF FITS file (from PSF builder or any compatible PSF)
    reference_image : str, optional
        Path to reference FITS image to copy WCS and dimensions from
    star_table : astropy.Table, optional
        Table with star positions and fluxes (xcenter, ycenter, Net columns)
    """

    def __init__(self, psf_file, reference_image=None, star_table=None):
        # Load PSF
        with fits.open(psf_file) as hdul:
            self.psf = hdul[0].data.astype(np.float64)
            self.psf_header = hdul[0].header

        self.psf_file = psf_file
        self.psf_size = self.psf.shape[0]
        self.psf_half = self.psf_size // 2

        # Load reference image if provided
        self.ref_header = None
        self.wcs = None
        self.image_shape = None

        if reference_image is not None:
            self.load_reference_image(reference_image)

        # Load star table if provided
        self.star_table = star_table

        print(f"Loaded PSF from {psf_file}")
        print(f"PSF shape: {self.psf.shape}")
        print(f"PSF peak: {np.max(self.psf):.6f}")
        print(f"PSF sum: {np.sum(self.psf):.6f}")

    def load_reference_image(self, reference_image):
        """
        Load reference image to get WCS and dimensions.

        Parameters
        ----------
        reference_image : str
            Path to reference FITS image
        """
        with fits.open(reference_image) as hdul:
            self.ref_header = hdul[0].header.copy()
            self.image_shape = hdul[0].data.shape

            # Try to load WCS
            try:
                self.wcs = WCS(self.ref_header)
                print(f"Loaded WCS from {reference_image}")
            except Exception as e:
                warnings.warn(f"Could not load WCS: {e}")
                self.wcs = None

        print(f"Reference image shape: {self.image_shape}")

    def load_star_table(self, star_table):
        """
        Load star table with positions and fluxes.

        Parameters
        ----------
        star_table : str or astropy.Table
            Path to FITS table or Table object with xcenter, ycenter, Net columns
        """
        if isinstance(star_table, str):
            self.star_table = Table.read(star_table)
        else:
            self.star_table = star_table

        print(f"Loaded {len(self.star_table)} stars from catalog")

        # Verify required columns
        required = ['xcenter', 'ycenter', 'Net']
        missing = [col for col in required if col not in self.star_table.colnames]
        if missing:
            print('Trying alternative column names')
            try:
                self.star_table.rename_column('flux_fit','Net')
                self.star_table.rename_column('x_fit','xcenter')
                self.star_table.rename_column('y_fit','ycenter')
            except:
                raise ValueError(f"Star table missing required columns: {missing}")

    def add_star_to_image(self, image, x, y, flux, subpixel=True):
        """
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
        """
        ny, nx = image.shape

        # Integer pixel position
        ix = int(np.round(x))
        iy = int(np.round(y))

        # Calculate PSF placement bounds
        x_start = ix - self.psf_half
        x_end = ix + self.psf_half + 1
        y_start = iy - self.psf_half
        y_end = iy + self.psf_half + 1

        # Calculate overlap with image
        img_x_start = max(0, x_start)
        img_x_end = min(nx, x_end)
        img_y_start = max(0, y_start)
        img_y_end = min(ny, y_end)

        # Calculate corresponding PSF region
        psf_x_start = img_x_start - x_start
        psf_x_end = self.psf_size - (x_end - img_x_end)
        psf_y_start = img_y_start - y_start
        psf_y_end = self.psf_size - (y_end - img_y_end)

        # Check if star is within image bounds
        if img_x_end <= img_x_start or img_y_end <= img_y_start:
            return  # Star completely outside image

        # Get PSF section
        psf_section = self.psf[psf_y_start:psf_y_end, psf_x_start:psf_x_end].copy()

        # Apply subpixel shift if requested
        if subpixel:
            dx = x - ix
            dy = y - iy

            if abs(dx) > 0.01 or abs(dy) > 0.01:
                # Simple bilinear interpolation for subpixel shift
                # Create shifted PSF by weighted sum of neighboring pixels
                psf_section = self._shift_psf_subpixel(psf_section, dx, dy)

        # Scale PSF by flux and add to image
        # If PSF is already normalized to sum=1, flux is directly the total flux
        # If PSF is normalized to peak=1, we need to scale appropriately
        norm_method = self.psf_header.get('PSFNORM', 'peak')

        if norm_method == 'sum':
            # PSF sums to 1, so flux is the scaling factor
            scaled_psf = psf_section * flux
        elif norm_method == 'peak':
            # PSF peaks at 1, scale by flux/sum to preserve total flux
            psf_sum = np.sum(self.psf)
            scaled_psf = psf_section * (flux / psf_sum)
        else:
            # Unknown normalization, assume sum
            psf_sum = np.sum(self.psf)
            scaled_psf = psf_section * (flux / psf_sum)

        # Add to image
        image[img_y_start:img_y_end, img_x_start:img_x_end] += scaled_psf

    def _shift_psf_subpixel(self, psf, dx, dy):
        """
        Apply subpixel shift to PSF using bilinear interpolation.

        Parameters
        ----------
        psf : ndarray
            PSF array
        dx : float
            X shift in pixels (-0.5 to 0.5)
        dy : float
            Y shift in pixels (-0.5 to 0.5)

        Returns
        -------
        shifted : ndarray
            Shifted PSF
        """
        # Simple implementation using roll and weighted average
        # For better accuracy, could use scipy.ndimage.shift

        # Limit shift to reasonable range
        dx = np.clip(dx, -0.5, 0.5)
        dy = np.clip(dy, -0.5, 0.5)

        # Create shifted versions
        shifted = psf.copy()

        if abs(dx) > 0.01:
            # Shift in x
            if dx > 0:
                shifted = (1 - dx) * shifted + dx * np.roll(shifted, 1, axis=1)
            else:
                shifted = (1 + dx) * shifted - dx * np.roll(shifted, -1, axis=1)

        if abs(dy) > 0.01:
            # Shift in y
            if dy > 0:
                shifted = (1 - dy) * shifted + dy * np.roll(shifted, 1, axis=0)
            else:
                shifted = (1 + dy) * shifted - dy * np.roll(shifted, -1, axis=0)

        return shifted

    def generate_image(self, image_shape=None, background=0.0, subpixel=True, 
                      star_table=None, flux_column='Net'):
        """
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
        """
        # Determine image shape
        if image_shape is None:
            if self.image_shape is None:
                raise ValueError("Must provide image_shape or load reference image")
            image_shape = self.image_shape

        # Use provided star table or default
        if star_table is None:
            if self.star_table is None:
                raise ValueError("Must provide star_table or load one")
            star_table = self.star_table
        else:
            star_table = star_table

        # Verify flux column exists
        if flux_column not in star_table.colnames:
            raise ValueError(f"Flux column '{flux_column}' not found in star table. "
                           f"Available columns: {star_table.colnames}")

        # Create blank image
        image = np.full(image_shape, background, dtype=np.float64)

        # Add each star
        n_stars = len(star_table)
        n_added = 0
        total_flux = 0.0

        print(f"Generating image with {n_stars} stars...")

        for i, star in enumerate(star_table):
            x = star['xcenter']
            y = star['ycenter']
            flux = star[flux_column]

            # Skip stars with invalid flux
            if not np.isfinite(flux) or flux <= 0:
                continue

            self.add_star_to_image(image, x, y, flux, subpixel=subpixel)
            n_added += 1
            total_flux += flux

            # Progress indicator
            if (i + 1) % 100 == 0 or (i + 1) == n_stars:
                print(f"  Added {i + 1}/{n_stars} stars", end='\r')

        print(f"\nAdded {n_added} stars with total flux {total_flux:.2f}")

        return image

    def save_image(self, image, output_file, copy_header=True):
        """
        Save generated image to FITS file.

        Parameters
        ----------
        image : ndarray
            Image array to save
        output_file : str
            Output FITS filename
        copy_header : bool, optional
            Copy full header from reference image if available (default: True)
        """
        # Create FITS HDU with full header if available
        if copy_header and self.ref_header is not None:
            # Start with complete reference header
            hdu = fits.PrimaryHDU(image, header=self.ref_header.copy())
        else:
            # Create minimal header
            hdu = fits.PrimaryHDU(image)

        # Add/update metadata about the artificial image
        hdu.header['PSFFILE'] = (self.psf_file, 'PSF model used')
        hdu.header['PSFTYPE'] = (self.psf_header.get('PSFTYPE', 'UNKNOWN'), 'PSF model type')
        hdu.header['PSFNORM'] = (self.psf_header.get('PSFNORM', 'UNKNOWN'), 'PSF normalization')
        hdu.header['ARTIMG'] = (True, 'Artificial image generated from PSF')
        hdu.header['NSTARS'] = (len(self.star_table) if self.star_table is not None else 0, 
                               'Number of stars in catalog')

        # Add history
        hdu.header.add_history(f'Artificial image generated from PSF: {self.psf_file}')
        hdu.header.add_history(f'PSF type: {self.psf_header.get("PSFTYPE", "UNKNOWN")}')
        if self.star_table is not None:
            hdu.header.add_history(f'Number of stars: {len(self.star_table)}')

        # Save
        hdu.writeto(output_file, overwrite=True)
        print(f"Saved artificial image to {output_file}")


def generate_artificial_image(psf_file, star_table, reference_image, 
                              output_file=None, background=0.0, 
                              subpixel=True, flux_column='Net'):
    """
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
    flux_column : str, optional
        Column name for flux values (default: 'Net')

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
    """
    # Create generator
    generator = PSFImageGenerator(psf_file, reference_image=reference_image)

    # Load star table
    generator.load_star_table(star_table)

    # Generate image
    image = generator.generate_image(background=background, 
                                     subpixel=subpixel,
                                     flux_column=flux_column)

    # Save if output file specified
    if output_file is not None:
        generator.save_image(image, output_file)

    return image


def generate_residual_image(original_image, artificial_image, output_file=None):
    """
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
    """
    # Load images
    if isinstance(original_image, str):
        with fits.open(original_image) as hdul:
            orig = hdul[0].data
            orig_header = hdul[0].header.copy()
    else:
        orig = original_image
        orig_header = None

    if isinstance(artificial_image, str):
        with fits.open(artificial_image) as hdul:
            art = hdul[0].data
    else:
        art = artificial_image

    # Calculate residual
    residual = orig - art

    # Save if requested
    if output_file is not None:
        if orig_header is not None:
            # Use complete original header
            hdu = fits.PrimaryHDU(residual, header=orig_header)
        else:
            hdu = fits.PrimaryHDU(residual)

        # Add/update metadata
        hdu.header['RESIDUAL'] = (True, 'Residual image (original - model)')
        hdu.header.add_history('Residual image: original - artificial model')

        hdu.writeto(output_file, overwrite=True)
        print(f"Saved residual image to {output_file}")

    return residual


def doit(psf_file='my_psf_gaussian.fits',star_table='forced.fits',reference_image='data/SMC_c01_T08.N673.fits',root='test'):

    if root=='':
        root='test'

    # Method 1: Quick generation
    print("=== Quick Method ===")
    artificial = generate_artificial_image(
        psf_file=psf_file,
        star_table=star_table,
        reference_image=reference_image,
        output_file='%s_artificial.fits' % root,
        background=0.0,
        subpixel=True
        )

    # Method 2: More control
    # print("\n=== Advanced Method ===")

    # Load components
    # generator = PSFImageGenerator('psf_summed.fits', reference_image='science.fits')
    # generator.load_star_table('stars.fits')

    # Generate with custom settings
    # artificial = generator.generate_image(background=100.0, subpixel=True)
    #generator.save_image(artificial, 'artificial_summed.fits')

    # Generate residual
    print("\n=== Generating Residual ===")
    residual = generate_residual_image(
        reference_image, '%s_artificial.fits' % root,
        output_file='%s_residual.fits' % root
    )

    print(f"\nResidual statistics:")
    print(f"  Mean: {np.mean(residual):.2f}")
    print(f"  Std: {np.std(residual):.2f}")
    print(f"  Min: {np.min(residual):.2f}")
    print(f"  Max: {np.max(residual):.2f}")
    return



def steer(argv):
    '''
    ussage ImageBuilder [-h] [-out] original_image psf_file star_table
    '''

    ref_image=''
    psf_file=''
    star_table=''
    root=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            root=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown option :',argv)
        elif ref_image=='':
            ref_image=argv[i]
        elif psf_file=='':
            psf_file=argv[i]
        elif star_table=='':
            star_table=argv[i]
        else:
            print('Error: To many commands: ',argv)
        i+=1

    if root=='':
        root=ref_image.split('/')[-1]
        root=root.replace('.fits','')

    print('Reference image : %s' % ref_image)
    print('       PSF file : %s' % psf_file)
    print('     Star table : %s' % star_table)
    print('       Out root : %s' % root)


    doit(psf_file=psf_file,star_table=star_table,
         reference_image=ref_image,root=root)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
