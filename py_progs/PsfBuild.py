#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Build PSF from a selection of stars


Command line usage (if any):

    Usage: PsFBuild [-out root] image_file psf_stars

Description:  

Primary routines:

    doit

Notes:
                                       
History:

251210 ksl Coding begun
251222 ksl Added to kred

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt

"""
PSF Builder Module

Build Point Spread Functions from FITS images using star positions from astropy tables.
Supports Elliptical Gaussian, Elliptical Moffat, and empirical summed PSF methods.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import minimize
import warnings


class PSFBuilder:
    """
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
    """

    def __init__(self, fits_file, star_table, stamp_size=25, normalize='peak'):
        # Load FITS data
        if isinstance(fits_file, str):
            with fits.open(fits_file) as hdul:
                self.data = hdul[0].data.astype(np.float64)
                self.header = hdul[0].header
        else:
            self.data = fits_file[0].data.astype(np.float64)
            self.header = fits_file[0].header

        self.star_table = star_table
        self.stamp_size = stamp_size
        self.normalize = normalize
        self.stamps = None
        self.psf_gaussian = None
        self.psf_moffat = None
        self.psf_summed = None

    def extract_stamps(self, filter_outliers=True, sigma_clip=3.0):
        """
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
        """
        stamps = []
        half = self.stamp_size // 2
        ny, nx = self.data.shape

        for row in self.star_table:
            xc = int(np.round(row['xcenter']))
            yc = int(np.round(row['ycenter']))

            # Check boundaries
            if (xc - half < 0 or xc + half + 1 > nx or 
                yc - half < 0 or yc + half + 1 > ny):
                continue

            # Extract stamp
            stamp_data = self.data[yc-half:yc+half+1, xc-half:xc+half+1].copy()

            # Skip if stamp contains NaN or Inf
            if not np.all(np.isfinite(stamp_data)):
                continue

            stamps.append({
                'data': stamp_data,
                'xcenter': xc,
                'ycenter': yc,
                'x_offset': row['xcenter'] - xc,
                'y_offset': row['ycenter'] - yc,
                'peak': np.max(stamp_data)
            })

        # Filter outliers based on peak flux
        if filter_outliers and len(stamps) > 3:
            peaks = np.array([s['peak'] for s in stamps])
            median_peak = np.median(peaks)
            std_peak = np.std(peaks)

            stamps = [s for s in stamps 
                     if np.abs(s['peak'] - median_peak) < sigma_clip * std_peak]

        self.stamps = stamps
        print(f"Extracted {len(stamps)} valid stamps from {len(self.star_table)} stars")
        return stamps

    def _normalize_psf(self, psf_data, method=None):
        """
        Normalize PSF according to specified method.

        Parameters
        ----------
        psf_data : ndarray
            PSF array to normalize
        method : str, optional
            'peak', 'sum', or 'none'. If None, uses self.normalize

        Returns
        -------
        normalized : ndarray
            Normalized PSF
        """
        if method is None:
            method = self.normalize

        if method == 'peak':
            return psf_data / np.max(psf_data)
        elif method == 'sum':
            return psf_data / np.sum(psf_data)
        elif method == 'none':
            return psf_data
        else:
            raise ValueError(f"Unknown normalization method: {method}")

    def build_gaussian_psf(self):
        """
        Build elliptical Gaussian PSF model.

        Returns
        -------
        params : dict
            Fitted parameters: amplitude, x0, y0, sigma_x, sigma_y, theta, fwhm_x, fwhm_y
        model : ndarray
            2D array of the fitted PSF model (normalized)
        """
        if self.stamps is None:
            self.extract_stamps()

        # Prepare data for fitting
        all_data = []
        half = self.stamp_size // 2

        for stamp in self.stamps:
            for i in range(self.stamp_size):
                for j in range(self.stamp_size):
                    x = j - half - stamp['x_offset']
                    y = i - half - stamp['y_offset']
                    val = stamp['data'][i, j]
                    all_data.append([x, y, val])

        all_data = np.array(all_data)

        # Initial parameter estimates
        peak = np.max(all_data[:, 2])
        background = np.percentile(all_data[:, 2], 10)

        # Initial guess
        p0 = [peak - background, 0.0, 0.0, 2.0, 2.0, 0.0, background]

        # Fit
        def gaussian_2d(params, x, y):
            amp, x0, y0, sx, sy, theta, bg = params
            cos_t = np.cos(theta)
            sin_t = np.sin(theta)
            xr = (x - x0) * cos_t + (y - y0) * sin_t
            yr = -(x - x0) * sin_t + (y - y0) * cos_t
            return amp * np.exp(-0.5 * (xr**2 / sx**2 + yr**2 / sy**2)) + bg

        def residuals(params):
            model = gaussian_2d(params, all_data[:, 0], all_data[:, 1])
            return np.sum((all_data[:, 2] - model)**2)

        # Bounds
        bounds = [(0, None), (-2, 2), (-2, 2), (0.5, 10), (0.5, 10), 
                  (-np.pi, np.pi), (0, None)]

        result = minimize(residuals, p0, bounds=bounds, method='L-BFGS-B')

        params = {
            'amplitude': result.x[0],
            'x0': result.x[1],
            'y0': result.x[2],
            'sigma_x': result.x[3],
            'sigma_y': result.x[4],
            'theta': result.x[5],
            'background': result.x[6],
            'fwhm_x': 2.355 * result.x[3],
            'fwhm_y': 2.355 * result.x[4]
        }

        # Create model PSF (without background)
        y, x = np.mgrid[-half:half+1, -half:half+1]
        model_raw = result.x[0] * np.exp(-0.5 * (
            ((x - result.x[1]) * np.cos(result.x[5]) + (y - result.x[2]) * np.sin(result.x[5]))**2 / result.x[3]**2 +
            (-(x - result.x[1]) * np.sin(result.x[5]) + (y - result.x[2]) * np.cos(result.x[5]))**2 / result.x[4]**2
        ))

        # Normalize
        model = self._normalize_psf(model_raw)

        self.psf_gaussian = {'params': params, 'model': model}
        return params, model

    def build_moffat_psf(self):
        """
        Build elliptical Moffat PSF model.

        Returns
        -------
        params : dict
            Fitted parameters: amplitude, alpha, beta, ellipticity, theta, fwhm
        model : ndarray
            2D array of the fitted PSF model (normalized)
        """
        if self.stamps is None:
            self.extract_stamps()

        # Prepare data
        all_data = []
        half = self.stamp_size // 2

        for stamp in self.stamps:
            for i in range(self.stamp_size):
                for j in range(self.stamp_size):
                    x = j - half - stamp['x_offset']
                    y = i - half - stamp['y_offset']
                    val = stamp['data'][i, j]
                    all_data.append([x, y, val])

        all_data = np.array(all_data)

        peak = np.max(all_data[:, 2])
        background = np.percentile(all_data[:, 2], 10)

        # Initial guess: [amplitude, alpha, beta, ellipticity, theta, background]
        p0 = [peak - background, 2.5, 2.5, 0.1, 0.0, background]

        def moffat_2d(params, x, y):
            amp, alpha, beta, ell, theta, bg = params
            cos_t = np.cos(theta)
            sin_t = np.sin(theta)
            xr = (x * cos_t + y * sin_t) / (1 - ell)
            yr = (-x * sin_t + y * cos_t)
            r2 = xr**2 + yr**2
            return amp / (1 + r2 / alpha**2)**beta + bg

        def residuals(params):
            model = moffat_2d(params, all_data[:, 0], all_data[:, 1])
            return np.sum((all_data[:, 2] - model)**2)

        bounds = [(0, None), (0.5, 10), (1, 10), (0, 0.9), 
                  (-np.pi, np.pi), (0, None)]

        result = minimize(residuals, p0, bounds=bounds, method='L-BFGS-B')

        # Calculate FWHM
        alpha = result.x[1]
        beta = result.x[2]
        fwhm = 2 * alpha * np.sqrt(2**(1/beta) - 1)

        params = {
            'amplitude': result.x[0],
            'alpha': alpha,
            'beta': beta,
            'ellipticity': result.x[3],
            'theta': result.x[4],
            'background': result.x[5],
            'fwhm': fwhm
        }

        # Create model PSF (without background)
        y, x = np.mgrid[-half:half+1, -half:half+1]
        cos_t = np.cos(result.x[4])
        sin_t = np.sin(result.x[4])
        ell = result.x[3]
        xr = (x * cos_t + y * sin_t) / (1 - ell)
        yr = (-x * sin_t + y * cos_t)
        r2 = xr**2 + yr**2
        model_raw = result.x[0] / (1 + r2 / result.x[1]**2)**result.x[2]

        # Normalize
        model = self._normalize_psf(model_raw)

        self.psf_moffat = {'params': params, 'model': model}
        return params, model

    def build_summed_psf(self, subtract_background=True):
        """
        Build empirical PSF by summing star stamps.

        Parameters
        ----------
        subtract_background : bool, optional
            Subtract local background from each stamp (default: True)

        Returns
        -------
        psf : ndarray
            2D array of the summed PSF (normalized)
        """
        if self.stamps is None:
            self.extract_stamps()

        summed = np.zeros((self.stamp_size, self.stamp_size))

        for stamp in self.stamps:
            stamp_data = stamp['data'].copy()

            if subtract_background:
                # Estimate background from corners
                corners = np.concatenate([
                    stamp_data[:3, :3].flatten(),
                    stamp_data[:3, -3:].flatten(),
                    stamp_data[-3:, :3].flatten(),
                    stamp_data[-3:, -3:].flatten()
                ])
                bg = np.median(corners)
                stamp_data -= bg
                stamp_data = np.maximum(stamp_data, 0)  # Remove negative values

            summed += stamp_data

        # Normalize
        summed = self._normalize_psf(summed)

        self.psf_summed = summed
        return summed

    def build_all_psfs(self):
        """
        Build all three PSF models.

        Returns
        -------
        results : dict
            Dictionary containing all PSF models and parameters
        """
        print("Building Elliptical Gaussian PSF...")
        gauss_params, gauss_model = self.build_gaussian_psf()

        print("Building Elliptical Moffat PSF...")
        moffat_params, moffat_model = self.build_moffat_psf()

        print("Building Summed PSF...")
        summed_psf = self.build_summed_psf()

        return {
            'gaussian': {'params': gauss_params, 'model': gauss_model},
            'moffat': {'params': moffat_params, 'model': moffat_model},
            'summed': summed_psf,
            'n_stars': len(self.stamps),
            'normalization': self.normalize
        }

    def save_psfs(self, output_prefix='psf'):
        """
        Save PSF models to FITS files.

        Parameters
        ----------
        output_prefix : str, optional
            Prefix for output files (default: 'psf')
        """
        if self.psf_gaussian is not None:
            hdu = fits.PrimaryHDU(self.psf_gaussian['model'])
            # Use FITS-compliant 8-character keywords
            hdu.header['AMP'] = (self.psf_gaussian['params']['amplitude'], 'Fitted amplitude')
            hdu.header['X0'] = (self.psf_gaussian['params']['x0'], 'X centroid offset')
            hdu.header['Y0'] = (self.psf_gaussian['params']['y0'], 'Y centroid offset')
            hdu.header['SIGMAX'] = (self.psf_gaussian['params']['sigma_x'], 'Sigma X (pixels)')
            hdu.header['SIGMAY'] = (self.psf_gaussian['params']['sigma_y'], 'Sigma Y (pixels)')
            hdu.header['THETA'] = (self.psf_gaussian['params']['theta'], 'Position angle (radians)')
            hdu.header['BKGD'] = (self.psf_gaussian['params']['background'], 'Background level')
            hdu.header['FWHMX'] = (self.psf_gaussian['params']['fwhm_x'], 'FWHM X (pixels)')
            hdu.header['FWHMY'] = (self.psf_gaussian['params']['fwhm_y'], 'FWHM Y (pixels)')
            hdu.header['PSFTYPE'] = ('GAUSSIAN', 'PSF model type')
            hdu.header['PSFNORM'] = (self.normalize, 'Normalization method')
            hdu.writeto(f'{output_prefix}_gaussian_psf.fits', overwrite=True)
            print(f"Saved {output_prefix}_gaussian_psf.fits")

        if self.psf_moffat is not None:
            hdu = fits.PrimaryHDU(self.psf_moffat['model'])
            hdu.header['AMP'] = (self.psf_moffat['params']['amplitude'], 'Fitted amplitude')
            hdu.header['ALPHA'] = (self.psf_moffat['params']['alpha'], 'Moffat alpha parameter')
            hdu.header['BETA'] = (self.psf_moffat['params']['beta'], 'Moffat beta parameter')
            hdu.header['ELLIP'] = (self.psf_moffat['params']['ellipticity'], 'Ellipticity')
            hdu.header['THETA'] = (self.psf_moffat['params']['theta'], 'Position angle (radians)')
            hdu.header['BKGD'] = (self.psf_moffat['params']['background'], 'Background level')
            hdu.header['FWHM'] = (self.psf_moffat['params']['fwhm'], 'FWHM (pixels)')
            hdu.header['PSFTYPE'] = ('MOFFAT', 'PSF model type')
            hdu.header['PSFNORM'] = (self.normalize, 'Normalization method')
            hdu.writeto(f'{output_prefix}_moffat_psf.fits', overwrite=True)
            print(f"Saved {output_prefix}_moffat_psf.fits")

        if self.psf_summed is not None:
            hdu = fits.PrimaryHDU(self.psf_summed)
            hdu.header['NSTARS'] = (len(self.stamps), 'Number of stars used')
            hdu.header['PSFTYPE'] = ('SUMMED', 'PSF model type')
            hdu.header['PSFNORM'] = (self.normalize, 'Normalization method')
            hdu.writeto(f'{output_prefix}_summed_psf.fits', overwrite=True)
            print(f"Saved {output_prefix}_summed_psf.fits")


def quick_psf_build(fits_file, star_table, stamp_size=25, normalize='peak', output_prefix=None):
    """
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
    """
    builder = PSFBuilder(fits_file, star_table, stamp_size=stamp_size, normalize=normalize)
    results = builder.build_all_psfs()

    if output_prefix is not None:
        builder.save_psfs(output_prefix)

    return results



def do_one(image_file,psf_stars,prefix=''):


    if prefix=='':
        root=image_file.split('/')[-1]
        root=root.replace('.fits','')
        root=root.replace('.fz','')
        prefix=root


    # Load your star table (from aperture photometry)
    star_table = Table.read(psf_stars)

    # Build PSFs with consistent normalization
    results = quick_psf_build(image_file, star_table, 
                             stamp_size=25, 
                             normalize='peak',  # or 'sum' or 'none'
                             output_prefix=prefix)

    # Access results
    print("\n=== Elliptical Gaussian PSF ===")
    for key, val in results['gaussian']['params'].items():
        print(f"{key}: {val:.4f}")

    print("\n=== Elliptical Moffat PSF ===")
    for key, val in results['moffat']['params'].items():
        print(f"{key}: {val:.4f}")

    print(f"\n=== Summed PSF ===")
    print(f"Shape: {results['summed'].shape}")
    print(f"Peak value: {np.max(results['summed']):.6f}")
    print(f"Sum: {np.sum(results['summed']):.6f}")
    print(f"Built from {results['n_stars']} stars")
    print(f"Normalization: {results['normalization']}")




def steer(argv):
    '''
    This is generally just a steering routine

    Usage: PsFBuild [-out root] image_file psf_stars
    '''

    image_file=''
    star_file=''
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
            print('Error: Could not intepret commands: ',argv)
            return
        elif image_file=='':
            image_file=argv[i]
        elif star_file=='':
            star_file=argv[i]
        else:
            print('Error: Two many commands: ',argv)
            return

        i+=1

    if root=='':
        root=image_file.split('/')[-1]
        root=root.replace('.fits','')
        root=root.replace('.gz','')


    print('     Image: ',image_file)
    print(' PSF_stars: ',star_file)
    print('     Root : ',root)



    do_one(image_file=image_file,psf_stars=star_file,prefix=root)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
