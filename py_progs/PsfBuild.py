#!/usr/bin/env python
# coding: utf-8

"""Build PSF from star catalog

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

"""


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
from scipy.stats import rankdata
import warnings


def select_psf_stars(phot_table, max_stars=100, weights=None, verbose=True):
    """
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
    """

    if weights is None:
        weights = {'snr': 0.35, 'fwhm': 0.25, 'ecc': 0.20, 'bkg': 0.10, 'conc': 0.10}

    # Reference thresholds for diagnostics (not used for selection)
    ref_thresholds = {
        'snr_min': 20,
        'fwhm_tolerance': 0.3,
        'ecc_max': 0.2,
        'bkg_contam_max': 1.5,
        'concentration_min': 3.0
    }

    # First filter to valid stars only
    valid_mask = (
        np.isfinite(phot_table['SNR']) &
        np.isfinite(phot_table['FWHM']) &
        np.isfinite(phot_table['Eccentricity']) &
        np.isfinite(phot_table['BkgContam']) &
        np.isfinite(phot_table['Concentration']) &
        (phot_table['SNR'] > 0) &
        (phot_table['FWHM'] > 0)
    )

    n_total = len(phot_table)
    n_valid = np.sum(valid_mask)

    if n_valid == 0:
        print("Error: No valid stars available for PSF construction")
        return phot_table[valid_mask]  # Return empty table

    candidates = phot_table[valid_mask].copy()

    # Extract metrics
    snr = np.array(candidates['SNR'])
    fwhm = np.array(candidates['FWHM'])
    ecc = np.array(candidates['Eccentricity'])
    bkg = np.array(candidates['BkgContam'])
    conc = np.array(candidates['Concentration'])

    # Calculate FWHM deviation from median
    median_fwhm = np.median(fwhm)
    fwhm_dev = np.abs(fwhm - median_fwhm) / median_fwhm

    # Calculate percentile scores (0-1, higher is better)
    n = len(candidates)

    # SNR: higher is better (use log scale for ranking)
    snr_rank = rankdata(np.log10(snr)) / n

    # FWHM deviation: lower is better
    fwhm_rank = 1.0 - rankdata(fwhm_dev) / n

    # Eccentricity: lower is better (more circular)
    ecc_rank = 1.0 - rankdata(ecc) / n

    # Background contamination: lower is better
    bkg_rank = 1.0 - rankdata(bkg) / n

    # Concentration: higher is better (more point-like)
    conc_rank = rankdata(conc) / n

    # Calculate weighted quality score
    quality_score = (
        weights['snr'] * snr_rank +
        weights['fwhm'] * fwhm_rank +
        weights['ecc'] * ecc_rank +
        weights['bkg'] * bkg_rank +
        weights['conc'] * conc_rank
    )

    # Add quality score to table
    candidates['QualityScore'] = quality_score

    # Print diagnostic statistics
    if verbose:
        print("\n=== PSF Star Selection Diagnostics ===")
        print(f"Total sources: {n_total}")
        print(f"Valid sources (finite values, SNR>0): {n_valid}")

        # Calculate how many pass each threshold
        pass_snr = np.sum(snr >= ref_thresholds['snr_min'])
        pass_fwhm = np.sum(fwhm_dev <= ref_thresholds['fwhm_tolerance'])
        pass_ecc = np.sum(ecc <= ref_thresholds['ecc_max'])
        pass_bkg = np.sum(bkg <= ref_thresholds['bkg_contam_max'])
        pass_conc = np.sum(conc >= ref_thresholds['concentration_min'])
        pass_all = np.sum(
            (snr >= ref_thresholds['snr_min']) &
            (fwhm_dev <= ref_thresholds['fwhm_tolerance']) &
            (ecc <= ref_thresholds['ecc_max']) &
            (bkg <= ref_thresholds['bkg_contam_max']) &
            (conc >= ref_thresholds['concentration_min'])
        )

        print("\nMetric distributions and reference threshold pass rates:")
        print(f"  {'Metric':<15} {'Min':>10} {'Median':>10} {'Max':>10} {'Pass Ref':>12}")
        print(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
        print(f"  {'SNR':<15} {np.min(snr):>10.1f} {np.median(snr):>10.1f} {np.max(snr):>10.1f} {pass_snr:>5}/{n_valid} ({100*pass_snr/n_valid:>4.0f}%)")
        print(f"  {'FWHM':<15} {np.min(fwhm):>10.2f} {np.median(fwhm):>10.2f} {np.max(fwhm):>10.2f} {'--':>12}")
        print(f"  {'FWHM deviation':<15} {np.min(fwhm_dev):>10.2f} {np.median(fwhm_dev):>10.2f} {np.max(fwhm_dev):>10.2f} {pass_fwhm:>5}/{n_valid} ({100*pass_fwhm/n_valid:>4.0f}%)")
        print(f"  {'Eccentricity':<15} {np.min(ecc):>10.3f} {np.median(ecc):>10.3f} {np.max(ecc):>10.3f} {pass_ecc:>5}/{n_valid} ({100*pass_ecc/n_valid:>4.0f}%)")
        print(f"  {'BkgContam':<15} {np.min(bkg):>10.2f} {np.median(bkg):>10.2f} {np.max(bkg):>10.2f} {pass_bkg:>5}/{n_valid} ({100*pass_bkg/n_valid:>4.0f}%)")
        print(f"  {'Concentration':<15} {np.min(conc):>10.2f} {np.median(conc):>10.2f} {np.max(conc):>10.2f} {pass_conc:>5}/{n_valid} ({100*pass_conc/n_valid:>4.0f}%)")
        print(f"\n  Stars passing ALL reference thresholds: {pass_all}/{n_valid} ({100*pass_all/n_valid:.0f}%)")

        # Identify the most problematic metric(s)
        pass_rates = {
            'SNR': pass_snr / n_valid,
            'FWHM deviation': pass_fwhm / n_valid,
            'Eccentricity': pass_ecc / n_valid,
            'BkgContam': pass_bkg / n_valid,
            'Concentration': pass_conc / n_valid
        }
        sorted_metrics = sorted(pass_rates.items(), key=lambda x: x[1])

        if pass_all == 0:
            print("\n  ** Warning: No stars pass all reference thresholds **")
            print("  Most problematic metrics (lowest pass rates):")
            for metric, rate in sorted_metrics[:3]:
                if rate < 0.5:
                    print(f"    - {metric}: only {100*rate:.0f}% pass")

    # Sort by quality score and select top stars
    candidates.sort('QualityScore', reverse=True)
    psf_stars = candidates[:max_stars]

    if verbose:
        print(f"\nSelected {len(psf_stars)} PSF stars by quality score")
        print(f"  Weights: snr={weights['snr']}, fwhm={weights['fwhm']}, ecc={weights['ecc']}, bkg={weights['bkg']}, conc={weights['conc']}")
        print(f"  Quality score range: {psf_stars['QualityScore'][-1]:.3f} - {psf_stars['QualityScore'][0]:.3f}")

        # Extract selected star metrics
        sel_snr = np.array(psf_stars['SNR'])
        sel_fwhm = np.array(psf_stars['FWHM'])
        sel_ecc = np.array(psf_stars['Eccentricity'])
        sel_bkg = np.array(psf_stars['BkgContam'])
        sel_conc = np.array(psf_stars['Concentration'])
        sel_fwhm_dev = np.abs(sel_fwhm - median_fwhm) / median_fwhm

        # Compare selected vs full sample
        print("\nSelected stars vs full sample (median values):")
        print(f"  {'Metric':<15} {'Selected':>12} {'Full Sample':>12} {'Improvement':>12}")
        print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*12}")

        # SNR - higher is better
        sel_med = np.median(sel_snr)
        full_med = np.median(snr)
        improvement = (sel_med - full_med) / full_med * 100 if full_med > 0 else 0
        print(f"  {'SNR':<15} {sel_med:>12.1f} {full_med:>12.1f} {improvement:>+11.0f}%")

        # FWHM - report absolute value
        sel_med = np.median(sel_fwhm)
        full_med = np.median(fwhm)
        print(f"  {'FWHM':<15} {sel_med:>12.2f} {full_med:>12.2f} {'--':>12}")

        # FWHM deviation - lower is better
        sel_med = np.median(sel_fwhm_dev)
        full_med = np.median(fwhm_dev)
        improvement = (full_med - sel_med) / full_med * 100 if full_med > 0 else 0
        print(f"  {'FWHM deviation':<15} {sel_med:>12.3f} {full_med:>12.3f} {improvement:>+11.0f}%")

        # Eccentricity - lower is better
        sel_med = np.median(sel_ecc)
        full_med = np.median(ecc)
        improvement = (full_med - sel_med) / full_med * 100 if full_med > 0 else 0
        print(f"  {'Eccentricity':<15} {sel_med:>12.3f} {full_med:>12.3f} {improvement:>+11.0f}%")

        # BkgContam - lower is better
        sel_med = np.median(sel_bkg)
        full_med = np.median(bkg)
        improvement = (full_med - sel_med) / full_med * 100 if full_med > 0 else 0
        print(f"  {'BkgContam':<15} {sel_med:>12.2f} {full_med:>12.2f} {improvement:>+11.0f}%")

        # Concentration - higher is better
        sel_med = np.median(sel_conc)
        full_med = np.median(conc)
        improvement = (sel_med - full_med) / abs(full_med) * 100 if full_med != 0 else 0
        print(f"  {'Concentration':<15} {sel_med:>12.2f} {full_med:>12.2f} {improvement:>+11.0f}%")

    return psf_stars


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

        # Check if we have any stamps to work with
        if len(self.stamps) == 0:
            print('Error: No valid stamps extracted. Cannot build Gaussian PSF.')
            return None, None

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

        # Check if we have any stamps to work with
        if len(self.stamps) == 0:
            print('Error: No valid stamps extracted. Cannot build Moffat PSF.')
            return None, None

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

        # Check if we have any stamps to work with
        if len(self.stamps) == 0:
            print('Error: No valid stamps extracted. Cannot build summed PSF.')
            return None

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
            Dictionary containing all PSF models and parameters, or None if no stamps available
        """
        # Extract stamps if not already done
        if self.stamps is None:
            self.extract_stamps()

        # Check if we have any stamps to work with
        if len(self.stamps) == 0:
            print('Error: No valid stamps available. Cannot build PSFs.')
            return None

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

    if results is not None and output_prefix is not None:
        builder.save_psfs(output_prefix)

    return results



def do_one(image_file, star_file, prefix='', max_stars=1000, weights=None):
    """
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
    """

    if prefix == '':
        root = image_file.split('/')[-1]
        root = root.replace('.fits', '')
        root = root.replace('.fz', '')
        prefix = root

    # Load all stars table
    all_stars = Table.read(star_file)

    # Select PSF stars using quality scoring
    psf_stars = select_psf_stars(all_stars, max_stars=max_stars, weights=weights)

    # Check if we have any PSF stars
    if len(psf_stars) == 0:
        print('Error: No valid PSF stars found. Cannot build PSF.')
        return None

    # Write selected PSF stars
    psf_stars_out = '%s_psf_stars.fits' % prefix
    psf_stars.write(psf_stars_out, format='fits', overwrite=True)
    print('Wrote PSF stars: %s' % psf_stars_out)

    # Build PSFs with consistent normalization
    results = quick_psf_build(image_file, psf_stars,
                              stamp_size=25,
                              normalize='peak',
                              output_prefix=prefix)

    # Check if PSF build succeeded
    if results is None:
        print('Error: PSF build failed. No valid stamps could be extracted.')
        return None

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

    Usage: PsfBuild [-out root] image_file star_catalog
    '''

    image_file = ''
    star_file = ''
    root = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i][:4] == '-out':
            i += 1
            root = argv[i]
        elif argv[i][0] == '-':
            print('Error: Could not interpret commands: ', argv)
            return
        elif image_file == '':
            image_file = argv[i]
        elif star_file == '':
            star_file = argv[i]
        else:
            print('Error: Too many commands: ', argv)
            return

        i += 1

    if root == '':
        root = image_file.split('/')[-1]
        root = root.replace('.fits', '')
        root = root.replace('.gz', '')

    print('       Image: ', image_file)
    print('Star catalog: ', star_file)
    print('        Root: ', root)

    do_one(image_file=image_file, star_file=star_file, prefix=root)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
