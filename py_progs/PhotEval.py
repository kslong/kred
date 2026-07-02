#!/usr/bin/env python
"""PhotEval - Photometric consistency evaluation across overlapping frames

Space Telescope Science Institute

Command Line Usage
------------------

::

    PhotEval.py [-h] [-dir DIR] [-filter FILTER] [-exp EXPTIME]
                [-snr MIN_SNR] [-prob MIN_PROB] [-min_n MIN_N] [-o OUTPUT]

Reads all ``*<FILTER>*.smash.fits`` catalogs from ``DIR`` (default:
``TabPhot/``), filters on exposure time and quality, groups detections by
``Source_name``, and writes a summary table with per-source photometric
statistics.

Two sets of magnitude statistics are produced:

* **Raw** (``mag_*``): based on ``phot_mag``, which uses a fixed instrumental
  zero point of 28.
* **Corrected** (``magc_*``): ``phot_mag + MAGZERO - 28``, i.e. the true
  calibrated magnitude using the per-image zero point measured by MagZero.
  When images have different zero points these corrections matter; the
  corrected scatter (``magc_std``) is the meaningful repeatability metric.

The key diagnostic is the scatter (``magc_std``) relative to the expected
photon-noise floor (``mag_err_mean``) and the reduced chi^2 (``chi2_nu_c``).

Optional Arguments
------------------

-h
    Print this help and exit.

-dir DIR
    Directory containing the ``*.smash.fits`` catalogs (default: ``TabPhot``).

-filter FILTER
    Filter string used to select files, e.g. ``r``, ``N662``, ``N673``.
    Files matching ``*<FILTER>*req.smash.fits`` are preferred; if none
    exist for that filter, ``*<FILTER>*req.gaia.fits`` are used instead.
    Default: ``r``.

-exp EXPTIME
    Exposure time to select, in seconds (e.g. ``30``).  If omitted all
    exposure times are used.

-snr MIN_SNR
    Minimum SNR for a detection to be included (default: 10).

-prob MIN_PROB
    Minimum SMASH star probability for inclusion (default: 0.5).

-min_n MIN_N
    Minimum number of qualifying detections required to include a source
    in the output (default: 3).

-o OUTPUT
    Output FITS filename (default: ``phot_eval_<FILTER>_<EXPTIME>.fits``
    in the current directory).

-zp_table FILE
    FITS table of empirical zero points (e.g. ``zeropoints.fits`` produced
    by CalcZeroPoint).  Must contain ``Filename`` and ``zp_calc`` columns.
    When supplied, a second set of corrected magnitudes is computed using
    ``zp_calc`` in place of the header ``MAGZERO``, giving additional output
    columns ``magc_emp_*`` and ``chi2_nu_emp`` for direct comparison.

Output Columns
--------------

Source_name, RA, Dec, n_detect,
mag_mean, mag_wmean, mag_median, mag_std,           (raw, ZP = 28)
magc_mean, magc_wmean, magc_median, magc_std,       (corrected, ZP = MAGZERO)
magc_emp_mean, magc_emp_wmean, magc_emp_median, magc_emp_std,  (empirical ZP, if -zp_table given)
mag_err_mean, snr_mean, chi2_nu, chi2_nu_c, chi2_nu_emp

Examples
--------

Evaluate r-band 30-second exposures::

    PhotEval.py -filter r -exp 30 -snr 10 -min_n 5

Evaluate Ha (N662) band, all exposure times::

    PhotEval.py -filter N662 -min_n 3 -o ha_eval.fits

Compare header MAGZERO vs empirical ZP from CalcZeroPoint::

    PhotEval.py -filter r -zp_table zeropoints.fits -o phot_eval_r_emp.fits
"""

import sys
import os
import numpy as np
import pandas as pd
from glob import glob
from astropy.table import Table
from astropy.io import fits


def steer(argv):
    tabphot_dir = 'TabPhot'
    filter_str = 'r'
    exptime = None
    snr_min = 10.0
    prob_min = 0.5
    min_n = 3
    outfile = None
    zp_table_file = None

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-dir':
            tabphot_dir = argv[i + 1]
            i += 2
        elif argv[i] == '-filter':
            filter_str = argv[i + 1]
            i += 2
        elif argv[i] == '-exp':
            exptime = float(argv[i + 1])
            i += 2
        elif argv[i] == '-snr':
            snr_min = float(argv[i + 1])
            i += 2
        elif argv[i] == '-prob':
            prob_min = float(argv[i + 1])
            i += 2
        elif argv[i] == '-min_n':
            min_n = int(argv[i + 1])
            i += 2
        elif argv[i] == '-o':
            outfile = argv[i + 1]
            i += 2
        elif argv[i] == '-zp_table':
            zp_table_file = argv[i + 1]
            i += 2
        else:
            print(f'Unknown argument: {argv[i]}')
            return

    if outfile is None:
        exp_tag = f'_{int(exptime)}s' if exptime is not None else ''
        outfile = f'phot_eval_{filter_str}{exp_tag}.fits'

    zp_lookup = {}
    if zp_table_file is not None:
        try:
            from astropy.table import Table as _Table
            zpt = _Table.read(zp_table_file)
            for row in zpt:
                fname = str(row['Filename']).strip()
                zp_lookup[fname] = float(row['zp_calc'])
            print(f'Loaded empirical ZPs for {len(zp_lookup)} files from {zp_table_file}')
        except Exception as e:
            print(f'Warning: could not load -zp_table {zp_table_file}: {e}')

    do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile, zp_lookup)


def do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile, zp_lookup=None):
    pattern_smash = os.path.join(tabphot_dir, f'*{filter_str}*req.smash.fits')
    pattern_gaia  = os.path.join(tabphot_dir, f'*{filter_str}*req.gaia.fits')
    files = sorted(glob(pattern_smash))
    if files:
        pattern_used = pattern_smash
        catalog_used = 'smash'
    else:
        files = sorted(glob(pattern_gaia))
        pattern_used = pattern_gaia
        catalog_used = 'gaia'
    if not files:
        print(f'No files matching {pattern_smash} or {pattern_gaia}')
        return

    print(f'Reading {len(files)} catalog files matching *{filter_str}*req.{catalog_used}.fits ...')

    per_row_cols = ['Source_name', 'phot_mag', 'Net', 'ErrNet', 'SNR', 'RA', 'Dec', 'prob']

    have_emp = zp_lookup is not None and len(zp_lookup) > 0
    n_emp_matched = 0

    chunks = []
    for f in files:
        try:
            with fits.open(f) as hdul:
                hdr = hdul[1].header
                t = Table(hdul[1].data)
            t = t[[c for c in per_row_cols if c in t.colnames]]
            # file-level metadata: header first, column fallback for old files
            t['Exptime'] = float(hdr['EXPTIME']) if 'EXPTIME' in hdr else (
                float(np.unique(t['Exptime'])[0]) if 'Exptime' in t.colnames else np.nan)
            t['MAGZERO'] = float(hdr['MAGZERO']) if 'MAGZERO' in hdr else (
                float(np.unique(t['MAGZERO'])[0]) if 'MAGZERO' in t.colnames else 28.0)
            t['Filter']  = str(hdr['FILTER']) if 'FILTER' in hdr else (
                str(t['Filter'][0]) if 'Filter' in t.colnames else '')
            if have_emp:
                basename = os.path.basename(f)
                zp_emp = zp_lookup.get(basename, np.nan)
                if np.isfinite(zp_emp):
                    n_emp_matched += 1
                t['ZP_EMP'] = float(zp_emp)
            chunks.append(t.to_pandas())
        except Exception as e:
            print(f'  Warning: could not read {f}: {e}')

    if have_emp:
        print(f'Empirical ZP matched for {n_emp_matched} of {len(files)} files')

    if not chunks:
        print('No data loaded.')
        return

    data = pd.concat(chunks, ignore_index=True)

    # decode bytes columns that pandas loads as objects
    if data['Source_name'].dtype == object:
        data['Source_name'] = data['Source_name'].apply(
            lambda x: x.decode() if isinstance(x, bytes) else x
        )

    print(f'Total rows loaded: {len(data):,}')

    # quality filters
    if exptime is not None:
        data = data[np.isclose(data['Exptime'], exptime, atol=0.5)]
        print(f'After Exptime={exptime}s cut: {len(data):,}')

    data = data[data['SNR'] >= snr_min]
    data = data[data['Net'] > 0]
    if 'prob' in data.columns:
        # SMASH prob == 99.99 is a sentinel for "no value"; exclude those
        data = data[(data['prob'] >= prob_min) & (data['prob'] < 99.0)]

    print(f'After quality cuts (SNR>={snr_min}, prob>={prob_min}): {len(data):,}')

    # per-detection photon-noise magnitude error: sigma_m = 2.5/ln(10) * ErrNet/Net
    data['mag_err'] = (2.5 / np.log(10)) * data['ErrNet'] / data['Net']
    data['ivar'] = 1.0 / data['mag_err'].clip(lower=1e-6) ** 2

    # zero-point-corrected magnitude: phot_mag uses ZP=28; true mag shifts by (MAGZERO-28)
    if 'MAGZERO' in data.columns:
        data['mag_corr'] = data['phot_mag'] + data['MAGZERO'] - 28.0
        have_magzero = True
    else:
        have_magzero = False

    # empirical ZP correction (only if -zp_table supplied and matched)
    have_emp = 'ZP_EMP' in data.columns and data['ZP_EMP'].notna().any()
    if have_emp:
        data['mag_emp'] = data['phot_mag'] + data['ZP_EMP'] - 28.0

    # --- vectorized groupby aggregation ---
    grp = data.groupby('Source_name', sort=False)

    agg_dict = dict(
        n_detect=('phot_mag', 'count'),
        mag_mean=('phot_mag', 'mean'),
        mag_median=('phot_mag', 'median'),
        mag_std=('phot_mag', 'std'),
        mag_err_mean=('mag_err', 'mean'),
        snr_mean=('SNR', 'mean'),
        ra_mean=('RA', 'mean'),
        dec_mean=('Dec', 'mean'),
    )
    if have_magzero:
        agg_dict.update(
            magc_mean=('mag_corr', 'mean'),
            magc_median=('mag_corr', 'median'),
            magc_std=('mag_corr', 'std'),
        )
    if have_emp:
        agg_dict.update(
            magc_emp_mean=('mag_emp', 'mean'),
            magc_emp_median=('mag_emp', 'median'),
            magc_emp_std=('mag_emp', 'std'),
        )
    agg = grp.agg(**agg_dict).reset_index()

    # weighted means and reduced chi^2 for both raw and corrected magnitudes
    def _weighted_stats(df, mag_col, ivar_col):
        grp2 = df.groupby('Source_name', sort=False)
        w = grp2.agg(
            w_sum=(ivar_col, 'sum'),
            w_mag_sum=(mag_col, lambda x: (x * df.loc[x.index, ivar_col]).sum()),
        ).reset_index()
        w['wmean'] = w['w_mag_sum'] / w['w_sum']
        df2 = df.merge(w[['Source_name', 'wmean']], on='Source_name', how='left')
        df2['chi2c'] = (df2[mag_col] - df2['wmean']) ** 2 * df2[ivar_col]
        chi2 = df2.groupby('Source_name', sort=False).agg(
            chi2_sum=('chi2c', 'sum'),
        ).reset_index()
        return w[['Source_name', 'wmean']], chi2

    wmean_raw, chi2_raw = _weighted_stats(data, 'phot_mag', 'ivar')
    wmean_raw = wmean_raw.rename(columns={'wmean': 'mag_wmean'})
    chi2_raw = chi2_raw.rename(columns={'chi2_sum': 'chi2_sum_raw'})

    stats = agg.merge(wmean_raw, on='Source_name')
    stats = stats.merge(chi2_raw, on='Source_name')
    stats['chi2_nu'] = np.where(
        stats['n_detect'] > 1,
        stats['chi2_sum_raw'] / (stats['n_detect'] - 1),
        np.nan,
    )

    if have_magzero:
        wmean_corr, chi2_corr = _weighted_stats(data, 'mag_corr', 'ivar')
        wmean_corr = wmean_corr.rename(columns={'wmean': 'magc_wmean'})
        chi2_corr = chi2_corr.rename(columns={'chi2_sum': 'chi2_sum_corr'})
        stats = stats.merge(wmean_corr, on='Source_name')
        stats = stats.merge(chi2_corr, on='Source_name')
        stats['chi2_nu_c'] = np.where(
            stats['n_detect'] > 1,
            stats['chi2_sum_corr'] / (stats['n_detect'] - 1),
            np.nan,
        )

    if have_emp:
        wmean_emp, chi2_emp = _weighted_stats(data, 'mag_emp', 'ivar')
        wmean_emp = wmean_emp.rename(columns={'wmean': 'magc_emp_wmean'})
        chi2_emp = chi2_emp.rename(columns={'chi2_sum': 'chi2_sum_emp'})
        stats = stats.merge(wmean_emp, on='Source_name')
        stats = stats.merge(chi2_emp, on='Source_name')
        stats['chi2_nu_emp'] = np.where(
            stats['n_detect'] > 1,
            stats['chi2_sum_emp'] / (stats['n_detect'] - 1),
            np.nan,
        )

    # keep only sources with enough detections
    stats = stats[stats['n_detect'] >= min_n].copy()
    print(f'Sources with >= {min_n} detections: {len(stats):,}')

    stats = stats.rename(columns={'ra_mean': 'RA', 'dec_mean': 'Dec'})

    out_cols = ['Source_name', 'RA', 'Dec', 'n_detect',
                'mag_mean', 'mag_wmean', 'mag_median', 'mag_std']
    if have_magzero:
        out_cols += ['magc_mean', 'magc_wmean', 'magc_median', 'magc_std']
    if have_emp:
        out_cols += ['magc_emp_mean', 'magc_emp_wmean', 'magc_emp_median', 'magc_emp_std']
    out_cols += ['mag_err_mean', 'snr_mean', 'chi2_nu']
    if have_magzero:
        out_cols += ['chi2_nu_c']
    if have_emp:
        out_cols += ['chi2_nu_emp']

    out = Table.from_pandas(stats[out_cols])

    for col in ['mag_mean', 'mag_wmean', 'mag_std', 'mag_err_mean',
                'magc_mean', 'magc_wmean', 'magc_std',
                'magc_emp_mean', 'magc_emp_wmean', 'magc_emp_std']:
        if col in out.colnames:
            out[col].format = '.4f'
    for col in ['chi2_nu', 'chi2_nu_c', 'chi2_nu_emp']:
        if col in out.colnames:
            out[col].format = '.3f'

    out.write(outfile, overwrite=True)
    print(f'Wrote {outfile}')

    _print_summary(out)
    return out


def _print_summary(t):
    n = len(t)
    if n == 0:
        return
    print(f'\n--- Summary ({n:,} sources) ---')
    print(f'  mag_std  : median={np.median(t["mag_std"]):.4f}  '
          f'90th%={np.percentile(t["mag_std"], 90):.4f}  (raw, ZP=28)')
    if 'magc_std' in t.colnames:
        print(f'  magc_std    : median={np.median(t["magc_std"]):.4f}  '
              f'90th%={np.percentile(t["magc_std"], 90):.4f}  (header MAGZERO)')
    if 'magc_emp_std' in t.colnames:
        print(f'  magc_emp_std: median={np.median(t["magc_emp_std"]):.4f}  '
              f'90th%={np.percentile(t["magc_emp_std"], 90):.4f}  (empirical ZP)')
    print(f'  mag_err     : median={np.median(t["mag_err_mean"]):.4f}  (photon noise)')
    good = t[~np.isnan(t['chi2_nu'])]
    if len(good):
        print(f'  chi2_nu     : median={np.median(good["chi2_nu"]):.3f}  '
              f'90th%={np.percentile(good["chi2_nu"], 90):.3f}  (raw)')
    if 'chi2_nu_c' in t.colnames:
        goodc = t[~np.isnan(t['chi2_nu_c'])]
        if len(goodc):
            print(f'  chi2_nu_c   : median={np.median(goodc["chi2_nu_c"]):.3f}  '
                  f'90th%={np.percentile(goodc["chi2_nu_c"], 90):.3f}  (header MAGZERO)')
    if 'chi2_nu_emp' in t.colnames:
        goode = t[~np.isnan(t['chi2_nu_emp'])]
        if len(goode):
            print(f'  chi2_nu_emp : median={np.median(goode["chi2_nu_emp"]):.3f}  '
                  f'90th%={np.percentile(goode["chi2_nu_emp"], 90):.3f}  (empirical ZP)')
    print(f'  n_detect : median={np.median(t["n_detect"]):.0f}  '
          f'max={t["n_detect"].max()}')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
