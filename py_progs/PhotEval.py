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
    Files matching ``*<FILTER>*req.smash.fits`` are read.
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

Output Columns
--------------

Source_name, RA, Dec, n_detect,
mag_mean, mag_wmean, mag_median, mag_std,       (raw, ZP = 28)
magc_mean, magc_wmean, magc_median, magc_std,   (corrected, ZP = MAGZERO)
mag_err_mean, snr_mean, chi2_nu, chi2_nu_c

Examples
--------

Evaluate r-band 30-second exposures::

    PhotEval.py -filter r -exp 30 -snr 10 -min_n 5

Evaluate Ha (N662) band, all exposure times::

    PhotEval.py -filter N662 -min_n 3 -o ha_eval.fits
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
        else:
            print(f'Unknown argument: {argv[i]}')
            return

    if outfile is None:
        exp_tag = f'_{int(exptime)}s' if exptime is not None else ''
        outfile = f'phot_eval_{filter_str}{exp_tag}.fits'

    do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile)


def do_eval(tabphot_dir, filter_str, exptime, snr_min, prob_min, min_n, outfile):
    pattern = os.path.join(tabphot_dir, f'*{filter_str}*req.smash.fits')
    files = sorted(glob(pattern))
    if not files:
        print(f'No files matching {pattern}')
        return

    print(f'Reading {len(files)} catalog files matching *{filter_str}*req.smash.fits ...')

    keep_cols = ['Source_name', 'phot_mag', 'Net', 'ErrNet', 'SNR',
                 'RA', 'Dec', 'Exptime', 'MAGZERO', 'prob', 'Filter', 'Filename']

    chunks = []
    for f in files:
        try:
            t = Table.read(f)
            # keep only the columns we need so vstack stays light
            t = t[[c for c in keep_cols if c in t.colnames]]
            chunks.append(t.to_pandas())
        except Exception as e:
            print(f'  Warning: could not read {f}: {e}')

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

    # keep only sources with enough detections
    stats = stats[stats['n_detect'] >= min_n].copy()
    print(f'Sources with >= {min_n} detections: {len(stats):,}')

    stats = stats.rename(columns={'ra_mean': 'RA', 'dec_mean': 'Dec'})

    out_cols = ['Source_name', 'RA', 'Dec', 'n_detect',
                'mag_mean', 'mag_wmean', 'mag_median', 'mag_std']
    if have_magzero:
        out_cols += ['magc_mean', 'magc_wmean', 'magc_median', 'magc_std']
    out_cols += ['mag_err_mean', 'snr_mean', 'chi2_nu']
    if have_magzero:
        out_cols += ['chi2_nu_c']

    out = Table.from_pandas(stats[out_cols])

    for col in ['mag_mean', 'mag_wmean', 'mag_std', 'mag_err_mean',
                'magc_mean', 'magc_wmean', 'magc_std']:
        if col in out.colnames:
            out[col].format = '.4f'
    for col in ['chi2_nu', 'chi2_nu_c']:
        if col in out.colnames:
            out[col].format = '.3f'

    out.write(outfile, overwrite=True)
    print(f'Wrote {outfile}')

    _print_summary(out)


def _print_summary(t):
    n = len(t)
    if n == 0:
        return
    print(f'\n--- Summary ({n:,} sources) ---')
    print(f'  mag_std  : median={np.median(t["mag_std"]):.4f}  '
          f'90th%={np.percentile(t["mag_std"], 90):.4f}  (raw, ZP=28)')
    if 'magc_std' in t.colnames:
        print(f'  magc_std : median={np.median(t["magc_std"]):.4f}  '
              f'90th%={np.percentile(t["magc_std"], 90):.4f}  (corrected)')
    print(f'  mag_err  : median={np.median(t["mag_err_mean"]):.4f}  (photon noise)')
    good = t[~np.isnan(t['chi2_nu'])]
    if len(good):
        print(f'  chi2_nu  : median={np.median(good["chi2_nu"]):.3f}  '
              f'90th%={np.percentile(good["chi2_nu"], 90):.3f}  (raw)')
    if 'chi2_nu_c' in t.colnames:
        goodc = t[~np.isnan(t['chi2_nu_c'])]
        if len(goodc):
            print(f'  chi2_nu_c: median={np.median(goodc["chi2_nu_c"]):.3f}  '
                  f'90th%={np.percentile(goodc["chi2_nu_c"], 90):.3f}  (corrected)')
    print(f'  n_detect : median={np.median(t["n_detect"]):.0f}  '
          f'max={t["n_detect"].max()}')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
