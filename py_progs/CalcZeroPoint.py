#!/usr/bin/env python
"""CalcZeroPoint - Calculate photometric zero points from TabPhot catalogs

Space Telescope Science Institute

Command Line Usage
------------------

::

    CalcZeroPoint.py [-h] [-dir DIR] [-filter FILTER] [-snr MIN_SNR]
                     [-ref COL] [-sigma N] [-mag_lo LO] [-mag_hi HI]
                     [-o OUTPUT]

For each ``*.smash.fits`` or ``*.gaia.fits`` catalog in ``DIR`` (default:
``TabPhot``), computes a best-fit photometric zero point by comparing
instrumental magnitudes (``phot_mag``, which use a fixed ZP = 28) to
reference catalog magnitudes (``R`` for SMASH, ``G`` for Gaia).

The zero point for a single star is::

    zp_i = 28 + (ref_mag_i - phot_mag_i)

After quality cuts, the distribution of ``zp_i`` is sigma-clipped and
summarised.  The output table has one row per input file.

Optional Arguments
------------------

-h
    Print this help and exit.

-dir DIR
    Directory containing the ``*.smash.fits`` / ``*.gaia.fits`` catalogs
    (default: ``TabPhot``).

-filter FILTER
    Only process files whose name contains FILTER (e.g. ``r``, ``N662``).
    Default: process all catalog files found.

-snr MIN_SNR
    Minimum SNR for a detection to be used (default: 20).

-ref COL
    Reference magnitude column to use (overrides the per-file default of
    ``R`` for SMASH and ``G`` for Gaia).

-sigma N
    Number of sigma for iterative sigma-clipping (default: 3).

-mag_lo LO
    Faint magnitude limit for reference stars (default: 21).
    Excludes faint stars where photon noise dominates.

-mag_hi HI
    Bright magnitude limit for reference stars (default: 15).
    Excludes stars that may be saturated.

-o OUTPUT
    Output FITS filename (default: ``zeropoints.fits``).

Output Columns
--------------

Filename, Filter, Exptime, Catalog, ref_col, MAGZERO,
zp_calc, zp_wmean, zp_std, zp_err, zp_mad,
n_stars, n_total, delta_zp

where:

* ``zp_calc``  -- sigma-clipped median of individual ZP estimates
* ``zp_wmean`` -- inverse-variance-weighted mean ZP
* ``zp_std``   -- std of clipped ZP estimates (scatter / repeatability)
* ``zp_err``   -- standard error on zp_calc (zp_std / sqrt(n_stars))
* ``zp_mad``   -- median absolute deviation of clipped estimates
* ``n_stars``  -- number of stars used after sigma-clipping
* ``n_total``  -- number of stars passing initial quality cuts
* ``delta_zp`` -- zp_calc - MAGZERO (should be ~0 if header ZP is good)

Examples
--------

Process all SMASH r-band files::

    CalcZeroPoint.py -filter r -snr 20

Process all filters, bright stars only::

    CalcZeroPoint.py -mag_hi 13 -mag_lo 18 -o zp_bright.fits
"""

import sys
import os
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.stats import sigma_clip, mad_std


def steer(argv):
    tabphot_dir = 'TabPhot'
    filter_str = None
    snr_min = 20.0
    ref_col_override = None
    n_sigma = 3.0
    mag_lo = 21.0
    mag_hi = 15.0
    outfile = 'zeropoints.fits'

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
        elif argv[i] == '-snr':
            snr_min = float(argv[i + 1])
            i += 2
        elif argv[i] == '-ref':
            ref_col_override = argv[i + 1]
            i += 2
        elif argv[i] == '-sigma':
            n_sigma = float(argv[i + 1])
            i += 2
        elif argv[i] == '-mag_lo':
            mag_lo = float(argv[i + 1])
            i += 2
        elif argv[i] == '-mag_hi':
            mag_hi = float(argv[i + 1])
            i += 2
        elif argv[i] == '-o':
            outfile = argv[i + 1]
            i += 2
        else:
            print(f'Unknown argument: {argv[i]}')
            return

    pattern_smash = os.path.join(tabphot_dir, '*.smash.fits')
    pattern_gaia = os.path.join(tabphot_dir, '*.gaia.fits')
    files = sorted(glob(pattern_smash) + glob(pattern_gaia))

    if filter_str:
        files = [f for f in files if filter_str in os.path.basename(f)]

    if not files:
        print(f'No catalog files found in {tabphot_dir}')
        return

    print(f'Processing {len(files)} catalog files ...')

    rows = []
    for f in files:
        row = do_one(f, snr_min, ref_col_override, n_sigma, mag_hi, mag_lo)
        if row is not None:
            rows.append(row)
            print(f'  {os.path.basename(f):55s}  '
                  f'zp_calc={row["zp_calc"]:7.4f}  '
                  f'MAGZERO={row["MAGZERO"]:7.4f}  '
                  f'delta={row["delta_zp"]:+.4f}  '
                  f'std={row["zp_std"]:.4f}  '
                  f'n={row["n_stars"]}')

    if not rows:
        print('No results computed.')
        return

    out = Table(rows=rows)

    for col in ['zp_calc', 'zp_wmean', 'zp_std', 'zp_err', 'zp_mad',
                'MAGZERO', 'delta_zp']:
        out[col].format = '.4f'

    out.write(outfile, overwrite=True)
    print(f'\nWrote {outfile}  ({len(out)} rows)')
    _print_summary(out)


def _decode(val):
    """Return a plain str from bytes, numpy bytes, numpy str, or str."""
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode().strip()
    # numpy string scalars print with np.str_('...') wrapper via str(); use item() to unwrap
    if hasattr(val, 'item'):
        val = val.item()
    return str(val).strip()


def do_one(filepath, snr_min, ref_col_override, n_sigma, mag_hi, mag_lo):
    try:
        t = Table.read(filepath)
    except Exception as e:
        print(f'  Warning: cannot read {filepath}: {e}')
        return None

    # determine catalog type and reference column
    catalog = _decode(t['Catalog'][0]) if 'Catalog' in t.colnames else 'unknown'
    if ref_col_override:
        ref_col = ref_col_override
    elif catalog.upper() == 'SMASH':
        ref_col = 'R'
    else:
        ref_col = 'G'

    if ref_col not in t.colnames:
        print(f'  Warning: {filepath} has no column {ref_col!r}, skipping')
        return None

    # header-level metadata (same for all rows in one file)
    magzero_hdr = float(np.unique(t['MAGZERO'])[0]) if 'MAGZERO' in t.colnames else 28.0
    filt = _decode(np.unique(t['Filter'])[0]) if 'Filter' in t.colnames else ''
    exptime = float(np.unique(t['Exptime'])[0]) if 'Exptime' in t.colnames else np.nan

    # quality cuts
    mask = (t['SNR'] >= snr_min) & (t['Net'] > 0)
    if 'prob' in t.colnames:
        prob = np.array(t['prob'])
        mask &= (prob >= 0.5) & (prob < 99.0)

    ref_mag = np.ma.filled(t[ref_col], fill_value=np.nan).astype(float)
    mask &= np.isfinite(ref_mag)
    mask &= (ref_mag >= mag_hi) & (ref_mag <= mag_lo)

    t_good = t[mask]
    ref_good = ref_mag[mask]
    n_total = int(np.sum(mask))

    if n_total < 5:
        print(f'  Warning: {filepath} has only {n_total} stars after cuts, skipping')
        return None

    # per-star ZP estimates
    zp_vals = 28.0 + (ref_good - np.array(t_good['phot_mag'], dtype=float))

    # sigma-clip
    clipped = sigma_clip(zp_vals, sigma=n_sigma, maxiters=5, masked=True)
    zp_use = zp_vals[~clipped.mask]
    n_stars = len(zp_use)

    if n_stars < 3:
        print(f'  Warning: {filepath} has only {n_stars} stars after clipping, skipping')
        return None

    zp_calc = float(np.median(zp_use))
    zp_std = float(np.std(zp_use, ddof=1))
    zp_err = zp_std / np.sqrt(n_stars)
    zp_mad = float(mad_std(zp_use))

    # inverse-variance-weighted mean using photon-noise mag errors
    net_use = np.array(t_good['Net'], dtype=float)[~clipped.mask]
    err_use = np.array(t_good['ErrNet'], dtype=float)[~clipped.mask]
    mag_err = (2.5 / np.log(10)) * np.abs(err_use) / np.abs(net_use)
    mag_err = np.clip(mag_err, 1e-6, None)
    ivar = 1.0 / mag_err ** 2
    zp_wmean = float(np.sum(ivar * zp_use) / np.sum(ivar))

    return {
        'Filename': os.path.basename(filepath),
        'Filter': filt,
        'Exptime': exptime,
        'Catalog': catalog,
        'ref_col': ref_col,
        'MAGZERO': magzero_hdr,
        'zp_calc': zp_calc,
        'zp_wmean': zp_wmean,
        'zp_std': zp_std,
        'zp_err': zp_err,
        'zp_mad': zp_mad,
        'n_stars': n_stars,
        'n_total': n_total,
        'delta_zp': zp_calc - magzero_hdr,
    }


def _print_summary(t):
    print(f'\n--- Summary ---')
    for filt in np.unique(t['Filter']):
        tf = t[t['Filter'] == filt]
        print(f'  Filter {_decode(filt)!r}  ({len(tf)} files):')
        print(f'    zp_calc  : mean={np.mean(tf["zp_calc"]):.4f}  '
              f'std={np.std(tf["zp_calc"]):.4f}  '
              f'range=[{np.min(tf["zp_calc"]):.4f}, {np.max(tf["zp_calc"]):.4f}]')
        print(f'    delta_zp : mean={np.mean(tf["delta_zp"]):.4f}  '
              f'std={np.std(tf["delta_zp"]):.4f}')
        print(f'    zp_std   : median={np.median(tf["zp_std"]):.4f}  (per-file scatter)')
        print(f'    n_stars  : median={np.median(tf["n_stars"]):.0f}')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
