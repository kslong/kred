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
    Only process files whose Filter header value starts with FILTER
    (e.g. ``r``, ``N662``, ``N673``).  The filter name is read from
    the ``Filter`` column inside each file, not from the filename.
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

Filename, Root, Field, Filter, Exptime, Catalog, ref_col, MAGZERO,
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

If the catalog files carry a ``FIELD`` header keyword (written by MefPhot),
``Summary/{field}_mef.tab`` is updated automatically with columns
``ZP_{catalog}_{band}``, ``ZP_std_{catalog}_{band}``, ``ZP_n_{catalog}_{band}``
(e.g. ``ZP_smash_r``, ``ZP_gaia_g``).  Only rows for MEF roots processed in
this run are touched; other rows and pre-existing ZP columns are left unchanged.
"""

import sys
import os
import numpy as np
from collections import defaultdict
from glob import glob
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table
from astropy.stats import sigma_clip, mad_std


def run(filter_str=None, tabphot_dir='TabPhot', snr_min=20.0,
        ref_col_override=None, n_sigma=3.0, mag_lo=21.0, mag_hi=15.0,
        outfile='zeropoints.fits'):
    """Process all TabPhot catalogs for one filter and return the per-file results Table.

    Writes ``outfile`` (default ``zeropoints.fits``) and updates
    ``Summary/{field}_mef.tab`` as side effects.  Returns None if no files
    are found or no results are produced.

    Parameters
    ----------
    filter_str : str or None
        Filter to select (e.g. ``'r'``, ``'N662'``).  None means all filters.
    tabphot_dir : str
        Directory containing ``*.smash.fits`` / ``*.gaia.fits`` catalogs.
    snr_min : float
        Minimum SNR for a star to be included (default 20).
    ref_col_override : str or None
        Override the reference magnitude column (default: auto from catalog).
    n_sigma : float
        Sigma-clipping threshold (default 3).
    mag_lo, mag_hi : float
        Faint and bright magnitude limits for reference stars (default 21, 15).
    outfile : str or None
        Output FITS filename.  Pass None to skip writing.
    """
    pattern_smash = os.path.join(tabphot_dir, '*.smash.fits')
    pattern_gaia  = os.path.join(tabphot_dir, '*.gaia.fits')
    files = sorted(glob(pattern_smash) + glob(pattern_gaia))

    if filter_str:
        files = [f for f in files if _file_filter(f) == filter_str]

    if not files:
        print(f'CalcZeroPoint: no catalog files found in {tabphot_dir}'
              + (f' for filter {filter_str}' if filter_str else ''))
        return None

    print(f'CalcZeroPoint: processing {len(files)} catalog files'
          + (f' (filter={filter_str})' if filter_str else '') + ' ...')

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
        print('CalcZeroPoint: no results computed.')
        return None

    out = Table(rows=rows)
    for col in ['zp_calc', 'zp_wmean', 'zp_std', 'zp_err', 'zp_mad',
                'MAGZERO', 'delta_zp']:
        out[col].format = '.4f'

    if outfile:
        out.write(outfile, overwrite=True)
        print(f'\nCalcZeroPoint: wrote {outfile}  ({len(out)} rows)')
    _print_summary(out)

    by_field = defaultdict(list)
    for r in rows:
        if r['Field']:
            by_field[r['Field']].append(r)
    for field, field_rows in by_field.items():
        update_mef_tab(field_rows, field)

    return out


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

    run(filter_str=filter_str, tabphot_dir=tabphot_dir, snr_min=snr_min,
        ref_col_override=ref_col_override, n_sigma=n_sigma,
        mag_lo=mag_lo, mag_hi=mag_hi, outfile=outfile)


def _file_filter(filepath):
    """Return the short filter name from the ext 1 header (FILTER) or legacy Filter column."""
    try:
        with fits.open(filepath, memmap=True) as hdul:
            hdr = hdul[1].header
            if 'FILTER' in hdr:
                return _decode(hdr['FILTER']).split()[0]
            return _decode(hdul[1].data['Filter'][0]).split()[0]
    except Exception:
        return None


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
        with fits.open(filepath) as hdul:
            hdr = hdul[1].header
            t = Table(hdul[1].data)
    except Exception as e:
        print(f'  Warning: cannot read {filepath}: {e}')
        return None

    # file-level metadata: header takes priority, columns are fallback for old files
    def _hdr_or_col(key, col, default):
        if key in hdr:
            return hdr[key]
        if col in t.colnames:
            return np.unique(t[col])[0]
        return default

    catalog     = _decode(_hdr_or_col('CATALOG', 'Catalog', 'unknown'))
    filt        = _decode(_hdr_or_col('FILTER',  'Filter',  ''))
    magzero_hdr = float(_hdr_or_col('MAGZERO', 'MAGZERO', 28.0))
    exptime     = float(_hdr_or_col('EXPTIME',  'Exptime',  np.nan))
    field       = _decode(hdr['FIELD']) if 'FIELD' in hdr else ''

    # MEF root: always derive from the TabPhot output filename — it is named
    # after the MEF root regardless of whether MefPhot ran on MEF or PREP files.
    stem = os.path.basename(filepath)
    for sfx in ('.smash.fits', '.gaia.fits', '.fits'):
        if stem.endswith(sfx):
            stem = stem[:-len(sfx)]
            break
    root = stem

    # determine reference column
    if ref_col_override:
        ref_col = ref_col_override
    elif catalog.upper() == 'SMASH':
        ref_col = 'R'
    else:
        ref_col = 'G'

    if ref_col not in t.colnames:
        print(f'  Warning: {filepath} has no column {ref_col!r}, skipping')
        return None

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
        'Root':     root,
        'Field':    field,
        'Filter':   filt,
        'Exptime':  exptime,
        'Catalog':  catalog,
        'ref_col':  ref_col,
        'MAGZERO':  magzero_hdr,
        'zp_calc':  zp_calc,
        'zp_wmean': zp_wmean,
        'zp_std':   zp_std,
        'zp_err':   zp_err,
        'zp_mad':   zp_mad,
        'n_stars':  n_stars,
        'n_total':  n_total,
        'delta_zp': zp_calc - magzero_hdr,
    }


def update_mef_tab(rows, field):
    """Add or update ZP columns in Summary/{field}_mef.tab from CalcZeroPoint results.

    Columns are named ZP_{catalog}_{band}, ZP_std_{catalog}_{band},
    ZP_n_{catalog}_{band} (e.g. ZP_smash_r, ZP_gaia_g).  Only rows whose
    Root appears in this run are touched; all other rows are left unchanged.
    New columns are initialised to NaN / 0 before filling.
    """
    tab_path = f'Summary/{field}_mef.tab'
    if not os.path.isfile(tab_path):
        print(f'CalcZeroPoint: {tab_path} not found, skipping summary update')
        return

    mef = asc.read(tab_path)

    # group rows by (catalog, ref_col) — each combination gets its own column set
    combos = defaultdict(list)
    for r in rows:
        key = (_decode(r['Catalog']).lower(), r['ref_col'].lower())
        combos[key].append(r)

    for (catalog, ref_col), combo_rows in combos.items():
        col_zp  = f'ZP_{catalog}_{ref_col}'
        col_std = f'ZP_std_{catalog}_{ref_col}'
        col_n   = f'ZP_n_{catalog}_{ref_col}'

        if col_zp not in mef.colnames:
            mef[col_zp]  = np.full(len(mef), np.nan)
        if col_std not in mef.colnames:
            mef[col_std] = np.full(len(mef), np.nan)
        if col_n not in mef.colnames:
            mef[col_n]   = np.zeros(len(mef), dtype=int)
        mef[col_zp].info.format  = '.3f'
        mef[col_std].info.format = '.3f'

        # one TabPhot file per MEF root; use its within-file stats directly
        n_updated = 0
        for r in combo_rows:
            root = r['Root']
            if not root:
                continue
            idx = np.where(np.array(mef['Root']) == root)[0]
            if len(idx) == 0:
                print(f'  Warning: root {root!r} not found in {tab_path}')
                continue
            mef[col_zp][idx[0]]  = float(r['zp_calc'])
            mef[col_std][idx[0]] = float(r['zp_std'])
            mef[col_n][idx[0]]   = int(r['n_stars'])
            n_updated += 1

        print(f'  {col_zp}: updated {n_updated} of {len(mef)} rows')

    asc.write(mef, tab_path, format='fixed_width_two_line', overwrite=True)
    print(f'Wrote updated {tab_path}')


def _print_summary(t):
    print(f'\n--- Summary ---')
    for filt in np.unique(t['Filter']):
        tf = t[t['Filter'] == filt]
        short = _decode(filt).split()[0]
        for exptime in np.unique(tf['Exptime']):
            te = tf[tf['Exptime'] == exptime]
            print(f'  Filter {short!r}  exptime={exptime:.0f}s  ({len(te)} files):')
            print(f'    zp_calc  : mean={np.mean(te["zp_calc"]):.4f}  '
                  f'std={np.std(te["zp_calc"]):.4f}  '
                  f'range=[{np.min(te["zp_calc"]):.4f}, {np.max(te["zp_calc"]):.4f}]')
            print(f'    delta_zp : mean={np.mean(te["delta_zp"]):.4f}  '
                  f'std={np.std(te["delta_zp"]):.4f}')
            print(f'    zp_std   : median={np.median(te["zp_std"]):.4f}  (per-file scatter)')
            print(f'    n_stars  : median={np.median(te["n_stars"]):.0f}')


if __name__ == '__main__':
    steer(sys.argv)
