#!/usr/bin/env python
# coding: utf-8

"""CheckPhot - Photometric consistency check for a DECam field

Space Telescope Science Institute

Synopsis
--------

Runs the complete photometric consistency workflow for one field and
produces two concise summary tables: one for intra-filter consistency
(are all exposures of the same filter on the same flux scale?) and one
for inter-filter consistency (does the same star have the same DN in the
r-band image as in the emission-line images?).

Command Line Usage
------------------

::

    CheckPhot.py LMC_c42
    CheckPhot.py -np 16 LMC_c42
    CheckPhot.py -no_mefphot LMC_c42        # skip MefPhot if TabPhot/ already populated
    CheckPhot.py -redo LMC_c42              # force reprocessing of all TabPhot files
    CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits   # reprint summary only
    CheckPhot.py -config myconfig.txt LMC_c42   # alternate image-group config

Flags
-----

``-np N``
    Number of parallel MefPhot workers (default 8).

``-no_mefphot``
    Skip the MefPhot step; use existing ``TabPhot/`` files.

``-redo``
    Reprocess already-existing CalcZeroPoint/PhotEval outputs.

``-o DIR``
    Write all outputs to DIR (default: ``CheckPhot``).

``-config FILE``
    Alternate image-group config file (default: ``DeMCELS_images.txt``
    in current directory, then ``$KRED/config/``).

``-summary FITS``
    Read an existing ``*_phot_check.fits`` and reprint the formatted
    summary without running any photometry.

``-ref COL``
    Additional ZP reference column for CalcZeroPoint (may be repeated).

``-h``
    Print this help.

Description
-----------

Exposures are grouped by *Image label* read from ``config/DeMCELS_images.txt``
(e.g. ``N662`` = 400-800 s N662 images, ``N662_s`` = 30-60 s N662 images).
Files in ``TabPhot/`` are assigned to groups by matching the ``FILTER``
header keyword and ``EXPTIME`` (within 5 s).  Four steps then run in sequence:

1. **MefPhot** — forced aperture photometry at SMASH and Gaia positions.
   Already-processed files are skipped unless ``-redo`` is given.

2. **CalcZeroPoint** — per-file zero-point fit vs catalog for each Image group.
   Measures how far each file's header MAGZERO deviates from the measured ZP.

3. **PhotEval** — groups detections of the same star across overlapping frames
   and computes ``chi2_nu_c``: the reduced chi-squared of the MAGZERO-corrected
   magnitudes relative to photon noise.

4. **PhotAnal** — cross-matches the same star between two Image groups and
   measures the scatter in ``(magc_1 - magc_2)``, which sets the minimum
   stellar residual in CleanStars output.  Pairs are formed automatically
   within each tier (long vs short exposures).

Output
------

All outputs are written to ``CheckPhot/`` (or the directory given via ``-o``).
All filenames include the field name so multiple fields share the directory
without collision.

``{field}_phot_check.fits`` — single FITS file with two binary table extensions:

**INTRA** (ext 1) — one row per Image group.  Zero-point comparison columns
are named ``{stat}_{ref}`` where ``{ref}`` encodes the catalog and reference
band (e.g. ``smash_r``, ``gaia_g``).  Consistency columns: ``N_sources``,
``chi2_nu_c_med``, ``chi2_nu_c_emp_med``, ``sigma_c_magzero``,
``sigma_c_gaia_g``, ``sigma_c_smash_r``.

**INTER** (ext 2) — one row per Image-group pair: ``Image1``, ``Image2``,
``Filter1``, ``Filter2``, ``N_stars``, ``mean_delta_mag``, ``sigma_c_magzero``,
``sigma_c_gaia_g``, ``sigma_c_smash_r``, ``sigma_residual``, ``residual_pct``,
``color_term_b``, ``N_pairs``, ``std_dzp_magzero``, ``std_dzp_gaia_g``,
``std_dzp_smash_r``.

``sigma_c_magzero`` and ``sigma_c_gaia_g`` / ``sigma_c_smash_r`` name the
same concept in both extensions — per-star scatter after applying a given
ZP reference — so the INTRA and INTER tables can be compared directly.

Intermediate files in ``CheckPhot/``:

- ``{field}_zeropoints_{label}.fits``         — CalcZeroPoint per-file ZP table
- ``{field}_phot_eval_{label}.fits``          — PhotEval per-source chi2 table
- ``{field}_filter_compare_{l1}_{l2}.fits``   — PhotAnal matched-star table

Interpreting the tables
-----------------------

Intra-filter:
  ``chi2_nu_c_emp_med ≈ 1`` means frames are mutually consistent at the
  photon-noise level.  ``std_dzp_{ref}`` is the MAD scatter of per-file
  MAGZERO deviations (robust to outliers); > 0.06 mag suggests the
  empirical ZP correction may help.  ``sigma_c_magzero`` vs
  ``sigma_c_smash_r`` shows whether MAGZERO or the SMASH calibration
  gives better intra-group consistency.

Inter-filter:
  ``residual_pct`` is the minimum percentage of stellar flux that will
  remain as a residual after CleanStars continuum subtraction.
  ``sigma_residual`` is the scatter after removing the color-term trend;
  this is the irreducible floor set by stellar color diversity.
  ``std_dzp_magzero`` measures image-to-image stability of the
  inter-filter flux ratio (the inter-filter analogue of ``std_dzp`` in
  the INTRA table).

History
-------

260701 ksl  Written
260702 ksl  Split ZP stats by (catalog, ref_col); add -ref flag
260703 ksl  Image-group approach from DeMCELS_images.txt; -summary flag; column renames
"""

import sys
import os
import numpy as np
from glob import glob
from collections import defaultdict
from astropy.table import Table
from astropy.stats import mad_std

import MefPhot
import CalcZeroPoint
import PhotEval
import PhotAnal


_CONFIG_NAME = 'DeMCELS_images.txt'

# Approximate central wavelength [Å] for display / sort ordering
_FILTER_WAVE = {
    'N501': 5010, 'N540': 5400,
    'r':    6400,
    'N662': 6620, 'N673': 6730, 'N708': 7080,
}


def _fwave(filt):
    """Return sort key for a filter name (strip _s suffix first)."""
    return _FILTER_WAVE.get(filt.rstrip('_s').rstrip('s') if filt.endswith('_s')
                            else filt, 9999)


def _best_label(vals, labels, tol=0.003):
    """Return the label(s) for the minimum value; '≈equal' if all within tol."""
    finite = [(v, l) for v, l in zip(vals, labels) if np.isfinite(v)]
    if not finite:
        return '?'
    mn = min(v for v, _ in finite)
    winners = [l for v, l in finite if v - mn < tol]
    if len(winners) == len(finite):
        return '≈equal'
    return '/'.join(winners)


def _print_summary(intra, inter, field):
    """Print formatted INTRA and INTER comparison tables."""
    bar = '=' * 80

    # ------------------------------------------------------------------ INTRA
    print(f'\n{bar}')
    print(f'  INTRA-FILTER SUMMARY  —  {field}')
    print(f'  sigma_c: per-star rms after ZP correction [mag]  '
          f'(chi2_emp ~ 1 → photon-noise limited)')
    print(bar)
    print(f'  {"Image":<10} {"Filt":<5} {"N_src":>8}  {"chi2_emp":>8}  '
          f'{"MagZero":>8} {"GaiaG":>8} {"SmashR":>8}  best')
    print(f'  {"-"*10} {"-"*5} {"-"*8}  {"-"*8}  {"-"*8} {"-"*8} {"-"*8}  -------')

    prev_tier = None
    for row in intra:
        img  = str(row['Image'])
        tier = 'short' if img.endswith('_s') else 'long'
        if prev_tier is not None and tier != prev_tier:
            print()
        prev_tier = tier
        mz   = float(row['sigma_c_magzero'])  if 'sigma_c_magzero'  in intra.colnames else np.nan
        gg   = float(row['sigma_c_gaia_g'])   if 'sigma_c_gaia_g'   in intra.colnames else np.nan
        sm   = float(row['sigma_c_smash_r'])  if 'sigma_c_smash_r'  in intra.colnames else np.nan
        chi  = float(row['chi2_nu_c_emp_med'])if 'chi2_nu_c_emp_med'in intra.colnames else np.nan
        ns   = int(row['N_sources'])          if 'N_sources'         in intra.colnames else 0
        best = _best_label([mz, gg, sm], ['MagZero', 'GaiaG', 'SmashR'])
        print(f'  {img:<10} {str(row["Filter"]):<5} {ns:>8}  '
              f'{chi:>8.3f}  {mz:>8.4f} {gg:>8.4f} {sm:>8.4f}  {best}')

    # ------------------------------------------------------------------ INTER
    print(f'\n{bar}')
    print(f'  INTER-FILTER SUMMARY  —  {field}')
    print(f'  sigma_c: scatter of (magc_F1 − magc_F2) per star [mag]')
    print(f'  resid_MZ: scatter after color-term fit (MAGZERO) = continuum subtraction floor')
    print(bar)
    hdr_a = (f'  {"Pair":<22} {"N_stars":>8}  '
             f'{"---- sigma_c ----":^28}  {"resid_MZ":>8}  {"resid%":>6}  '
             f'{"---- std_dzp ----":^28}')
    hdr_b = (f'  {"":22} {"":>8}  '
             f'{"MagZero":>8} {"GaiaG":>8} {"SmashR":>8}  {"":>8}  {"":>6}  '
             f'{"MagZero":>8} {"GaiaG":>8} {"SmashR":>8}')
    div   = (f'  {"-"*22} {"-"*8}  '
             f'{"-"*8} {"-"*8} {"-"*8}  {"-"*8}  {"-"*6}  '
             f'{"-"*8} {"-"*8} {"-"*8}')
    print(hdr_a)
    print(hdr_b)
    print(div)

    def _pair_sort_key(row):
        f1, f2 = str(row['Filter1']), str(row['Filter2'])
        i1, i2 = str(row['Image1']),  str(row['Image2'])
        tier   = 1 if (i1.endswith('_s') or i2.endswith('_s')) else 0
        is_bb  = 0 if (f1 == 'r' or f2 == 'r') else 1
        return (tier, is_bb, min(_fwave(f1), _fwave(f2)), max(_fwave(f1), _fwave(f2)))

    prev_tier = None
    for row in sorted(inter, key=_pair_sort_key):
        i1, i2 = str(row['Image1']), str(row['Image2'])
        f1, f2 = str(row['Filter1']), str(row['Filter2'])
        tier   = 'short' if (i1.endswith('_s') or i2.endswith('_s')) else 'long'
        if prev_tier is not None and tier != prev_tier:
            print()
        prev_tier = tier

        # r-band first; otherwise shorter wavelength first
        if f2 == 'r' or (_fwave(f2) < _fwave(f1) and f1 != 'r'):
            i1, i2, f1, f2 = i2, i1, f2, f1
        label = f'{i1} × {i2}'

        mz  = float(row['sigma_c_magzero'])   if 'sigma_c_magzero'   in inter.colnames else np.nan
        gg  = float(row['sigma_c_gaia_g'])    if 'sigma_c_gaia_g'    in inter.colnames else np.nan
        sm  = float(row['sigma_c_smash_r'])   if 'sigma_c_smash_r'   in inter.colnames else np.nan
        res = float(row['sigma_residual'])     if 'sigma_residual'    in inter.colnames else np.nan
        pct = float(row['residual_pct'])       if 'residual_pct'      in inter.colnames else np.nan
        dmz = float(row['std_dzp_magzero'])   if 'std_dzp_magzero'   in inter.colnames else np.nan
        dgg = float(row['std_dzp_gaia_g'])    if 'std_dzp_gaia_g'    in inter.colnames else np.nan
        dsm = float(row['std_dzp_smash_r'])   if 'std_dzp_smash_r'   in inter.colnames else np.nan
        ns  = int(row['N_stars'])              if 'N_stars'           in inter.colnames else 0

        print(f'  {label:<22} {ns:>8}  '
              f'{mz:>8.4f} {gg:>8.4f} {sm:>8.4f}  '
              f'{res:>8.4f}  {pct:>6.1f}%  '
              f'{dmz:>8.4f} {dgg:>8.4f} {dsm:>8.4f}')
    print()


def summarize(fits_file):
    """Read a ``*_phot_check.fits`` file and print the comparison summary.

    Can be called as::

        CheckPhot.py -summary CheckPhot/LMC_c42_phot_check.fits
    """
    from astropy.io import fits as _fits
    from astropy.table import Table as _Table

    with _fits.open(fits_file) as hdul:
        field = hdul[0].header.get('FIELD', os.path.basename(fits_file))
        intra = _Table(hdul['INTRA'].data)
        inter = _Table(hdul['INTER'].data)

    _print_summary(intra, inter, field)


def _find_config():
    """Locate DeMCELS_images.txt: current directory first, then $KRED/config."""
    if os.path.isfile(_CONFIG_NAME):
        return _CONFIG_NAME
    kred = os.environ.get('KRED', '')
    candidate = os.path.join(kred, 'config', _CONFIG_NAME)
    if os.path.isfile(candidate):
        return candidate
    raise FileNotFoundError(
        f'Cannot find {_CONFIG_NAME} in current directory or $KRED/config/')


def _parse_image_config(config_file):
    """Parse DeMCELS_images.txt.

    Returns {label: {'filter': str, 'exptimes': set of float, 'tier': str}}
    where tier='short' if label ends with '_s', else 'long'.
    """
    groups = {}
    with open(config_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('-') or line.lower().startswith('image'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            label, filt, exptime = parts[0], parts[1], float(parts[2])
            if label not in groups:
                groups[label] = {
                    'filter':   filt,
                    'exptimes': set(),
                    'tier':     'short' if label.endswith('_s') else 'long',
                }
            groups[label]['exptimes'].add(exptime)
    return groups


def _assign_files_to_groups(tabphot_files, image_groups):
    """Read FILTER and EXPTIME from each TabPhot header and assign to Image groups.

    Returns (group_smash, group_gaia) where each is {label: [filepath, ...]}.
    """
    from astropy.io import fits as _fits

    group_smash = {label: [] for label in image_groups}
    group_gaia  = {label: [] for label in image_groups}
    n_skip = 0

    for f in tabphot_files:
        try:
            with _fits.open(f, memmap=False) as hdul:
                hdr  = hdul[1].header
                data = hdul[1].data
                # Header takes priority; fall back to first-row column value
                # for old-format files that store metadata as per-row columns
                if 'FILTER' in hdr:
                    # DECam FILTER header is "N662 DECam c0009 6620.0 172.0";
                    # take only the first token (the filter label)
                    filt = str(hdr['FILTER']).strip().split()[0]
                elif data is not None and 'Filter' in data.dtype.names:
                    filt = str(data['Filter'][0]).strip().split()[0]
                else:
                    filt = ''
                if 'EXPTIME' in hdr:
                    exptime = float(hdr['EXPTIME'])
                elif data is not None and 'Exptime' in data.dtype.names:
                    exptime = float(data['Exptime'][0])
                else:
                    exptime = -1
        except Exception:
            n_skip += 1
            continue

        matched = False
        for label, ginfo in image_groups.items():
            if (ginfo['filter'] == filt and
                    any(abs(exptime - e) < 5.0 for e in ginfo['exptimes'])):
                if f.endswith('.smash.fits'):
                    group_smash[label].append(f)
                else:
                    group_gaia[label].append(f)
                matched = True
                break
        if not matched:
            n_skip += 1

    if n_skip:
        print(f'CheckPhot: {n_skip} TabPhot files not matched to any Image group')
    return group_smash, group_gaia


def _auto_image_pairs(active_groups):
    """All (label1, label2) pairs in the same tier with different filters."""
    from itertools import combinations

    tier_labels = {}
    for label, ginfo in active_groups.items():
        tier_labels.setdefault(ginfo['tier'], []).append(label)

    pairs = []
    for labels in tier_labels.values():
        for l1, l2 in combinations(sorted(labels), 2):
            if active_groups[l1]['filter'] != active_groups[l2]['filter']:
                pairs.append((l1, l2))
    return pairs


def _build_zp_lookup(zp_tab):
    """Build the zp_lookup dict expected by PhotEval.do_eval().

    Keys are TabPhot file basenames (e.g. 'c4d_..._req.smash.fits');
    values are the empirical zp_calc from CalcZeroPoint.
    """
    if zp_tab is None or len(zp_tab) == 0:
        return {}
    lookup = {}
    for row in zp_tab:
        fname = os.path.basename(str(row['Filename']).strip())
        try:
            lookup[fname] = float(row['zp_calc'])
        except (KeyError, ValueError):
            pass
    return lookup


def _zp_stats_by_key(zp_tab):
    """Split a CalcZeroPoint result table by (catalog, ref_col) and compute stats.

    Returns dict mapping key (e.g. 'smash_r', 'gaia_g') to
    {'N_files', 'N_zp', 'mean_dzp', 'std_dzp'}.
    """
    if zp_tab is None or len(zp_tab) == 0:
        return {}

    cats   = [str(c).lower() for c in zp_tab['Catalog']]
    refs   = [str(r).lower() for r in zp_tab['ref_col']]
    deltas = np.array(zp_tab['delta_zp'], dtype=float)
    nstars = np.array(zp_tab['n_stars'])

    buckets = defaultdict(list)
    for i, (cat, ref) in enumerate(zip(cats, refs)):
        buckets[f'{cat}_{ref}'].append(i)

    out = {}
    for key, idxs in buckets.items():
        d  = deltas[idxs]
        ns = nstars[idxs]
        ok = np.isfinite(d)
        if ok.sum() == 0:
            continue
        out[key] = {
            'N_files':  int(ok.sum()),
            'N_zp':     int(np.nansum(ns[ok])),
            'mean_dzp': float(np.mean(d[ok])),
            'std_dzp':  float(mad_std(d[ok])) if ok.sum() > 1 else np.nan,
        }
    return out


def check_phot(field, mef_dir='DECam_MEF', tabphot_dir='TabPhot',
               config_file=None, np_proc=8, redo=False,
               run_mefphot=True, checkdir=None, extra_refs=None):
    """Run the full photometric consistency check for *field*.

    Parameters
    ----------
    field : str
        Field name, e.g. ``'LMC_c42'``.
    mef_dir : str
        Directory containing ``{field}/mef/*.fits.fz`` (default ``DECam_MEF``).
    tabphot_dir : str
        Directory where MefPhot writes its output (default ``TabPhot``).
    config_file : str or None
        Path to DeMCELS_images.txt.  Searched in current directory then
        $KRED/config/ if not given.
    np_proc : int
        MefPhot parallel workers (default 8).
    redo : bool
        Reprocess existing TabPhot files (default False).
    run_mefphot : bool
        If False, skip MefPhot and assume TabPhot/ is already populated.
    checkdir : str or None
        Directory for all CheckPhot outputs (intermediate and summary).
        Default: ``CheckPhot``.
    extra_refs : list of str or None
        Additional catalog column names to compare against in CalcZeroPoint.

    Returns
    -------
    intra : astropy.table.Table
        One row per Image group with intra-filter consistency statistics.
    inter : astropy.table.Table
        One row per Image-group pair with inter-filter consistency statistics.
    """
    if config_file is None:
        config_file = _find_config()
    if checkdir is None:
        checkdir = 'CheckPhot'
    os.makedirs(checkdir, exist_ok=True)

    # Parse image groups from config
    image_groups = _parse_image_config(config_file)
    print(f'CheckPhot: parsed {len(image_groups)} Image groups from {config_file}')

    # Assign existing TabPhot files to Image groups
    all_tabphot = (sorted(glob(os.path.join(tabphot_dir, '*.smash.fits'))) +
                   sorted(glob(os.path.join(tabphot_dir, '*.gaia.fits'))))
    group_smash, group_gaia = _assign_files_to_groups(all_tabphot, image_groups)

    # Active groups: those with at least one file present
    active_groups = {lb: gi for lb, gi in image_groups.items()
                     if group_smash.get(lb) or group_gaia.get(lb)}
    if not active_groups:
        print('CheckPhot: no TabPhot files matched any Image group — '
              'run MefPhot first or check the config file.')
        return None, None

    print(f'CheckPhot: {len(active_groups)} active Image groups for {field}:')
    for lb, gi in active_groups.items():
        ns = len(group_smash.get(lb, []))
        ng = len(group_gaia.get(lb, []))
        print(f'  {lb:12s}  filter={gi["filter"]:6s}  tier={gi["tier"]:5s}  '
              f'smash={ns}  gaia={ng}')

    image_pairs = _auto_image_pairs(active_groups)
    print(f'CheckPhot: {len(image_pairs)} inter-image pairs:',
          ', '.join(f'{a}×{b}' for a, b in image_pairs))

    # ------------------------------------------------------------------
    # Step 1: MefPhot
    # ------------------------------------------------------------------
    if run_mefphot:
        mef_glob = f'{mef_dir}/{field}/mef/*.fits.fz'
        mef_files = sorted(glob(mef_glob))
        if not mef_files:
            print(f'CheckPhot: no MEF files found matching {mef_glob}')
        else:
            base = ['MefPhot.py', '-r', '6', '-np', str(np_proc)]
            if redo:
                base.append('-redo')

            print(f'\n{"="*60}')
            print(f'MefPhot SMASH  ({len(mef_files)} files)')
            print(f'{"="*60}')
            MefPhot.steer(base + ['-cat', 'smash'] + mef_files)

            print(f'\n{"="*60}')
            print(f'MefPhot Gaia   ({len(mef_files)} files)')
            print(f'{"="*60}')
            MefPhot.steer(base + mef_files)

    # ------------------------------------------------------------------
    # Step 2: CalcZeroPoint — collect ZP stats split by (catalog, ref_col)
    # ------------------------------------------------------------------
    _NAMED_REFS = {'G': 'gaia_g', 'R': 'smash_r'}
    extra_ref_list = list(extra_refs or [])
    refs_to_run = (
        [None, 'G', 'R'] +
        [r for r in extra_ref_list if r not in _NAMED_REFS]
    )

    # all_zp[label][key] = {N_files, N_zp, mean_dzp, std_dzp}
    all_zp              = {lb: {} for lb in active_groups}
    all_keys            = []
    all_zp_lookup       = {}   # label -> auto zp_lookup for PhotEval
    all_named_zp_lookup = {}   # label -> {'gaia_g': {...}, 'smash_r': {...}}

    for lb, gi in active_groups.items():
        filt = gi['filter']
        all_named_zp_lookup[lb] = {}
        lb_files = sorted(group_smash.get(lb, []) + group_gaia.get(lb, []))
        for idx, ref_override in enumerate(refs_to_run):
            rstr = f'ref_col={ref_override}' if ref_override else 'auto'
            print(f'\n{"="*60}')
            print(f'CalcZeroPoint  image={lb}  filter={filt}  ({rstr})')
            print(f'{"="*60}')
            ref_tag = f'_{ref_override.lower()}' if ref_override else ''
            zp_out  = f'{checkdir}/{field}_zeropoints_{lb}{ref_tag}.fits'
            zp_tab  = CalcZeroPoint.run(files=lb_files,
                                        ref_col_override=ref_override,
                                        outfile=zp_out)
            if idx == 0:
                all_zp_lookup[lb] = _build_zp_lookup(zp_tab)
            if ref_override in _NAMED_REFS:
                all_named_zp_lookup[lb][_NAMED_REFS[ref_override]] = _build_zp_lookup(zp_tab)
            for key, stats in _zp_stats_by_key(zp_tab).items():
                if key not in all_zp[lb]:
                    all_zp[lb][key] = stats
                if key not in all_keys:
                    all_keys.append(key)

    # ------------------------------------------------------------------
    # Step 3: PhotEval per Image group
    # ------------------------------------------------------------------
    all_eval = {}
    for lb, gi in active_groups.items():
        filt = gi['filter']
        print(f'\n{"="*60}')
        print(f'PhotEval       image={lb}  filter={filt}')
        print(f'{"="*60}')
        named_lkp  = all_named_zp_lookup.get(lb) or None
        eval_files = group_smash.get(lb) or group_gaia.get(lb) or []
        eval_tab = PhotEval.do_eval(
            tabphot_dir=tabphot_dir,
            filter_str=filt,
            exptime=None,
            snr_min=10,
            prob_min=0.5,
            min_n=3,
            outfile=f'{checkdir}/{field}_phot_eval_{lb}.fits',
            zp_lookup=all_zp_lookup.get(lb),
            extra_zp_lookups=named_lkp,
            files=eval_files if eval_files else None,
        )
        er = {'N_sources': 0,
              'chi2_nu_c_med': np.nan, 'chi2_nu_c_90': np.nan,
              'chi2_nu_c_emp_med': np.nan, 'chi2_nu_c_emp_90': np.nan,
              'sigma_c_magzero': np.nan,
              'sigma_c_gaia_g':  np.nan,
              'sigma_c_smash_r': np.nan}
        if eval_tab is not None and len(eval_tab) > 0:
            if 'chi2_nu_c' in eval_tab.colnames:
                c = np.array(eval_tab['chi2_nu_c'], dtype=float)
                c = c[np.isfinite(c)]
                er['N_sources']     = len(c)
                er['chi2_nu_c_med'] = float(np.median(c))         if len(c) else np.nan
                er['chi2_nu_c_90']  = float(np.percentile(c, 90)) if len(c) else np.nan
            if 'chi2_nu_emp' in eval_tab.colnames:
                ce = np.array(eval_tab['chi2_nu_emp'], dtype=float)
                ce = ce[np.isfinite(ce)]
                er['chi2_nu_c_emp_med'] = float(np.median(ce))         if len(ce) else np.nan
                er['chi2_nu_c_emp_90']  = float(np.percentile(ce, 90)) if len(ce) else np.nan

            # sigma_c: median per-star scatter for bright stars (photon noise < 10 mmag)
            # Bright-star cut makes sigma_c ≈ calibration floor, not photon noise
            err_arr = np.array(eval_tab['mag_err_mean'], dtype=float)
            bright  = np.isfinite(err_arr) & (err_arr < 0.010)

            def _sc(col):
                if col not in eval_tab.colnames:
                    return np.nan
                v = np.array(eval_tab[col], dtype=float)[bright]
                v = v[np.isfinite(v)]
                return float(np.median(v)) if len(v) >= 10 else np.nan

            er['sigma_c_magzero'] = _sc('magc_std')
            er['sigma_c_gaia_g']  = _sc('magc_std_gaia_g')
            er['sigma_c_smash_r'] = _sc('magc_std_smash_r')

        all_eval[lb] = er

    # ------------------------------------------------------------------
    # Build intra-filter table  (one row per Image group)
    # ------------------------------------------------------------------
    colnames = ['Image', 'Filter']
    for key in all_keys:
        colnames += [f'N_files_{key}', f'N_zp_{key}',
                     f'mean_dzp_{key}', f'std_dzp_{key}']
    colnames += ['N_sources', 'chi2_nu_c_med', 'chi2_nu_c_90',
                 'chi2_nu_c_emp_med', 'chi2_nu_c_emp_90',
                 'sigma_c_magzero', 'sigma_c_gaia_g', 'sigma_c_smash_r']

    intra_rows = []
    for lb, gi in active_groups.items():
        filt = gi['filter']
        row = {'Image': lb, 'Filter': filt}
        for key in all_keys:
            ks = all_zp[lb].get(key, {})
            row[f'N_files_{key}']  = ks.get('N_files',  0)
            row[f'N_zp_{key}']     = ks.get('N_zp',     0)
            row[f'mean_dzp_{key}'] = ks.get('mean_dzp', np.nan)
            row[f'std_dzp_{key}']  = ks.get('std_dzp',  np.nan)
        er = all_eval.get(lb, {})
        row['N_sources']         = er.get('N_sources',         0)
        row['chi2_nu_c_med']     = er.get('chi2_nu_c_med',     np.nan)
        row['chi2_nu_c_90']      = er.get('chi2_nu_c_90',      np.nan)
        row['chi2_nu_c_emp_med'] = er.get('chi2_nu_c_emp_med', np.nan)
        row['chi2_nu_c_emp_90']  = er.get('chi2_nu_c_emp_90',  np.nan)
        row['sigma_c_magzero']   = er.get('sigma_c_magzero',   np.nan)
        row['sigma_c_gaia_g']    = er.get('sigma_c_gaia_g',    np.nan)
        row['sigma_c_smash_r']   = er.get('sigma_c_smash_r',   np.nan)
        intra_rows.append(row)

    intra = Table(rows=intra_rows, names=colnames)
    for key in all_keys:
        for col in (f'mean_dzp_{key}', f'std_dzp_{key}'):
            intra[col].format = '.4f'
    for col in ('chi2_nu_c_med', 'chi2_nu_c_90',
                'chi2_nu_c_emp_med', 'chi2_nu_c_emp_90',
                'sigma_c_magzero', 'sigma_c_gaia_g', 'sigma_c_smash_r'):
        intra[col].format = '.4f'

    _intra_desc = {
        'Image':              'Image-group label from DeMCELS_images.txt (e.g. N662, r_s)',
        'Filter':             'Filter name',
        'N_files_magzero':    'Files with MAGZERO-based ZP measurement',
        'N_zp_magzero':       'Stars used for MAGZERO-based ZP',
        'mean_dzp_magzero':   'Mean(measured_ZP - header MAGZERO) across files [mag]',
        'std_dzp_magzero':    'MAD scatter of (measured_ZP - MAGZERO) across files [mag]',
        'N_files_gaia_g':     'Files with Gaia-G-based ZP measurement',
        'N_zp_gaia_g':        'Stars used for Gaia-G-based ZP',
        'mean_dzp_gaia_g':    'Mean(ZP_GaiaG - header MAGZERO) across files [mag]',
        'std_dzp_gaia_g':     'MAD scatter of (ZP_GaiaG - MAGZERO) across files [mag]',
        'N_files_smash_r':    'Files with SMASH-R-based ZP measurement',
        'N_zp_smash_r':       'Stars used for SMASH-R-based ZP',
        'mean_dzp_smash_r':   'Mean(ZP_SmashR - header MAGZERO) across files [mag]',
        'std_dzp_smash_r':    'MAD scatter of (ZP_SmashR - MAGZERO) across files [mag]',
        'N_sources':          'Unique sources detected across all frames',
        'chi2_nu_c_med':      'Median reduced chi-squared of MAGZERO-corrected mag across frames',
        'chi2_nu_c_90':       '90th-percentile reduced chi-squared (MAGZERO)',
        'chi2_nu_c_emp_med':  'Median reduced chi-squared using empirical ZP',
        'chi2_nu_c_emp_90':   '90th-percentile reduced chi-squared (empirical ZP)',
        'sigma_c_magzero':    'Median per-star scatter in MAGZERO-corrected mag, bright stars [mag]',
        'sigma_c_gaia_g':     'Median per-star scatter using Gaia G ZP, bright stars [mag]',
        'sigma_c_smash_r':    'Median per-star scatter using SMASH R ZP, bright stars [mag]',
    }
    for col, desc in _intra_desc.items():
        if col in intra.colnames:
            intra[col].description = desc

    # ------------------------------------------------------------------
    # Step 4: PhotAnal per Image-group pair → inter-filter table
    # ------------------------------------------------------------------
    # Combined named ZP lookups across all Image groups for PhotAnal
    combined_named = {}
    for lb in active_groups:
        for name, lkup in all_named_zp_lookup.get(lb, {}).items():
            combined_named.setdefault(name, {}).update(lkup)
    named_for_anal = combined_named if combined_named else None

    # Pre-load gaia data for each unique Image group that appears in a pair.
    # Passing pre-filtered gaia files avoids loading the wrong exptime.
    unique_labels = list(dict.fromkeys(lb for pair in image_pairs for lb in pair))
    filter_data_cache = {}
    for lb in unique_labels:
        filt = active_groups[lb]['filter']
        lb_gaia = group_gaia.get(lb, [])
        if not lb_gaia:
            print(f'CheckPhot: no gaia files for Image group {lb} — '
                  f'skipping PhotAnal pairs involving {lb}')
            continue
        print(f'\n{"="*60}')
        print(f'PhotAnal: loading {lb} ({filt}) data...')
        print(f'{"="*60}')
        filter_data_cache[lb] = PhotAnal.load_filter_data(
            lb_gaia, filt, extra_zp_lookups=named_for_anal)

    inter_rows = []
    for l1, l2 in image_pairs:
        if l1 not in filter_data_cache or l2 not in filter_data_cache:
            print(f'CheckPhot: skipping {l1}×{l2} (no gaia data)')
            continue
        f1 = active_groups[l1]['filter']
        f2 = active_groups[l2]['filter']
        print(f'\n{"="*60}')
        print(f'PhotAnal       {l1} ({f1}) vs {l2} ({f2})')
        print(f'{"="*60}')
        result = PhotAnal.do_compare(
            group_gaia.get(l1, []), f1, f2,
            outfile=f'{checkdir}/{field}_filter_compare_{l1}_{l2}.fits',
            extra_zp_lookups=named_for_anal,
            _preloaded=(filter_data_cache[l1], filter_data_cache[l2]),
        )

        row = {'Image1': l1, 'Image2': l2, 'Filter1': f1, 'Filter2': f2}
        if result is not None and len(result) > 0:
            m = result.meta
            std_res = float(m.get('STD_RES', np.nan))
            row['N_stars']            = len(result)
            row['mean_delta_mag']     = float(m.get('MEAN_DM',  np.nan))
            row['sigma_c_magzero']    = float(m.get('STD_DM',   np.nan))
            row['sigma_residual']     = std_res
            row['residual_pct']       = float((10 ** (abs(std_res) / 2.5) - 1) * 100) \
                                        if np.isfinite(std_res) else np.nan
            row['color_term_b']       = float(m.get('COLOR_B',  np.nan))
            row['N_pairs']            = int(m.get('N_PAIRS',    0))
            row['std_dzp_magzero']    = float(m.get('STD_PDM',  np.nan))
            row['sigma_c_gaia_g']     = float(m.get('SDTM_GA',  np.nan))
            row['sigma_resid_gaia_g'] = float(m.get('STRE_GA',  np.nan))
            row['std_dzp_gaia_g']     = float(m.get('SPDT_GA',  np.nan))
            row['sigma_c_smash_r']    = float(m.get('SDTM_SM',  np.nan))
            row['sigma_resid_smash_r']= float(m.get('STRE_SM',  np.nan))
            row['std_dzp_smash_r']    = float(m.get('SPDT_SM',  np.nan))
        else:
            for k in ('N_stars', 'N_pairs'):
                row[k] = 0
            for k in ('mean_delta_mag', 'sigma_c_magzero', 'sigma_residual',
                      'residual_pct', 'color_term_b', 'std_dzp_magzero',
                      'sigma_c_gaia_g', 'sigma_resid_gaia_g', 'std_dzp_gaia_g',
                      'sigma_c_smash_r', 'sigma_resid_smash_r', 'std_dzp_smash_r'):
                row[k] = np.nan
        inter_rows.append(row)

    inter = Table(rows=inter_rows,
                  names=['Image1', 'Image2', 'Filter1', 'Filter2', 'N_stars',
                         'mean_delta_mag',
                         'sigma_c_magzero', 'sigma_residual', 'residual_pct',
                         'color_term_b', 'N_pairs', 'std_dzp_magzero',
                         'sigma_c_gaia_g', 'sigma_resid_gaia_g', 'std_dzp_gaia_g',
                         'sigma_c_smash_r', 'sigma_resid_smash_r', 'std_dzp_smash_r'])
    for col in ('mean_delta_mag',
                'sigma_c_magzero', 'sigma_residual', 'residual_pct',
                'color_term_b', 'std_dzp_magzero',
                'sigma_c_gaia_g', 'sigma_resid_gaia_g', 'std_dzp_gaia_g',
                'sigma_c_smash_r', 'sigma_resid_smash_r', 'std_dzp_smash_r'):
        inter[col].format = '.4f'

    _inter_desc = {
        'Image1':              'Reference Image group (from DeMCELS_images.txt)',
        'Image2':              'Comparison Image group (from DeMCELS_images.txt)',
        'Filter1':             'Filter of Image1',
        'Filter2':             'Filter of Image2',
        'N_stars':             'Stars matched between Image1 and Image2',
        'mean_delta_mag':      'Mean(magc_Image1 - magc_Image2) across matched stars [mag]',
        'sigma_c_magzero':     'MAD scatter of (magc_F1 - magc_F2) using MAGZERO [mag]',
        'sigma_residual':      'Scatter after fitting and removing color term (MAGZERO) [mag]',
        'residual_pct':        'sigma_residual expressed as percentage of flux',
        'color_term_b':        'Color term coefficient b in delta_mag = a + b*(G-R)',
        'N_pairs':             'Number of (Image1 file, Image2 file) pairs compared',
        'std_dzp_magzero':     'MAD scatter of per-pair mean delta_mag using MAGZERO [mag]',
        'sigma_c_gaia_g':      'MAD scatter of (magc_F1 - magc_F2) using Gaia G ZP [mag]',
        'sigma_resid_gaia_g':  'Scatter after color term removal using Gaia G ZP [mag]',
        'std_dzp_gaia_g':      'MAD scatter of per-pair mean delta_mag using Gaia G ZP [mag]',
        'sigma_c_smash_r':     'MAD scatter of (magc_F1 - magc_F2) using SMASH R ZP [mag]',
        'sigma_resid_smash_r': 'Scatter after color term removal using SMASH R ZP [mag]',
        'std_dzp_smash_r':     'MAD scatter of per-pair mean delta_mag using SMASH R ZP [mag]',
    }
    for col, desc in _inter_desc.items():
        if col in inter.colnames:
            inter[col].description = desc

    _print_summary(intra, inter, field)

    # ------------------------------------------------------------------
    # Write single FITS with two named table extensions
    # ------------------------------------------------------------------
    from astropy.io import fits as _fits
    primary = _fits.PrimaryHDU()
    primary.header['FIELD'] = (field, 'Field name')
    intra_hdu      = _fits.table_to_hdu(intra)
    intra_hdu.name = 'INTRA'
    inter_hdu      = _fits.table_to_hdu(inter)
    inter_hdu.name = 'INTER'
    outfits = f'{checkdir}/{field}_phot_check.fits'
    _fits.HDUList([primary, intra_hdu, inter_hdu]).writeto(outfits, overwrite=True)
    print(f'Wrote {outfits}  (ext INTRA + INTER)')

    return intra, inter


def steer(argv):
    field        = None
    np_proc      = 8
    redo         = False
    run_mefphot  = True
    checkdir     = None
    config_file  = None
    extra_refs   = []
    summary_file = None

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(__doc__)
            return
        elif arg == '-summary':
            i += 1
            summary_file = argv[i]
        elif arg == '-np':
            i += 1
            np_proc = int(argv[i])
        elif arg == '-redo':
            redo = True
        elif arg == '-no_mefphot':
            run_mefphot = False
        elif arg == '-ref':
            i += 1
            extra_refs.append(argv[i])
        elif arg == '-o':
            i += 1
            checkdir = argv[i]
        elif arg == '-config':
            i += 1
            config_file = argv[i]
        elif arg.startswith('-'):
            print(f'CheckPhot: unknown flag {arg}')
            return
        else:
            field = arg
        i += 1

    if summary_file is not None:
        summarize(summary_file)
        return

    if field is None:
        print('CheckPhot: field name required')
        print(__doc__)
        return

    check_phot(field, np_proc=np_proc, redo=redo,
               run_mefphot=run_mefphot, checkdir=checkdir,
               config_file=config_file,
               extra_refs=extra_refs or None)


if __name__ == '__main__':
    try:
        if len(sys.argv) > 1:
            steer(sys.argv)
        else:
            print(__doc__)
    except KeyboardInterrupt:
        print('\nCheckPhot: interrupted')
