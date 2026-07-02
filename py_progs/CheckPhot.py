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
    CheckPhot.py -no_mefphot LMC_c42     # skip MefPhot if TabPhot/ already populated
    CheckPhot.py -redo LMC_c42           # force reprocessing of all TabPhot files

Flags
-----

-np N           Number of parallel MefPhot workers (default 8)
-redo           Reprocess already-existing TabPhot files
-no_mefphot     Skip the MefPhot step
-o ROOT         Root name for output FITS tables
                (default: Summary/{field}_phot_check)
-h              Print this help

Description
-----------

Four steps are run in sequence:

1. **MefPhot** (twice) — forced aperture photometry at SMASH positions
   (for PhotEval) and at Gaia positions (for PhotAnal).  Already-processed
   files are skipped automatically unless ``-redo`` is given.

2. **CalcZeroPoint** — per-file zero-point fit vs catalog for each filter.
   Measures how far each file's header MAGZERO deviates from the measured ZP.

3. **PhotEval** — groups detections of the same star across overlapping
   frames and computes ``chi2_nu_c``: the reduced chi-squared of the
   MAGZERO-corrected magnitudes relative to photon noise.

4. **PhotAnal** — cross-matches the same star between two filters and
   measures the scatter in (mag_filter1 − mag_filter2), which sets the
   minimum stellar residual in CleanStars output.

Output
------

``{ROOT}_intra.fits`` — one row per filter:

  Filter, N_files, N_stars_zp, mean_delta_zp, std_delta_zp,
  N_sources, chi2_nu_c_med, chi2_nu_c_90

``{ROOT}_inter.fits`` — one row per filter pair:

  Filter1, Filter2, N_stars, mean_delta_mag,
  sigma_total, sigma_residual, residual_pct, color_term_b

Interpreting the tables
-----------------------

Intra-filter:
  ``chi2_nu_c_med ≈ 1`` means frames are mutually consistent at the
  photon-noise level; no further calibration is needed.
  ``std_delta_zp`` is the MAD-based scatter of per-file MAGZERO deviations
  from the measured zero point (robust to outlier exposures).

Inter-filter:
  ``residual_pct`` is the minimum percentage of stellar flux that will
  remain as a residual after CleanStars continuum subtraction,
  assuming the best-fit mean scale factor between the two filters.
  ``sigma_residual`` is the scatter after removing the color-term trend;
  this is the irreducible floor set by stellar color diversity.

History
-------

260701 ksl  Written
"""

import sys
import os
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.stats import mad_std

import MefPhot
import CalcZeroPoint
import PhotEval
import PhotAnal


FILTERS_DEFAULT    = ['r', 'N662', 'N673', 'N708']
PAIRS_DEFAULT      = [('r', 'N662'), ('r', 'N673'),
                      ('N708', 'N662'), ('N708', 'N673')]


def check_phot(field, mef_dir='DECam_MEF', tabphot_dir='TabPhot',
               filters=None, filter_pairs=None, np_proc=8, redo=False,
               run_mefphot=True, outroot=None):
    """Run the full photometric consistency check for *field*.

    Parameters
    ----------
    field : str
        Field name, e.g. ``'LMC_c42'``.
    mef_dir : str
        Directory containing ``{field}/mef/*.fits.fz`` (default ``DECam_MEF``).
    tabphot_dir : str
        Directory where MefPhot writes its output (default ``TabPhot``).
    filters : list of str or None
        Filters to check individually.  Default: r, N662, N673, N708.
    filter_pairs : list of (str, str) or None
        Filter pairs for inter-filter comparison.
        Default: (r,N662), (r,N673), (N708,N662), (N708,N673).
    np_proc : int
        MefPhot parallel workers (default 8).
    redo : bool
        Reprocess existing TabPhot files (default False).
    run_mefphot : bool
        If False, skip MefPhot and assume TabPhot/ is already populated.
    outroot : str or None
        Root path for output FITS tables.
        Default: ``Summary/{field}_phot_check``.

    Returns
    -------
    intra : astropy.table.Table
        One row per filter with intra-filter consistency statistics.
    inter : astropy.table.Table
        One row per filter pair with inter-filter consistency statistics.
    """
    if filters is None:
        filters = FILTERS_DEFAULT
    if filter_pairs is None:
        filter_pairs = PAIRS_DEFAULT
    if outroot is None:
        os.makedirs('Summary', exist_ok=True)
        outroot = f'Summary/{field}_phot_check'

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
    # Step 2 & 3: CalcZeroPoint + PhotEval per filter → intra-filter table
    # ------------------------------------------------------------------
    intra_rows = []
    for filt in filters:
        print(f'\n{"="*60}')
        print(f'CalcZeroPoint  filter={filt}')
        print(f'{"="*60}')
        zp_tab = CalcZeroPoint.run(filter_str=filt, tabphot_dir=tabphot_dir)

        print(f'\n{"="*60}')
        print(f'PhotEval       filter={filt}')
        print(f'{"="*60}')
        eval_tab = PhotEval.do_eval(
            tabphot_dir=tabphot_dir,
            filter_str=filt,
            exptime=None,
            snr_min=10,
            prob_min=0.5,
            min_n=3,
            outfile=f'phot_eval_{filt}.fits',
        )

        row = {'Filter': filt}

        if zp_tab is not None and len(zp_tab) > 0:
            delta = np.array(zp_tab['delta_zp'], dtype=float)
            ok = np.isfinite(delta)
            row['N_files']      = int(np.sum(ok))
            row['N_stars_zp']   = int(np.nansum(np.array(zp_tab['n_stars'])[ok]))
            row['mean_delta_zp'] = float(np.mean(delta[ok]))
            row['std_delta_zp']  = float(mad_std(delta[ok])) if ok.sum() > 1 else np.nan
        else:
            row['N_files']       = 0
            row['N_stars_zp']    = 0
            row['mean_delta_zp'] = np.nan
            row['std_delta_zp']  = np.nan

        if eval_tab is not None and len(eval_tab) > 0 and 'chi2_nu_c' in eval_tab.colnames:
            c = np.array(eval_tab['chi2_nu_c'], dtype=float)
            c = c[np.isfinite(c)]
            row['N_sources']     = len(c)
            row['chi2_nu_c_med'] = float(np.median(c))    if len(c) else np.nan
            row['chi2_nu_c_90']  = float(np.percentile(c, 90)) if len(c) else np.nan
        else:
            row['N_sources']     = 0
            row['chi2_nu_c_med'] = np.nan
            row['chi2_nu_c_90']  = np.nan

        intra_rows.append(row)

    intra = Table(rows=intra_rows,
                  names=['Filter', 'N_files', 'N_stars_zp',
                         'mean_delta_zp', 'std_delta_zp',
                         'N_sources', 'chi2_nu_c_med', 'chi2_nu_c_90'])
    for col in ('mean_delta_zp', 'std_delta_zp', 'chi2_nu_c_med', 'chi2_nu_c_90'):
        intra[col].format = '.4f'

    # ------------------------------------------------------------------
    # Step 4: PhotAnal per filter pair → inter-filter table
    # ------------------------------------------------------------------
    gaia_files = sorted(glob(os.path.join(tabphot_dir, '*.gaia.fits')))

    inter_rows = []
    for f1, f2 in filter_pairs:
        print(f'\n{"="*60}')
        print(f'PhotAnal       {f1} vs {f2}')
        print(f'{"="*60}')
        result = PhotAnal.do_compare(
            gaia_files, f1, f2,
            outfile=f'filter_compare_{f1}_{f2}.fits',
        )

        row = {'Filter1': f1, 'Filter2': f2}
        if result is not None and len(result) > 0:
            m = result.meta
            std_res = float(m.get('STD_RES', np.nan))
            row['N_stars']       = len(result)
            row['mean_delta_mag'] = float(m.get('MEAN_DM',  np.nan))
            row['sigma_total']   = float(m.get('STD_DM',   np.nan))
            row['sigma_residual'] = std_res
            row['residual_pct']  = float((10 ** (abs(std_res) / 2.5) - 1) * 100) \
                                   if np.isfinite(std_res) else np.nan
            row['color_term_b']  = float(m.get('COLOR_B',  np.nan))
        else:
            row['N_stars']        = 0
            row['mean_delta_mag'] = np.nan
            row['sigma_total']    = np.nan
            row['sigma_residual'] = np.nan
            row['residual_pct']   = np.nan
            row['color_term_b']   = np.nan
        inter_rows.append(row)

    inter = Table(rows=inter_rows,
                  names=['Filter1', 'Filter2', 'N_stars',
                         'mean_delta_mag', 'sigma_total',
                         'sigma_residual', 'residual_pct', 'color_term_b'])
    for col in ('mean_delta_mag', 'sigma_total', 'sigma_residual',
                'residual_pct', 'color_term_b'):
        inter[col].format = '.4f'

    # ------------------------------------------------------------------
    # Print concise summary
    # ------------------------------------------------------------------
    bar = '=' * 65
    print(f'\n\n{bar}')
    print(f'  PHOTOMETRY CHECK: {field}')
    print(bar)

    print('\nIntra-filter  (chi2_nu_c ~ 1 → frames consistent at photon noise)')
    print(f'  std_delta_zp = MAD scatter of per-file MAGZERO offsets (robust)\n')
    intra.pprint(max_lines=-1, max_width=-1)

    print('\nInter-filter  (residual_pct = minimum stellar residual in CleanStars output)')
    print(f'  sigma_residual = scatter after color correction [irreducible floor]\n')
    inter.pprint(max_lines=-1, max_width=-1)
    print()

    # ------------------------------------------------------------------
    # Write output tables
    # ------------------------------------------------------------------
    intra_file = f'{outroot}_intra.fits'
    inter_file = f'{outroot}_inter.fits'
    intra.write(intra_file, overwrite=True)
    inter.write(inter_file, overwrite=True)
    print(f'Wrote {intra_file}')
    print(f'Wrote {inter_file}')

    return intra, inter


def steer(argv):
    field        = None
    np_proc      = 8
    redo         = False
    run_mefphot  = True
    outroot      = None

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(__doc__)
            return
        elif arg == '-np':
            i += 1
            np_proc = int(argv[i])
        elif arg == '-redo':
            redo = True
        elif arg == '-no_mefphot':
            run_mefphot = False
        elif arg == '-o':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print(f'CheckPhot: unknown flag {arg}')
            return
        else:
            field = arg
        i += 1

    if field is None:
        print('CheckPhot: field name required')
        print(__doc__)
        return

    check_phot(field, np_proc=np_proc, redo=redo,
               run_mefphot=run_mefphot, outroot=outroot)


if __name__ == '__main__':
    try:
        if len(sys.argv) > 1:
            steer(sys.argv)
        else:
            print(__doc__)
    except KeyboardInterrupt:
        print('\nCheckPhot: interrupted')
