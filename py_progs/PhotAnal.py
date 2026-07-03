#!/usr/bin/env python
# coding: utf-8

"""Compare stellar photometry between two filters to assess inter-filter flux consistency

Space Telescope Science Institute

Synopsis
--------

PhotAnal measures the same stars in two different filters and reports how
consistently they are scaled relative to each other.  This is the key
diagnostic for CleanStars continuum subtraction: if a star has DN_r in the
r-band image and DN_N662 in the Hα image, CleanStars computes

    ha_clean = N662 - r_pure

and a star cancels only if DN_N662 = DN_r.  After MefPrep has placed all
images on the mag-28 scale, this condition is equivalent to asking whether
the star has the same magc (= phot_mag + MAGZERO − 28) in both filters.
PhotAnal measures that difference and its scatter.

Command Line Usage
------------------

::

    PhotAnal.py -filter1 r -filter2 N662 TabPhot/*.gaia.fits
    PhotAnal.py -filter1 r -filter2 N708 TabPhot/*.gaia.fits
    PhotAnal.py -filter1 r -filter2 N662 -o result.fits -no_plot TabPhot/*.gaia.fits

Flags
-----

-filter1 NAME   Reference filter (e.g. r)
-filter2 NAME   Comparison filter (e.g. N662, N673, N708)
-o FILE         Output FITS table (default: filter_compare_<f1>_<f2>.fits)
-no_plot        Skip the diagnostic plot
-h              Print this help

Description
-----------

Both filter sets must be Gaia-matched TabPhot files (.gaia.fits).  Run::

    MefPhot.py -r 6 DECam_MEF/LMC_c42/*.fits

for all filters before running PhotAnal.

For each star detected in both filter1 and filter2:

    delta_mag = magc_filter1 − magc_filter2
              = (phot_mag_1 + MAGZERO_1 − 28) − (phot_mag_2 + MAGZERO_2 − 28)

This is the quantity that determines the stellar residual after continuum
subtraction.  A non-zero mean indicates a systematic inter-filter offset
(correctable by adjusting MAGZERO for one filter).  The scatter around
the best-fit color correction is the irreducible floor:

    typical stellar residual ≈ σ_residual × ln(10) / 2.5
                              ≈ σ_residual × 0.92    (in fractional DN units)

Output
------

Printed summary: N_stars, mean Δmag, total scatter, color term b, residual
scatter, and typical stellar residual percentage.

filter_compare_<f1>_<f2>.fits: per-star table with Source_name, magc in
each filter, delta_mag, Gaia G−R color, delta_mag after color correction.

filter_compare_<f1>_<f2>.png: Δmag vs G−R color with best-fit line.

Notes
-----

The comparison uses per-source median magc aggregated across all files,
so it is not sensitive to intra-filter frame-to-frame scatter (that is
diagnosed separately by PhotEval).

History
-------

251123 ksl  Prototype coded
260701 ksl  Rewritten as inter-filter consistency tool
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, join
from astropy.io import fits
from astropy.stats import sigma_clipped_stats


def load_filter_data(filenames, filter_name, extra_zp_lookups=None):
    """
    Load all .gaia.fits TabPhot files matching filter_name.

    Returns (aggregate, chunks) where aggregate is an Astropy Table with one
    row per star: Source_name, magc (median across files using header MAGZERO),
    G, R (Gaia magnitudes), and optionally magc_{name} for each entry in
    extra_zp_lookups.  chunks is the list of per-file tables (before
    aggregation) needed for per-image-pair statistics.

    Parameters
    ----------
    extra_zp_lookups : dict of {name: {basename: zp_calc}}, optional
        Additional ZP sources (e.g. {'gaia_g': {...}, 'smash_r': {...}}).
        Files not in a lookup get NaN for that ZP source.
    """
    extra_safes = []
    if extra_zp_lookups:
        extra_safes = [n.replace('-', '_') for n in extra_zp_lookups]

    chunks = []
    n_files = 0

    for fname in filenames:
        basename = os.path.basename(fname)
        try:
            with fits.open(fname) as hdul:
                hdr = dict(hdul[1].header)
            tab = Table.read(fname)
        except Exception as e:
            print(f'PhotAnal: cannot read {basename}: {e}')
            continue

        # Filter name: header first, column fallback
        filt = hdr.get('FILTER', hdr.get('filter', None))
        if filt is None and 'Filter' in tab.colnames:
            filt = str(tab['Filter'][0]).strip()
        if filt is None:
            continue
        filt = str(filt).split()[0]
        if filt != filter_name:
            continue

        # MAGZERO: header first, column fallback
        magzero = hdr.get('MAGZERO', None)
        if magzero is None and 'MAGZERO' in tab.colnames:
            magzero = float(np.nanmedian(tab['MAGZERO']))
        if magzero is None:
            print(f'PhotAnal: no MAGZERO in {basename}, skipping')
            continue

        # Quality: drop saturated pixels and non-finite phot_mag
        mask = np.ones(len(tab), dtype=bool)
        if 'Max' in tab.colnames:
            mask &= np.array(tab['Max'], dtype=float) < 45000
        if 'phot_mag' in tab.colnames:
            mask &= np.isfinite(np.array(tab['phot_mag'], dtype=float))
        tab = tab[mask]
        if len(tab) == 0:
            continue

        phot = np.array(tab['phot_mag'], dtype=float)
        tab['magc'] = phot + float(magzero) - 28.0

        # Extra ZP sources — NaN for files not in the lookup
        if extra_zp_lookups:
            for (name, lkup), safe in zip(extra_zp_lookups.items(), extra_safes):
                zp_emp = float(lkup.get(basename, np.nan))
                tab[f'magc_{safe}'] = phot + zp_emp - 28.0 if np.isfinite(zp_emp) else np.nan

        keep = ['Source_name', 'magc'] + [f'magc_{s}' for s in extra_safes
                                           if f'magc_{s}' in tab.colnames]
        for col in ('G', 'R'):
            if col in tab.colnames:
                keep.append(col)
        chunks.append(tab[keep])
        n_files += 1

    if not chunks:
        print(f'PhotAnal: no usable files found for filter "{filter_name}"')
        return None, []

    print(f'PhotAnal: {n_files} files loaded for filter {filter_name}')

    # Per-source median magc using pandas groupby (much faster than astropy group_by
    # on million-row tables)
    import pandas as pd

    # Build a flat pandas DataFrame from the per-file chunks
    frames = []
    for c in chunks:
        d = {'Source_name': [str(x) for x in c['Source_name']],
             'magc':        np.array(c['magc'], dtype=float)}
        for col in ('G', 'R'):
            if col in c.colnames:
                d[col] = np.array(c[col], dtype=float)
        for s in extra_safes:
            col = f'magc_{s}'
            if col in c.colnames:
                d[col] = np.array(c[col], dtype=float)
        frames.append(pd.DataFrame(d))

    df = pd.concat(frames, ignore_index=True)

    agg_dict = {'magc': 'median'}
    for col in ('G', 'R'):
        if col in df.columns:
            agg_dict[col] = 'median'
    for s in extra_safes:
        col = f'magc_{s}'
        if col in df.columns:
            agg_dict[col] = 'median'

    agg = df.groupby('Source_name', sort=False).agg(agg_dict).reset_index()
    print(f'PhotAnal: {len(agg):,} unique sources for filter {filter_name}')

    return Table.from_pandas(agg), chunks


def _compute_pair_stats(chunks1, chunks2, min_stars=20, mag_lim=19.0,
                        magc_col='magc'):
    """Compute the scatter of per-image-pair mean Δmag.

    For each (filter1 file, filter2 file) pair, matches stars by Source_name,
    computes the mean Δmag for that pair, then returns the MAD scatter of
    those per-pair means.  Only stars with R <= mag_lim are used so that
    the pair mean is driven by high-S/N detections.

    Parameters
    ----------
    magc_col : str
        Column name to use for calibrated magnitude (default 'magc' = MAGZERO;
        use e.g. 'magc_gaia_g' to compute pair scatter with Gaia G ZP).

    Returns (N_pairs, mean_pair_delta_mag, std_pair_delta_mag).
    """
    from astropy.stats import mad_std as _mad_std

    def _to_dict(tab):
        col = magc_col if magc_col in tab.colnames else 'magc'
        names = np.array([str(n) for n in tab['Source_name']])
        magcs = np.array(tab[col], dtype=float)
        ok = np.isfinite(magcs)
        if 'R' in tab.colnames:
            r = np.array(tab['R'], dtype=float)
            ok &= np.isfinite(r) & (r >= 15.0) & (r <= mag_lim)
        return dict(zip(names[ok], magcs[ok]))

    dicts1 = [_to_dict(t) for t in chunks1]
    dicts2 = [_to_dict(t) for t in chunks2]
    sets2  = [set(d) for d in dicts2]

    pair_means = []
    for d1 in dicts1:
        if not d1:
            continue
        s1 = set(d1)
        for d2, s2 in zip(dicts2, sets2):
            common = s1 & s2
            if len(common) < min_stars:
                continue
            delta = np.array([d1[n] - d2[n] for n in common])
            ok = np.isfinite(delta)
            if ok.sum() < min_stars:
                continue
            pair_means.append(float(np.mean(delta[ok])))

    if len(pair_means) < 2:
        return 0, np.nan, np.nan
    arr = np.array(pair_means)
    return len(arr), float(np.mean(arr)), float(_mad_std(arr))


def do_compare(filenames, filter1, filter2, outfile=None, do_plot=True,
               extra_zp_lookups=None, _preloaded=None):
    """
    Cross-match filter1 and filter2 photometry and report inter-filter scatter.

    Parameters
    ----------
    filenames : list of str
        Paths to .gaia.fits TabPhot files (both filters mixed together).
    filter1 : str
        Reference filter name (e.g. 'r').
    filter2 : str
        Comparison filter name (e.g. 'N662').
    outfile : str, optional
        Output FITS table path.  Default: filter_compare_<f1>_<f2>.fits.
    do_plot : bool
        If True, save a PNG diagnostic plot.
    extra_zp_lookups : dict of {name: {basename: zp_calc}}, optional
        Additional ZP sources to evaluate in parallel (e.g. from CheckPhot's
        CalcZeroPoint runs).  For each name the same sigma metrics are
        computed and stored in the output meta as STD_DM_{name},
        STD_RE_{name}, N_PR_{name}, STD_PD_{name}.
    _preloaded : tuple of two (Table, list) pairs, optional
        Pre-loaded (tab, chunks) from load_filter_data() for (filter1, filter2).
        When supplied, load_filter_data() is skipped (avoids redundant I/O when
        a filter appears in multiple pairs).
    """
    if _preloaded is not None:
        # Copy so in-place renames below don't corrupt the caller's cache
        (tab1_raw, chunks1), (tab2_raw, chunks2) = _preloaded
        tab1, tab2 = tab1_raw.copy(), tab2_raw.copy()
    else:
        print(f'\nPhotAnal: loading {filter1} data...')
        tab1, chunks1 = load_filter_data(filenames, filter1,
                                         extra_zp_lookups=extra_zp_lookups)
        print(f'PhotAnal: loading {filter2} data...')
        tab2, chunks2 = load_filter_data(filenames, filter2,
                                         extra_zp_lookups=extra_zp_lookups)

    if tab1 is None or tab2 is None:
        return None

    # Rename before joining to avoid column-name clashes
    tab1.rename_column('magc', 'magc_1')
    tab2.rename_column('magc', 'magc_2')
    for col in ('G', 'R'):
        if col in tab1.colnames:
            tab1.rename_column(col, f'{col}_1')
        if col in tab2.colnames:
            tab2.rename_column(col, f'{col}_2')

    matched = join(tab1, tab2, keys='Source_name', join_type='inner')
    print(f'PhotAnal: {len(tab1)} {filter1} stars, '
          f'{len(tab2)} {filter2} stars, {len(matched)} matched')

    if len(matched) < 10:
        print('PhotAnal: too few matched stars for meaningful statistics')
        return None

    matched['delta_mag'] = np.array(matched['magc_1'], dtype=float) - \
                           np.array(matched['magc_2'], dtype=float)

    # General quality cut: reject wildly discrepant Δmag
    finite = np.isfinite(matched['delta_mag']) & (np.abs(matched['delta_mag']) < 5.0)
    if 'R_1' in matched.colnames:
        r_arr = np.array(matched['R_1'], dtype=float)
        finite &= np.isfinite(r_arr) & (r_arr > 15.0) & (r_arr < 21.5)
    matched = matched[finite]

    if len(matched) < 10:
        print('PhotAnal: too few stars after quality cuts')
        return None

    delta = np.array(matched['delta_mag'], dtype=float)
    _, mean_d, std_d = sigma_clipped_stats(delta, sigma=3.0)

    # Color correction: fit delta_mag = a + b*(G−R)
    has_color = ('G_1' in matched.colnames and 'R_1' in matched.colnames)
    color = None
    a, b = mean_d, 0.0
    std_res = std_d
    if has_color:
        color = np.array(matched['G_1'], dtype=float) - \
                np.array(matched['R_1'], dtype=float)
        ok = np.isfinite(color) & np.isfinite(delta)
        if np.sum(ok) > 20:
            coeffs = np.polyfit(color[ok], delta[ok], 1)
            b, a = float(coeffs[0]), float(coeffs[1])
            residual = delta[ok] - (a + b * color[ok])
            _, _, std_res = sigma_clipped_stats(residual, sigma=3.0)
        else:
            has_color = False

    # Fractional stellar residual from σ (linearised: 10^(σ/2.5) − 1 ≈ 0.92σ)
    def _frac_pct(sigma):
        return (10 ** (abs(sigma) / 2.5) - 1.0) * 100.0

    # Per-image-pair scatter: how consistently are the two filters matched image-to-image?
    print(f'PhotAnal: computing per-image-pair statistics...')
    n_pairs, mean_pair_dm, std_pair_dm = _compute_pair_stats(chunks1, chunks2)

    # ---- Print summary ----
    print()
    print(f'--- Inter-filter consistency: {filter1} vs {filter2} ---')
    print(f'  N stars matched          : {len(matched)}')
    print(f'  Mean Δmag ({filter1}−{filter2})'.ljust(28) +
          f': {mean_d:+.4f} mag')
    print(f'  Scatter (total)          : σ = {std_d:.4f} mag  '
          f'→ {_frac_pct(std_d):.1f}% typical stellar residual')
    if has_color:
        print(f'  Color term b (G−R)       : {b:+.4f} mag/mag')
        print(f'  Scatter (after color)    : σ = {std_res:.4f} mag  '
              f'→ {_frac_pct(std_res):.1f}% typical stellar residual  [floor]')
    if n_pairs >= 2:
        print(f'  Image-pair scatter       : N_pairs={n_pairs}  '
              f'std_pair_dmag={std_pair_dm:.4f} mag  '
              f'(image-to-image ZP variation between filters)')
    print()

    # ---- Output table ----
    if outfile is None:
        outfile = f'filter_compare_{filter1}_{filter2}.fits'

    out = Table()
    out['Source_name'] = matched['Source_name']
    out['magc_1'] = np.array(matched['magc_1'], dtype=float)
    out['magc_2'] = np.array(matched['magc_2'], dtype=float)
    out['delta_mag'] = delta
    if has_color and color is not None:
        out['G_R'] = color
        ok_full = np.isfinite(color) & np.isfinite(delta)
        corr = np.full(len(matched), np.nan)
        corr[ok_full] = delta[ok_full] - (a + b * color[ok_full])
        out['delta_mag_corrected'] = corr
    if 'R_1' in matched.colnames:
        out['R_Gaia'] = np.array(matched['R_1'], dtype=float)
    out.meta['FILTER1']   = filter1
    out.meta['FILTER2']   = filter2
    out.meta['MEAN_DM']   = round(float(mean_d),       5)
    out.meta['STD_DM']    = round(float(std_d),        5)
    out.meta['COLOR_B']   = round(float(b),            5)
    out.meta['STD_RES']   = round(float(std_res),      5)
    out.meta['N_PAIRS']   = int(n_pairs)
    out.meta['MEAN_PDM']  = round(float(mean_pair_dm), 5) if np.isfinite(mean_pair_dm) else np.nan
    out.meta['STD_PDM']   = round(float(std_pair_dm),  5) if np.isfinite(std_pair_dm)  else np.nan

    # ---- Extra ZP sources: compute the same sigma metrics ----
    if extra_zp_lookups:
        # Collect G−R colors (same for all ZP sources)
        color_ok = None
        if has_color:
            color_ok = np.isfinite(color) & np.isfinite(delta)

        for name, _ in extra_zp_lookups.items():
            safe = name.replace('-', '_')
            c1 = f'magc_{safe}_1'
            c2 = f'magc_{safe}_2'
            if c1 not in matched.colnames or c2 not in matched.colnames:
                continue
            de = np.array(matched[c1], dtype=float) - np.array(matched[c2], dtype=float)
            ok_e = np.isfinite(de) & (np.abs(de) < 5.0)
            if 'R_1' in matched.colnames:
                ok_e &= np.isfinite(np.array(matched['R_1'], dtype=float))
            if ok_e.sum() < 10:
                continue

            _, mean_e, std_e = sigma_clipped_stats(de[ok_e], sigma=3.0)

            # Color correction using same G−R colors
            std_res_e = std_e
            if has_color and color_ok is not None:
                ok_ce = ok_e & color_ok
                if ok_ce.sum() > 20:
                    cf = np.polyfit(color[ok_ce], de[ok_ce], 1)
                    res_e = de[ok_ce] - (cf[1] + cf[0] * color[ok_ce])
                    _, _, std_res_e = sigma_clipped_stats(res_e, sigma=3.0)

            np_e, _, std_pdm_e = _compute_pair_stats(chunks1, chunks2,
                                                      magc_col=f'magc_{safe}')
            print(f'  [{name}]  sigma_total={std_e:.4f}  '
                  f'sigma_residual={std_res_e:.4f}  '
                  f'std_pair_dmag={std_pdm_e:.4f}')

            # FITS header keys: 8-char limit; encode first 2 chars of safe name
            tag = safe[:2].upper()  # 'GA' for gaia_g, 'SM' for smash_r
            out.meta[f'SDTM_{tag}'] = round(float(std_e),     5)
            out.meta[f'STRE_{tag}'] = round(float(std_res_e), 5)
            out.meta[f'NPR_{tag}']  = int(np_e)
            out.meta[f'SPDT_{tag}'] = round(float(std_pdm_e), 5) \
                                      if np.isfinite(std_pdm_e) else np.nan

    out.write(outfile, overwrite=True)
    print(f'PhotAnal: wrote {outfile}')

    # ---- Plot ----
    if do_plot:
        figname = outfile.replace('.fits', '.png')
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        ax = axes[0]
        if has_color and color is not None:
            ok = np.isfinite(color) & np.isfinite(delta)
            ax.scatter(color[ok], delta[ok], s=2, alpha=0.2, color='steelblue',
                       rasterized=True)
            xlim = np.nanpercentile(color[ok], [1, 99])
            xfit = np.linspace(xlim[0], xlim[1], 100)
            ax.plot(xfit, a + b * xfit, 'r-', lw=1.5,
                    label=f'Δmag = {a:+.3f} {b:+.3f}·(G−R)')
            ax.axhline(mean_d, color='gray', ls='--', lw=1,
                       label=f'mean = {mean_d:+.3f}')
            ax.set_xlabel('Gaia G − R (mag)')
            ax.legend(fontsize=8)
        else:
            ax.hist(delta[np.isfinite(delta)], bins=60, color='steelblue',
                    alpha=0.7)
            ax.set_xlabel(f'Δmag ({filter1}−{filter2})')
        ax.set_ylabel(f'Δmag ({filter1}−{filter2})')
        ax.set_title(f'{filter1} − {filter2}  σ={std_d:.3f} mag')

        ax = axes[1]
        if has_color and color is not None:
            ok = np.isfinite(color) & np.isfinite(delta)
            residual_all = delta - (a + b * color)
            ax.hist(residual_all[ok & np.isfinite(residual_all)], bins=60,
                    color='steelblue', alpha=0.7)
            ax.axvline(0, color='r', lw=1)
            ax.set_xlabel(f'Δmag − color fit  (σ = {std_res:.3f} mag)')
            ax.set_title(f'Residual after color correction  '
                         f'→ {_frac_pct(std_res):.1f}% stellar residual [floor]')
        else:
            ax.hist(delta[np.isfinite(delta)] - mean_d, bins=60,
                    color='steelblue', alpha=0.7)
            ax.axvline(0, color='r', lw=1)
            ax.set_xlabel(f'Δmag − mean  (σ = {std_d:.3f} mag)')
            ax.set_title(f'Scatter  → {_frac_pct(std_d):.1f}% stellar residual')
        ax.set_ylabel('N stars')

        plt.suptitle(f'Inter-filter consistency: {filter1} vs {filter2}  '
                     f'(N={len(matched)})', fontsize=11)
        plt.tight_layout()
        plt.savefig(figname, dpi=150)
        plt.close()
        print(f'PhotAnal: wrote {figname}')

    return out


def steer(argv):
    filenames = []
    filter1 = None
    filter2 = None
    outfile = None
    do_plot = True

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(__doc__)
            return
        elif arg == '-filter1':
            i += 1
            filter1 = argv[i]
        elif arg == '-filter2':
            i += 1
            filter2 = argv[i]
        elif arg == '-o':
            i += 1
            outfile = argv[i]
        elif arg == '-no_plot':
            do_plot = False
        elif arg.startswith('-'):
            print(f'PhotAnal: unknown flag {arg}')
            print(__doc__)
            return
        else:
            filenames.append(arg)
        i += 1

    if filter1 is None or filter2 is None:
        print('PhotAnal: -filter1 and -filter2 are required')
        print(__doc__)
        return
    if not filenames:
        print('PhotAnal: no TabPhot files specified')
        print(__doc__)
        return

    do_compare(filenames, filter1, filter2, outfile=outfile, do_plot=do_plot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
