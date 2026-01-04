#!/usr/bin/env python
# coding: utf-8

"""Quantify in various ways how good star subtraction

Space Telescope Science Institute

Synopsis
--------

Quantify in various ways how good star subtraction
is between two images, in which star fluxes have been
calculated using MefPhot

Command Line Usage
------------------

::

    usage: QuantifySubtraction.py orig_file subtracted_file [more pairs...]

Description
-----------

Primary routines:

    doit - Process a single pair of files
    do_many - Process multiple pairs and create comparison table

Primary Routines
----------------

doit - Process a single pair of files
    do_many - Process multiple pairs and create comparison table

Notes
-----

History:

251229 ksl Coding begun
251229 ksl Refactored for modularity and result capture

Version History
---------------

251229 ksl
    Coding begun

251229 ksl
    Refactored for modularity and result capture

"""


import sys
from astropy.io import ascii,fits
from astropy.table import Table,join,vstack
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import binned_statistic, spearmanr

from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d


def calculate_metrics(orig, subtracted):
    """
    Calculate all subtraction quality metrics
    
    Returns a dictionary with summary statistics
    """
    # Fractional residual for each star
    fractional_residual = subtracted['Net'] / orig['Net']
    
    # Basic statistics
    metrics = {
        'median_frac_resid': np.median(fractional_residual),
        'mean_frac_resid': np.mean(fractional_residual),
        'std_frac_resid': np.std(fractional_residual),
        'rms_frac_resid': np.sqrt(np.mean(fractional_residual**2)),
        'n_stars': len(fractional_residual),
        'n_under_subtracted': int(np.sum(fractional_residual > 0.1)),
        'n_over_subtracted': int(np.sum(fractional_residual < -0.1)),
        'n_good': int(np.sum(np.abs(fractional_residual) < 0.1)),
    }
    
    # Brightness dependence
    log_flux = np.log10(orig['Net'])
    valid = np.isfinite(log_flux) & np.isfinite(fractional_residual)
    if np.sum(valid) > 10:
        corr, pval = spearmanr(log_flux[valid], fractional_residual[valid])
        metrics['brightness_correlation'] = corr
        metrics['brightness_corr_pval'] = pval
    else:
        metrics['brightness_correlation'] = np.nan
        metrics['brightness_corr_pval'] = np.nan
    
    # Color dependence (if available)
    valid_color = (~np.ma.getmaskarray(orig['B']) & 
                   ~np.ma.getmaskarray(orig['R']))
    if np.sum(valid_color) > 10:
        color = orig['B'][valid_color] - orig['R'][valid_color]
        frac_resid = fractional_residual[valid_color]
        corr, pval = spearmanr(color, frac_resid)
        metrics['color_correlation'] = corr
        metrics['color_corr_pval'] = pval
        metrics['mean_color'] = np.mean(color)
        metrics['std_color'] = np.std(color)
    else:
        metrics['color_correlation'] = np.nan
        metrics['color_corr_pval'] = np.nan
        metrics['mean_color'] = np.nan
        metrics['std_color'] = np.nan
    
    # Spatial variation (RMS across image)
    x = orig['xcenter']
    y = orig['ycenter']
    valid_spatial = np.isfinite(fractional_residual) & np.isfinite(x) & np.isfinite(y)
    if np.sum(valid_spatial) > 100:
        nbins = 10
        stat, xedges, yedges, _ = binned_statistic_2d(
            x[valid_spatial], y[valid_spatial], 
            np.abs(fractional_residual[valid_spatial]),
            statistic='std', bins=nbins
        )
        metrics['spatial_rms_range'] = np.nanmax(stat) - np.nanmin(stat)
        metrics['spatial_rms_mean'] = np.nanmean(stat)
    else:
        metrics['spatial_rms_range'] = np.nan
        metrics['spatial_rms_mean'] = np.nan
    
    return metrics


def create_combined_figure(orig, subtracted, outroot):
    """Create a single combined figure with all analysis types"""
    
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fractional_residual = subtracted['Net'] / orig['Net']
    log_flux = np.log10(orig['Net'])
    x = orig['xcenter']
    y = orig['ycenter']
    
    # === ROW 1: Brightness and color analysis ===
    
    # Fractional residual vs brightness
    ax1 = fig.add_subplot(gs[0, 0])
    valid = np.isfinite(log_flux) & np.isfinite(fractional_residual)
    ax1.scatter(log_flux[valid], fractional_residual[valid], 
           marker='.', alpha=0.05, s=1, c='steelblue', rasterized=True)
    ax1.axhline(0, color='red', linestyle='--', linewidth=1.5, zorder=5)
    
    bins = np.linspace(log_flux[valid].min(), log_flux[valid].max(), 25)
    medians, edges, _ = binned_statistic(log_flux[valid], fractional_residual[valid], 
                                      statistic='median', bins=bins)
    stds, _, _ = binned_statistic(log_flux[valid], fractional_residual[valid], 
                               statistic='std', bins=bins)
    bin_centers = (edges[:-1] + edges[1:]) / 2
    
    ax1.plot(bin_centers, medians, 'darkred', linewidth=2, zorder=10)
    ax1.fill_between(bin_centers, medians - stds, medians + stds, 
                alpha=0.3, color='red', zorder=8)
    ax1.set_ylabel('Fractional Residual', fontsize=10)
    ax1.set_xlabel('log₁₀(Net Flux)', fontsize=10)
    ax1.set_ylim(-1, 1)
    ax1.set_xlim(2, 6)
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.set_title('Subtraction vs Brightness', fontweight='bold')
    
    # Photometric uncertainty
    ax2 = fig.add_subplot(gs[0, 1])
    rel_err = orig['ErrNet'] / orig['Net']
    valid2 = np.isfinite(log_flux) & np.isfinite(rel_err)
    ax2.scatter(log_flux[valid2], rel_err[valid2], 
           marker='.', alpha=0.05, s=1, c='darkorange', rasterized=True)
    
    medians2, edges2, _ = binned_statistic(log_flux[valid2], rel_err[valid2], 
                                        statistic='median', bins=bins)
    bin_centers2 = (edges2[:-1] + edges2[1:]) / 2
    ax2.plot(bin_centers2, medians2, 'darkblue', linewidth=2, zorder=10)
    ax2.set_xlabel('log₁₀(Net Flux)', fontsize=10)
    ax2.set_ylabel('Relative Uncertainty', fontsize=10)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(2, 6)
    ax2.grid(True, alpha=0.3, linestyle=':')
    ax2.set_title('Photometric Uncertainty', fontweight='bold')
    
    # Color dependence
    ax3 = fig.add_subplot(gs[0, 2])
    valid_color = ~np.ma.getmaskarray(subtracted['B']) & ~np.ma.getmaskarray(subtracted['R'])
    
    if np.sum(valid_color) > 10:
        color = (subtracted['B'] - subtracted['R'])[valid_color]
        frac_resid = fractional_residual[valid_color]
        ax3.scatter(color, frac_resid, alpha=0.3, s=3)
        ax3.axhline(0, color='r', linestyle='--', linewidth=1.5)
        
        color_bins = np.linspace(color.min(), color.max(), 20)
        means, edges_c, _ = binned_statistic(color, frac_resid, statistic='median', bins=color_bins)
        bin_centers_c = (edges_c[:-1] + edges_c[1:]) / 2
        ax3.plot(bin_centers_c, means, 'r-', linewidth=2)
        ax3.set_ylim(-1, 1)
    else:
        ax3.text(0.5, 0.5, 'Insufficient\ncolor data', 
                ha='center', va='center', transform=ax3.transAxes)
    
    ax3.set_xlabel('B - R Color', fontsize=10)
    ax3.set_ylabel('Fractional Residual', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Subtraction vs Color', fontweight='bold')
    
    # === ROW 2: Spatial maps - Brightness before, after, and fractional residual ===
    
    valid_spatial = np.isfinite(fractional_residual) & np.isfinite(x) & np.isfinite(y)
    nbins = 50
    
    # Mean brightness map - ORIGINAL (before subtraction)
    ax4 = fig.add_subplot(gs[1, 0])
    valid_flux_orig = np.isfinite(log_flux) & np.isfinite(x) & np.isfinite(y)
    mean_flux_orig, xedges_orig, yedges_orig, _ = binned_statistic_2d(
        x[valid_flux_orig], y[valid_flux_orig], log_flux[valid_flux_orig],
        statistic='mean', bins=nbins
    )
    im1 = ax4.imshow(mean_flux_orig.T, origin='lower', cmap='viridis',
                extent=[xedges_orig[0], xedges_orig[-1], yedges_orig[0], yedges_orig[-1]],
                aspect='auto', interpolation='bilinear')
    plt.colorbar(im1, ax=ax4, label='log₁₀ Flux (orig)')
    ax4.set_xlabel('X (pixels)', fontsize=10)
    ax4.set_ylabel('Y (pixels)', fontsize=10)
    ax4.set_title('Spatial: Brightness (Original)', fontweight='bold')
    
    # Mean brightness map - SUBTRACTED (residual flux)
    ax5 = fig.add_subplot(gs[1, 1])
    log_flux_sub = np.log10(np.abs(subtracted['Net']) + 1)  # Add 1 to handle zeros/negatives
    valid_flux_sub = np.isfinite(log_flux_sub) & np.isfinite(x) & np.isfinite(y)
    mean_flux_sub, xedges_sub, yedges_sub, _ = binned_statistic_2d(
        x[valid_flux_sub], y[valid_flux_sub], log_flux_sub[valid_flux_sub],
        statistic='mean', bins=nbins
    )
    im2 = ax5.imshow(mean_flux_sub.T, origin='lower', cmap='viridis',
                extent=[xedges_sub[0], xedges_sub[-1], yedges_sub[0], yedges_sub[-1]],
                aspect='auto', interpolation='bilinear')
    plt.colorbar(im2, ax=ax5, label='log₁₀ |Flux| (sub)')
    ax5.set_xlabel('X (pixels)', fontsize=10)
    ax5.set_ylabel('Y (pixels)', fontsize=10)
    ax5.set_title('Spatial: Brightness (Residual)', fontweight='bold')
    
    # Fractional residual map - use binned median for better structure visibility
    ax6 = fig.add_subplot(gs[1, 2])
    
    # Use binned statistic instead of interpolation
    frac_resid_binned, xedges_fr, yedges_fr, _ = binned_statistic_2d(
        x[valid_spatial], y[valid_spatial], fractional_residual[valid_spatial],
        statistic='median', bins=nbins
    )
    
    # Auto-scale to data range, using robust percentiles to avoid outliers
    valid_bins = np.isfinite(frac_resid_binned)
    if np.sum(valid_bins) > 0:
        vmin_data = np.nanpercentile(frac_resid_binned[valid_bins], 5)
        vmax_data = np.nanpercentile(frac_resid_binned[valid_bins], 95)
        vmax_abs = max(abs(vmin_data), abs(vmax_data))
        
        # Use symmetric scale around zero, but adjusted to data
        norm = TwoSlopeNorm(vmin=-vmax_abs, vcenter=0, vmax=vmax_abs)
    else:
        norm = TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=0.1)
    
    im3 = ax6.imshow(frac_resid_binned.T, origin='lower', cmap='RdBu_r', norm=norm,
                     extent=[xedges_fr[0], xedges_fr[-1], yedges_fr[0], yedges_fr[-1]],
                     aspect='auto', interpolation='bilinear')
    plt.colorbar(im3, ax=ax6, label='Median Frac Resid')
    ax6.set_xlabel('X (pixels)', fontsize=10)
    ax6.set_ylabel('Y (pixels)', fontsize=10)
    ax6.set_title('Spatial: Fractional Residual', fontweight='bold')
    
    # === ROW 3: Summary statistics and color map ===
    
    # Text summary
    ax7 = fig.add_subplot(gs[2, :2])
    ax7.axis('off')
    
    metrics = calculate_metrics(orig, subtracted)
    summary_text = f"""
SUMMARY STATISTICS

Total Stars: {metrics['n_stars']}
Good Subtraction (|resid| < 0.1): {metrics['n_good']} ({100*metrics['n_good']/metrics['n_stars']:.1f}%)
Under-subtracted (resid > 0.1): {metrics['n_under_subtracted']} ({100*metrics['n_under_subtracted']/metrics['n_stars']:.1f}%)
Over-subtracted (resid < -0.1): {metrics['n_over_subtracted']} ({100*metrics['n_over_subtracted']/metrics['n_stars']:.1f}%)

Median Fractional Residual: {metrics['median_frac_resid']:.4f}
RMS Fractional Residual: {metrics['rms_frac_resid']:.4f}

Brightness Correlation: ρ = {metrics['brightness_correlation']:.3f} (p = {metrics['brightness_corr_pval']:.2e})
Color Correlation: ρ = {metrics['color_correlation']:.3f} (p = {metrics['color_corr_pval']:.2e})

Spatial RMS Range: {metrics['spatial_rms_range']:.4f}
    """
    
    ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes, 
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # Mean color map
    ax8 = fig.add_subplot(gs[2, 2])
    valid_color_spatial = (~np.ma.getmaskarray(orig['B']) & 
                          ~np.ma.getmaskarray(orig['R']) &
                          np.isfinite(x) & np.isfinite(y))
    
    if np.sum(valid_color_spatial) > 100:
        color_spatial = orig['B'][valid_color_spatial] - orig['R'][valid_color_spatial]
        mean_color, xedges, yedges, _ = binned_statistic_2d(
            x[valid_color_spatial], y[valid_color_spatial], color_spatial,
            statistic='mean', bins=nbins
        )
        im4 = ax8.imshow(mean_color.T, origin='lower', cmap='coolwarm',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                    aspect='auto', interpolation='bilinear')
        plt.colorbar(im4, ax=ax8, label='B-R')
        ax8.set_title('Spatial: Mean Color', fontweight='bold')
    else:
        ax8.text(0.5, 0.5, 'Insufficient\ncolor data', 
                ha='center', va='center', transform=ax8.transAxes)
        ax8.set_title('Spatial: Mean Color', fontweight='bold')
    
    ax8.set_xlabel('X (pixels)', fontsize=10)
    ax8.set_ylabel('Y (pixels)', fontsize=10)
    
    fig.suptitle(f'Star Subtraction Quality Analysis: {outroot}', 
                fontsize=14, fontweight='bold', y=0.995)
    
    return fig


def doit(orig_phot, subtracted_phot, outroot='', combined=True):
    """
    Analyze subtraction quality for a single pair of files
    
    Parameters:
    -----------
    orig_phot : str
        Path to original photometry table
    subtracted_phot : str
        Path to subtracted photometry table
    outroot : str
        Root name for output files (derived from subtracted_phot if empty)
    combined : bool
        If True, create a single combined figure
    
    Returns:
    --------
    metrics : dict
        Dictionary of calculated metrics
    """
    
    orig = Table.read(orig_phot)
    subtracted = Table.read(subtracted_phot)
    
    # Calculate all metrics
    metrics = calculate_metrics(orig, subtracted)
    
    # Add file information
    metrics['orig_file'] = orig_phot
    metrics['subtracted_file'] = subtracted_phot
    
    # Determine output root
    if outroot == '':
        outroot = subtracted_phot.split('/')[-1].replace('.fits', '')
    metrics['outroot'] = outroot
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Analysis: {outroot}")
    print(f"{'='*60}")
    print(f"Median fractional residual: {metrics['median_frac_resid']:.4f}")
    print(f"RMS fractional residual: {metrics['rms_frac_resid']:.4f}")
    print(f"Under-subtracted stars: {metrics['n_under_subtracted']}")
    print(f"Over-subtracted stars: {metrics['n_over_subtracted']}")
    print(f"Good subtraction: {metrics['n_good']}")
    
    # Create output directory
    os.makedirs('SubQual', exist_ok=True)
    
    # Generate combined plot
    if combined:
        fig = create_combined_figure(orig, subtracted, outroot)
        fig.savefig(f'SubQual/{outroot}.combined.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: SubQual/{outroot}.combined.png")
    
    return metrics


def do_many(file_pairs, output_table='SubQual/comparison_metrics.fits'):
    """
    Process multiple file pairs and create comparison table
    
    Parameters:
    -----------
    file_pairs : list of tuples
        List of (orig_phot, subtracted_phot) pairs
    output_table : str
        Path to save the comparison metrics table
    
    Returns:
    --------
    results : astropy.table.Table
        Table with metrics for all file pairs
    """
    
    all_metrics = []
    
    for i, (orig_phot, subtracted_phot) in enumerate(file_pairs):
        print(f"\nProcessing pair {i+1}/{len(file_pairs)}")
        try:
            metrics = doit(orig_phot, subtracted_phot, combined=True)
            all_metrics.append(metrics)
        except Exception as e:
            print(f"Error processing {orig_phot}, {subtracted_phot}: {e}")
            continue
    
    # Convert to astropy table
    if len(all_metrics) > 0:
        results = Table(all_metrics)
        
        # Save results
        os.makedirs(os.path.dirname(output_table), exist_ok=True)
        results.write(output_table, overwrite=True)
        print(f"\nSaved comparison table: {output_table}")
        
        # Print comparison summary
        print(f"\n{'='*80}")
        print("COMPARISON SUMMARY")
        print(f"{'='*80}")
        print(f"{'Outroot':<30} {'RMS Resid':>10} {'Color Corr':>11} {'N Good':>8} {'N Under':>8} {'N Over':>8}")
        print(f"{'-'*80}")
        for row in results:
            print(f"{row['outroot']:<30} {row['rms_frac_resid']:10.4f} {row['color_correlation']:11.3f} "
                  f"{row['n_good']:8d} {row['n_under_subtracted']:8d} {row['n_over_subtracted']:8d}")
        
        return results
    else:
        print("No successful analyses to compare")
        return None

from glob import glob
from astropy.table import Table,join,vstack

def get_pairs(xdir='TabPhot',prefix='LMC'):
    search_string='%s/%s*fits' % (xdir,prefix)
    all_files=glob(search_string)
    print(len(all_files))
    xtype=[]
    xroot=[]
    for one in all_files:
        one_type=one.split('.')[-2]
        one_root=one.replace('%s.fits' % one_type,'')
        xtype.append(one_type)
        xroot.append(one_root)
    xtab=Table([xroot,xtype,all_files],names=['Root','Type','Filename'])

    ha_sub=xtab[xtab['Type']=='ha_sub_r']
    s2_sub=xtab[xtab['Type']=='s2_sub_r']
    ha=xtab[xtab['Type']=='N662']
    s2=xtab[xtab['Type']=='N673']

    qha=join(ha_sub,ha,keys=['Root'])
    qs2=join(s2_sub,s2,keys=['Root'])
    xfinal=vstack([qha,qs2])
    xfinal.rename_columns(['Filename_1','Filename_2'],['Subtracted','Original'])
    xfinal.rename_columns(['Type_1','Type_2'],['Sub_type','Orig_Type'])
    return xfinal
        

def steer(argv):
    '''
    Steering routine for command line usage
    '''

    xdir='TabPhot'
    xprefix=''
    
    if len(argv) < 2:
        print(__doc__)
        return
    
    files = []
    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i]=='-dir':
            i+=1
            xdir=argv[i]
        elif argv[i]=='-pre':
            i+=1
            xprefix=argv[i]

        elif argv[i][0] == '-':
            print('Error: Could not interpret commands:', argv)
            return
        elif argv[i].count('.fits'):
            files.append(argv[i])
        else:
            print('Error: cannot parse command line:', argv)
            return
        i += 1
    
    # Group files into pairs (every two consecutive files)
    if xprefix!='':
        xpairs=get_pairs(xdir,xprefix)
        print(xpairs)
        file_pairs=[(xpairs['Original'][i],xpairs['Subtracted'][i]) for i in range(0,len(xpairs))]
        print(file_pairs)
        do_many(file_pairs)
    elif len(files) == 2:
        # Single pair
        doit(files[0], files[1])
    elif len(files) % 2 == 0 and len(files) > 2:
        # Multiple pairs
        file_pairs = [(files[i], files[i+1]) for i in range(0, len(files), 2)]
        do_many(file_pairs)
    else:
        print(f"Error: Expected pairs of files, got {len(files)} files")
        return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)        
    else:
        print(__doc__)
