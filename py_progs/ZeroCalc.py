#!/usr/bin/env python
# coding: utf-8

"""Calculate the magnitude zero point for an image

Space Telescope Science Institute

Synopsis
--------

Given one or more forced-photometry tables (output of MefPhot), fit a
linear model to determine the magnitude zero point and color term that
place instrumental magnitudes on the Gaia or SMASH photometric scale.
The derived zero point can be compared directly to the MAGZERO keyword
carried in the MEF file header.

Command Line Usage
------------------

::

    ZeroCalc.py [-h] [-R] [-G] [-color] [-smash] [-out ROOT] file1.fits file2.fits ...

    -h        Print this help and exit
    -R        Fit to the reference catalog R band (default)
    -G        Fit to the reference catalog G band
    -color    Add a color term to the fit (default: simple zero-point fit only)
    -smash    Input tables were produced with SMASH as the reference catalog
    -out ROOT Output table root name (default: MagZero)

Description
-----------

Two fit modes are available:

**Simple fit** (default)::

    m_ref = m_inst + c_0

**Color-corrected fit** (``-color``)::

    m_ref = m_inst + c_0 + c_1 * color

where the color predictor is chosen to be independent of the target band:

+--------+---------+-------+
| Target | Catalog | Color |
+========+=========+=======+
| R      | Gaia    | G−R   |
+--------+---------+-------+
| G      | Gaia    | B−R   |
+--------+---------+-------+
| R      | SMASH   | G−R   |
+--------+---------+-------+
| G      | SMASH   | U−R   |
+--------+---------+-------+

In all cases ``m_inst = 28 - 2.5 * log10(flux)`` and the fitted ``c_0``
gives the zero-point correction: ``ZP_derived = 28 + c_0``.

Output filenames encode the band, catalog, and fit mode::

    MagZero.<band>.gaia.txt          simple fit, Gaia (default)
    MagZero.<band>.gaia.color.txt    color-corrected fit, Gaia
    MagZero.<band>.smash.txt         simple fit, SMASH
    MagZero.<band>.smash.color.txt   color-corrected fit, SMASH

Each row in the summary table contains: Filter, Exptime, Root, MagZero
(= 28 + c_0), c_0, c_1, rms, HdrZero (pipeline MAGZERO from the MEF
header), Catalog, and Filename.

Diagnostic plots are written to ``FigZero/`` with matching suffixes::

    FigZero/<band>_<root>.gaia.png
    FigZero/<band>_<root>.gaia.color.png
    FigZero/<band>_<root>.smash.png
    FigZero/<band>_<root>.smash.color.png

Primary Routines
----------------

do_one
    Process a single photometry table and return fit results.

do_many
    Process multiple tables and accumulate results into a summary file.

Notes
-----

This routine was developed to assess the consistency of MAGZERO as
delivered by the DECam community pipeline.  Comparing the HdrZero column
(pipeline MAGZERO) with the MagZero column (28 + c_0) across many
exposures reveals systematic trends with filter, time, or CCD.  Running
with both Gaia and SMASH provides an additional cross-check because the
SMASH color term should be close to zero for r-band data.

Version History
---------------

251128 ksl
    Coding begun

251228 ksl
    Updated to allow fitting to the Gaia G band

260414 ksl
    Added SMASH support (-smash flag).
    Output filenames now include catalog suffix (.gaia / .smash).
    Plot axis labels and titles reflect the reference catalog used.

260504 ksl
    Add simple fit mode (default) and optional color-corrected fit (-color).
    Color predictor is now always independent of the target band:
    Gaia R: G-R, Gaia G: B-R, SMASH R: G-R, SMASH G: U-R.
    Output filenames include .color suffix when -color is used.

"""


import sys
import matplotlib.pyplot as plt
import os
from glob import glob
from astropy.table import Table, join,vstack
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import Table

def fit_magnitude_model(data, use_color=False):
    """
    Fit catalogued magnitudes to either a simple or color-corrected model.

    Simple (use_color=False):   m_ref = m_inst + c_0
    Color-corrected (use_color=True): m_ref = m_inst + c_0 + c_1 * color

    Parameters
    ----------
    data : astropy.table.Table
        Table with columns 'phot_mag', 'target_mag', and (if use_color)
        'target_color'.
    use_color : bool
        If True, include a color term in the fit. Default False.

    Returns
    -------
    dict with keys c_0, c_0_err, c_1, c_1_err, rms, table.
    c_1 and c_1_err are 0.0 when use_color=False.
    """
    mag_obs   = np.asarray(data['phot_mag'], dtype=float)
    target    = np.asarray(data['target_mag'], dtype=float)

    if use_color:
        color = np.asarray(data['target_color'], dtype=float)

        def model_wrapper(x_data, c_0, c_1):
            return x_data[0] + c_0 + c_1 * x_data[1]

        popt, pcov = curve_fit(model_wrapper,
                               np.vstack([mag_obs, color]),
                               target, p0=[0.0, 0.0])
        c_0_fit, c_1_fit = popt
        c_0_err, c_1_err = np.sqrt(np.diag(pcov))
        fitted = mag_obs + c_0_fit + c_1_fit * color
    else:
        def model_wrapper(mag, c_0):
            return mag + c_0

        popt, pcov = curve_fit(model_wrapper, mag_obs, target, p0=[0.0])
        c_0_fit  = float(popt[0])
        c_0_err  = float(np.sqrt(pcov[0, 0]))
        c_1_fit  = 0.0
        c_1_err  = 0.0
        fitted   = mag_obs + c_0_fit

    residuals = target - fitted
    rms = float(np.sqrt(np.mean(residuals**2)))

    data['mag_model'] = fitted
    data['residual']  = residuals

    return {'c_0': c_0_fit, 'c_0_err': c_0_err,
            'c_1': c_1_fit, 'c_1_err': c_1_err,
            'rms': rms, 'table': data}




def do_fig(xtab, band='R', outroot='', catalog='gaia', use_color=False):
    '''
    Plot results.  The top two panels plot the
    magnitudes as measured by aperstats, assuming
    a zeropoint of 28

    The bottom two panels plot the fitted
    fluxes
    '''

    # outdir='./Figs_phot%s' %  XDIR

    ref_label = 'SMASH' if catalog.lower() == 'smash' else 'Gaia'

    if band=='G':
        color_label='B-R'
        color_label='G-R'   # using this fits produced a better fit regardless
    else:
        color_label='G-R'

    # os.makedirs(outdir,exist_ok=True)
    plt.figure(1,(9,8))
    plt.clf()
    plt.subplot(2,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['G'],xtab['phot_mag'],marker='.',alpha=.05,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['G'],-xtab['phot_mag'],marker='.',alpha=.05,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label(color_label)
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0)
        plt.xlabel(f'{ref_label} G mag')
    else:
        plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05)
        plt.xlabel(f'{ref_label} R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')


    plt.ylim(14,22)
    plt.xlim(14,22) 



    plt.tight_layout()
    plt.subplot(2,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['R'],xtab['mag_model'],marker='.',alpha=.05,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['R'],-xtab['mag_model'],marker='.',alpha=.05,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label(color_label)
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0) 
    else:
        plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05)
    plt.xlabel(f'{ref_label} R mag')
    plt.ylabel('Corrected DECam mag')
    plt.plot([11,24],[11,24],'k-')
    plt.ylim(14,22)
    plt.xlim(14,22)



    plt.subplot(2,2,3)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']-xtab['R'],marker='.',alpha=.01,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']+xtab['R'],marker='.',alpha=.01,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label(color_label)
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)
    plt.xlabel(f'{ref_label} G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')

    plt.ylim(-2,2)
    plt.xlim(14,22)

    plt.subplot(2,2,4)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    under=xtab[xtab['phot_mag']>0]
    # plt.text(16,1.5,'Under %d Over %d' % (len(under),len(xtab)-len(under)))
    sc=plt.scatter(xtab['R'],xtab['mag_model']-xtab['R'],marker='.',alpha=.01,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['R'],xtab['mag_model']+xtab['R'],marker='.',alpha=.01,c=xtab['target_color'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label(color_label)
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)
    plt.xlabel(f'{ref_label} R mag')
    plt.ylabel('Corrected DECam mag')
    plt.plot([11,24],[0,0],'k-')
    plt.ylim(-2,2)
    plt.xlim(14,22)

    if outroot!='':
        plt.suptitle('Band %s photometry of %s (%s)' % (band, outroot, ref_label))

    os.makedirs('FigZero',exist_ok=True)

    plt.tight_layout()

    if outroot!='':
        cat_suffix   = '.smash' if catalog.lower() == 'smash' else '.gaia'
        color_suffix = '.color' if use_color else ''
        plt.savefig('FigZero/%s_%s%s%s.png' % (band, outroot, cat_suffix, color_suffix))
    plt.close()


def get_filter_from_filename(filename):
    words=filename.split('.')
    name=words[-2]
    print('The file type is %s',name)
    return name



def do_one(filename='TabPhot/c4d_241122_023910_ooi_N673_v1.fits', option='R',
           catalog='gaia', use_color=False):
    '''
    Process a single photometry table and return fit results.

    Parameters
    ----------
    filename : str
        Path to photometry FITS table (output of MefPhot).
    option : str
        Reference band to fit against: 'R' (default) or 'G'.
    catalog : str
        Reference catalog: 'gaia' (default) or 'smash'.
    use_color : bool
        If True, include an independent color term in the fit. Default False.
        Color predictor is chosen independent of the target band:
        Gaia R: G-R, Gaia G: B-R, SMASH R: G-R, SMASH G: U-R.
    '''
    try:
        xtab=Table.read(filename)
    except:
        try:
            xtab=Table.read(filename,format='ascii.fixed_width_two_line')
        except:
            print('Could not read %s' % filename)
            return

    if option == 'R':
        xtab['target_mag']   = xtab['R']
        xtab['target_color'] = xtab['G'] - xtab['R']   # G-R independent of R
    elif option == 'G':
        xtab['target_mag'] = xtab['G']
        if catalog.lower() == 'smash':
            xtab['target_color'] = xtab['U'] - xtab['R']   # U-R independent of G
        else:
            xtab['target_color'] = xtab['B'] - xtab['R']   # B-R independent of G
    else:
        print('Unknown band: %s' % option)
        return

    mask = (~xtab['target_color'].mask) & np.isfinite(xtab['target_color'])
    xtab = xtab[mask]
    xtab = xtab[xtab['Max'] < 45000]

    try:
        phot_zero=np.median(xtab['MAGZERO'])
    except:
        phot_zero=-99.

    try:
        xfilt=xtab['Filter'][0]
        word=xfilt.split()
        xfilt=word[0]
    except:
        xfilt=get_filter_from_filename(filename)

    try:
        xtime=xtab['Exptime'][0]
    except:
        xtime=-99.

    print('Filename :', filename)
    results = fit_magnitude_model(xtab[:30000], use_color=use_color)
    print(f"c_0 = {results['c_0']:.4f} +/- {results['c_0_err']:.4f}")
    if use_color:
        print(f"c_1 = {results['c_1']:.4f} +/- {results['c_1_err']:.4f}")
    print(f"RMS = {results['rms']:.4f}")
    fitted_table = results['table']

    outroot = filename.split('/')[-1].replace('.fits', '')
    for _sfx in ('.gaia', '.smash'):
        if outroot.endswith(_sfx):
            outroot = outroot[:-len(_sfx)]
            break

    do_fig(fitted_table, option, outroot, catalog=catalog, use_color=use_color)

    return 28.+results['c_0'], results['c_0'], results['c_1'], results['rms'], phot_zero, xfilt, xtime


def do_many(filenames, band='G', outroot='MagZero', catalog='gaia', use_color=False):
    zz=[]
    cc0=[]
    cc1=[]
    rrms=[]
    header_zero=[]
    xfilter=[]
    xtime=[]
    root=[]
    for one in filenames:
        try:
            zero,c_0,c_1,rms,hzero,xfilt,xt=do_one(one, band, catalog=catalog, use_color=use_color)
            zz.append(zero)
            cc0.append(c_0)
            cc1.append(c_1)
            rrms.append(rms)
            header_zero.append(hzero)
            xfilter.append(xfilt)
            xtime.append(xt)
            one_root=one.split('/')[-1].replace('.fits','')
            root.append(one_root)
        except Exception as e:
            print('Failed on %s' % (one))
            print(f'Exception: {e}')

    cat_label    = 'smash' if catalog.lower() == 'smash' else 'gaia'
    color_suffix = '.color' if use_color else ''
    xtab=Table([xfilter,xtime,root,zz,cc0,cc1,rrms,header_zero,
                [cat_label]*len(root),filenames],
               names=['Filter','Exptime','Root','MagZero','c_0','c_1','rms',
                      'HdrZero','Catalog','Filename'])
    xtab['MagZero'].format='.3f'
    xtab['c_0'].format='.3f'
    xtab['c_1'].format='.3f'
    xtab['rms'].format='.3f'
    xtab['HdrZero'].format='.3f'
    outfile = '%s.%s.%s%s.txt' % (outroot, band, cat_label, color_suffix)
    if os.path.isfile(outfile):
        qtab=Table.read(outfile,format='ascii.fixed_width_two_line')
        i=0
        select=[]
        while i<len(qtab):
            j=0
            name=qtab['Root'][i]
            while j<len(xtab):
                if name==xtab['Root'][j]:
                    break
                j+=1
            if j==len(xtab):
                select.append(i)
            i+=1
        if len(select)>0:
            xtab=vstack([qtab[select],xtab])

    xtab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return


def steer(argv):
    '''
    usage: ZeroCalc.py [-h] [-R] [-G] [-smash] file1.fits file2.fits ...
    '''

    filenames=[]
    band='R'
    outroot='MagZero'
    catalog='gaia'
    use_color=False

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:2]=='-G':
            band='G'
        elif argv[i][:2]=='-R':
            band='R'
        elif argv[i][:6]=='-smash':
            catalog='smash'
        elif argv[i][:6]=='-color':
            use_color=True
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        elif argv[i][0]=='-':
            print('Error: unknown switch: ',argv)
            return
        elif argv[i].count('fits'):
            filenames.append(argv[i])
        else:
            print('Error: Cannot parse command line : ',argv)
            return
        i+=1

    do_many(filenames, band, outroot, catalog=catalog, use_color=use_color)






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)

