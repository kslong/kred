#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Calculate the ZeroPoint for an image 
given one or more tables containing 
forced photometry based on the Gaia
catalog


Command line usage (if any):

    usage: ZeroCalc.py filename

Description:  

Primary routines:

    doit

Notes:

    This is one of the routines developed to see how
    consistent MAGZERO is as delivered by the community
    pipeline.
                                       
History:

251128 ksl Coding begun

'''

import sys
import matplotlib.pyplot as plt
import os
from glob import glob
from astropy.table import Table, join
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import Table

def fit_magnitude_model(data):
    """
    Fit catalogued magnitudes to the model: R = mag + c_0 + c_1 * (G-R)
    
    Parameters:
    -----------
    data : astropy.table.Table
        Table containing columns 'mag', 'R', and 'G-R'
    
    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'c_0': fitted zero-point offset
        - 'c_0_err': uncertainty in c_0
        - 'c_1': fitted color term coefficient
        - 'c_1_err': uncertainty in c_1
        - 'rms': RMS of residuals
        - 'table': input table with added 'R_model' and 'residual' columns
    """
    # Extract columns from the table
    mag_obs = data['phot_mag']
    R_cat = data['R']
    GR_color = data['G-R']
    
    # Define the model function
    def magnitude_model(mag, GR, c_0, c_1):
        return mag + c_0 + c_1 * GR
    
    # Wrapper for curve_fit
    def model_wrapper(x_data, c_0, c_1):
        mag, GR = x_data
        return magnitude_model(mag, GR, c_0, c_1)
    
    # Prepare data for fitting
    x_data = np.vstack([mag_obs, GR_color])
    
    # Perform the fit with initial guesses
    p0 = [0.0, 0.0]
    popt, pcov = curve_fit(model_wrapper, x_data, R_cat, p0=p0)
    
    # Extract fitted parameters and uncertainties
    c_0_fit, c_1_fit = popt
    c_0_err, c_1_err = np.sqrt(np.diag(pcov))
    
    # Calculate model values and residuals
    R_model = magnitude_model(mag_obs, GR_color, c_0_fit, c_1_fit)
    residuals = R_cat - R_model
    rms = np.sqrt(np.mean(residuals**2))
    
    # Add fitted values to the table
    data['R_model'] = R_model
    data['residual'] = residuals
    
    # Return results dictionary
    results = {
        'c_0': c_0_fit,
        'c_0_err': c_0_err,
        'c_1': c_1_fit,
        'c_1_err': c_1_err,
        'rms': rms,
        'table': data
    }
    
    return results




def do_fig(xtab,outroot=''):

    # outdir='./Figs_phot%s' %  XDIR

    # os.makedirs(outdir,exist_ok=True)
    plt.figure(1,(9,8))
    plt.clf()
    plt.subplot(2,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['G'],xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['G'],-xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label('G-R')
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0) 
    else:
        plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')


    plt.ylim(14,22)
    plt.xlim(14,22) 



    plt.tight_layout()
    plt.subplot(2,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['R'],xtab['R_model'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['R'],-xtab['R_model'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label('G-R')
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0) 
    else:
        plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('Corrected DECam mag')
    plt.plot([11,24],[11,24],'k-')
    plt.ylim(14,22)
    plt.xlim(14,22)  



    plt.subplot(2,2,3)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']-xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']+xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label('G-R')
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0) 
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')

    plt.ylim(-2,2) 
    plt.xlim(14,22) 

    plt.subplot(2,2,4)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    under=xtab[xtab['phot_mag']>0]
    # plt.text(16,1.5,'Under %d Over %d' % (len(under),len(xtab)-len(under)))
    sc=plt.scatter(xtab['R'],xtab['R_model']-xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['R'],xtab['R_model']+xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label('G-R')
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0) 
    plt.xlabel('Gaia R mag')
    plt.ylabel('Corrected DECam mag')
    plt.plot([11,24],[0,0],'k-')
    plt.ylim(-2,2) 
    plt.xlim(14,22)  

    if outroot!='':
        plt.suptitle(outroot)
        os.makedirs('FigZero',exist_ok=True)
        plt.suptitle(outroot)

    plt.tight_layout()

    if outroot!='':
        plt.savefig('FigZero/%s.png' % (outroot))




def do_one(filename='TabPhot/c4d_241122_023910_ooi_N673_v1.fits'):
    xtab=Table.read(filename)
    xtab['G-R']=xtab['G']-xtab['R']
    mask = (~xtab['G-R'].mask) & np.isfinite(xtab['G-R'])
    xtab = xtab[mask]

    xtab=xtab[xtab['Max']<45000]




    results=fit_magnitude_model(xtab[:30000])
    print(f"c_0 = {results['c_0']:.4f} ± {results['c_0_err']:.4f}")
    print(f"c_1 = {results['c_1']:.4f} ± {results['c_1_err']:.4f}")
    print(f"RMS = {results['rms']:.4f}")
    fitted_table = results['table']

    outroot=filename.split('/')[-1].replace('.fits','')

    do_fig(fitted_table,outroot)

    return 28.+results['c_0'],results['c_0'],results['c_1'],results['rms']


def do_many(filenames,outroot='MagZero'):
    zz=[]
    cc0=[]
    cc1=[]
    rrms=[]
    for one in filenames:
        zero,c_0,c_1,rms=do_one(one)
        zz.append(zero)
        cc0.append(c_0)
        cc1.append(c_1)
        rrms.append(rms)

    root=[]
    for one in filenames:
        one_root=one.split('/')[-1].replace('.fits','')
        root.append(one_root)
    xtab=Table([root,zz,cc0,cc1,rrms,filenames],names=['Root','MagZero','c_0','c_1','rms','Filename'])
    xtab['MagZero'].format='.3f'
    xtab['c_0'].format='.3f'
    xtab['c_1'].format='.3f'
    xtab['rms'].format='.3f'
    xtab.write('%s.txt' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    return


def steer(argv):

    filenames=[]
    
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: unknow switch:',argv)
            return
        else:
            filenames.append(argv[i])
        i+=1

    do_many(filenames)






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)

