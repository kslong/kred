#!/usr/bin/env python
# coding: utf-8

"""Analyze the forced photometry associated with various filters

Space Telescope Science Institute

Synopsis
--------

Analyze the forced photometry associated with various filters

Command Line Usage
------------------

::

    usage: PhotAnal.py filename

Description
-----------

Primary routines:

    doit

Primary Routines
----------------

doit

Notes
-----

History:

251123 ksl Coding begun

Version History
---------------

251123 ksl
    Coding begun

"""


import sys
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.table import Table,join



def remove_masked_rows_single_column(tab, column_name):
    """
    Remove rows where the specified column has masked values.
    
    Parameters
    ----------
    tab : astropy.table.Table
        Input table with masked column(s)
    column_name : str
        Name of the column to check for masked values
    
    Returns
    -------
    filtered_tab : astropy.table.Table
        Table with masked rows removed
    """
    # Get mask for the column (True where data is valid)
    if hasattr(tab[column_name], 'mask'):
        valid_mask = ~tab[column_name].mask
    else:
        # Column is not masked, keep all rows
        valid_mask = np.ones(len(tab), dtype=bool)
    
    return tab[valid_mask]




def get_stats(table):
    qtab=remove_masked_rows_single_column(table,'R')
    qtab['G-R']=qtab['G']-qtab['R']
    qtab['delta_mag_simple']=qtab['phot_mag_simple']-qtab['R']
    qtab['delta_mag']=qtab['phot_mag']-qtab['R']
    mag_min=17
    mag_max=20
    color_cut=0.1
    mask=( qtab['G-R']>color_cut) & (qtab['R']>mag_min) & (qtab['R']<mag_max)
    ztab=qtab[mask]
    # print(len(ztab))
    med_simple=np.median(ztab['delta_mag_simple'])
    med=np.median(ztab['delta_mag'])
    ave_simple=np.average(ztab['delta_mag_simple'])
    ave=np.average(ztab['delta_mag'])
    std_simple=np.std(ztab['delta_mag_simple'])
    std=np.std(ztab['delta_mag'])
    
    return med_simple,ave_simple,std_simple,med,ave,std
    




def do_all(filenames):
    xmed_simple=[]
    xmed=[]
    xave_simple=[]
    xstd_simple=[]
    xave=[]
    xstd=[]
    for one in filenames:
        xtab=ascii.read(one)
        qtab=remove_masked_rows_single_column(xtab,'R')
        med_simple,ave_simple,std_simple,med,ave,std,=get_stats(qtab)
        xmed_simple.append(med_simple)
        xmed.append(med)
        xave_simple.append(ave_simple)
        xstd_simple.append(std_simple)
        xave.append(ave)
        xstd.append(std)
    ptab=Table([filenames,xmed_simple,xave_simple,xstd_simple,xmed,xave,xstd],names=['filename','med_simple','ave_simple','std_simple','med','ave','std'])
    return ptab
        


def get_names(files):
    names=[]
    for one_file in files:
        word=one_file.split('_x_')
        xname=word[0].replace('TabPhot_','')
        xname=xname.split('/')[-1]
        xname=xname.replace('_phot','')
        names.append('%s' % xname)
    return names



def doit(xdir='DECam_PREP2/LMC_c42/T14'):

    xxdir=xdir.replace('/','-')

    xfile='TabPhot_%s/*_x_*.txt' % xxdir
    ifile='Image_Sum_%s.txt' % xxdir

    # xfile='TabPhot_DECam_PREP2-LMC_c42-T15/*_x_*.txt'
    # ifile='Image_Sum_DECam_PREP2-LMC_c42-T15.txt'

    # xfile='TabPhot_DECam_PREP2-LMC_c42-T14/*_x_*.txt'
    # ifile='Image_Sum_DECam_PREP2-LMC_c42-T14.txt'



    files=glob(xfile)
    if len(files)==0:
        print('Error: not files found in : ',xfile)
        return 

    try:
        imsum=ascii.read(ifile)
    except:
        pritnt('Could not raad: ',ifile)
        rerurn


    print('Evaluating statistics for %d files' % len(files))
    ptab=do_all(files)
    names=get_names(ptab['filename'])
    ptab['filename']=names

    print('Finished evaluating statistics')


    xname=[]
    for one_name in imsum['filename']:
        foo=one_name.split('/')[-1]
        xname.append(foo.replace('.fits',''))

    imsum['filename']=xname

    print('imsum')
    print(imsum[0:10])

    print('ptab')
    print(ptab[:10])

    all=join(imsum,ptab,join_type='left')

    all.info()




    s2=all[all['Filter']=='N673']
    ha=all[all['Filter']=='N662']
    r=all[all['Filter']=='r']
    n708=all[all['Filter']=='N708']



    ha=ha[ha['Exptime']>100]
    s2=s2[s2['Exptime']>100]
    n708=n708[n708['Exptime']>100]
    r=r[r['Exptime']>10]


    plt.figure(1,(12,6))
    plt.subplot(1,2,1)
    foo=plt.hist([ha['med_simple'],s2['med_simple'],n708['med_simple'],r['med_simple']],100,range=(-0.5,.1),cumulative=True,histtype='step',density=True,label=['Ha','SII','N708','R'])
    plt.legend()
    plt.subplot(1,2,2)
    foo=plt.hist([ha['med'],s2['med'],n708['med'],r['med']],100,range=(-0.5,.1),cumulative=True,histtype='step',density=True,label=['Ha','SII','N708','R'])
    plt.legend()
    plt.suptitle(xxdir)
    plt.tight_layout()
    plt.savefig('%s_phot.png' % xxdir )


    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: PhotAnal.py filename')




