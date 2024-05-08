#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Get flux corresponding to 1DN through various filters


Command line usage (if any):

    usage: ZeroPoint.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240505 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import random


from PhotCompare import get_gaia_spec


                                                            

def get_flux(xtab , wavelength=6563):
    '''
     Get the flux of a star at a particular wavelength
     
    '''
    
    i=0
    while xtab['WAVE'][i] < wavelength and i<len(xtab):
        i+=1
    #print(xtab['WAVE'][i])
    frac=(wavelength-xtab['WAVE'][i-1])/(xtab['WAVE'][i]-xtab['WAVE'][i-1])
    # print(frac)
    flux=(1-frac) * xtab['FLUX'][i-1]+frac*xtab['FLUX'][i]
    # print(flux)

    return flux



def tab_remove_bad(qtab,key='teff'):
    '''
    Return a table with masked values of a column removed
    '''
    try:
        filtered_indices= ~qtab[key].mask
    except:
        qtab.info()

    filtered_tab=qtab[filtered_indices]
    return filtered_tab
    
    

def gaia_choose(filename='TabPhot/xmatch_Gaia.LMC_c32_T07_LMC_c32_T07.N673.t800_phot.txt',Rmin=14,Rmax=18,tmin=8000,tmax=15000, nmax=1000):

    try:
        xtab=ascii.read(filename)
        # xtab.info()
    except:
        print('Error: Could not read %s' % filename)
      
    ftab=tab_remove_bad(xtab,'teff')
    
    ftab=ftab[ftab['R']>Rmin]
    ftab=ftab[ftab['R']<Rmax]
    
    ftab=ftab[ftab['teff']>tmin]
    ftab=ftab[ftab['teff']<tmax]
    
    
    if len(ftab)>nmax:
        
        print('Selecting %d from %d possibilities' % (nmax,len(ftab)))

        # Generate a list of indices corresponding to the list length
        indices = list(range(len(ftab)))

        # Randomly select num_elements_to_select indices without replacement
        selected_indices = random.sample(indices, nmax)
        # print(np.sort(selected_indices))
        ztab=ftab[selected_indices]
        return ztab
    else:
        print('There are only %d objects when wanted %d' % (len(ftab),nmax))
        return ftab
        
        

def get_gaia_flux(xtab,key='Ha',wave=6563):
    
    i=0
    select=[]
    ha_values=[]
    s2_values=[]
    while i<len(xtab):
        one=xtab[i]
        xtest='%s' % one['Source_name']
        spec_tab=get_gaia_spec(xtest)
        if len(spec_tab)>0:
            f_ha=get_flux(spec_tab,6563)
            ha_values.append(f_ha)
            f_s2=get_flux(spec_tab,6720)
            s2_values.append(f_s2)
            select.append(i)
            print('Succeeded for ', one['Source_name'],f_ha,f_s2)
        # else:
        #    print('Failed    for ',  one['Source_name'])
        i+=1
        
    print(len(select))
    ztab=xtab[select]
    ztab['Ha']=ha_values
    ztab['S2']=s2_values
    if len(ztab)>0:
        return ztab
    else:
        return None
        
        
def do_plot(ztab, outroot='foo'):
    plt.figure(1,(8,8))
    plt.clf()
    plt.subplot(3,2,1)
    plt.hist(ztab['Ha'])
    plt.subplot(3,2,3)
    plt.loglog(ztab['Ha'],ztab['Net'],'o')
    plt.subplot(3,2,5)
    plt.hist(ztab['S2']/ztab['Net'],range=[0,5e-20])
    plt.subplot(3,2,2)
    plt.hist(ztab['S2'])
    plt.subplot(3,2,4)
    plt.loglog(ztab['S2'],ztab['Net'],'o')
    plt.subplot(3,2,6)
    plt.hist(ztab['S2']/ztab['Net'],range=[0,5e-20])
    plt.tight_layout()
    plt.savefig('%s.png' % outroot)

def do_one(xmatch_file,Rmin=14,Rmax=18,tmin=8000,tmax=15000, nmax=1000):


    ztab=gaia_choose(xmatch_file,Rmin,Rmax,tmin,tmax, nmax)
    final=get_gaia_flux(xtab=ztab,key='Ha',wave=6563)

    words=xmatch_file.split('/')
    outroot='Figs_phot/%s' % (words[-1])
    do_plot(final,outroot)
        
        
    
def steer(argv):


    i=1
    files=[]
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: Unknown switch in command line: ',argv)
        else:
            files.append(argv[i])

        i+=1

    for one in files:
        print(one)
        do_one(one)









# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
