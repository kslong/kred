#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Produce pure emision (and continuum) images via 
simple subtraction


Command line usage (if any):

    usage: CleanStars.py [-all] [-bsub] field [tiles]

    where 
        -all will cause CleanStars to be run on all 16 tiles in a field
            and if -all is not given, then one or more tiles should be listed
        -bsub will search for inputs in the DECam_SWARP2 directories which
            have "better" background matching, while if it is absence
            one will use the data in the DECam_SWARP directories, which use
            the background form the overal fields




Description:  

    This routine produces simple images that
    have been subtracted to remove the coinuum
    (or from the r-band continuum the emissionlines)

    It must be run from the top level kred directory

    The routine is currently hardwired to use the
    longest exposures of a particular filter type.

    
    The routine is NOT sophisticated.  

    For creating the r-band results, we first subtract the two emission
    line images from the r-band images; this should produce an image
    with the emission lines removed.  We then subtract this from the emission
    line images.  There is some scaling involved.

    For creating the N708 results, we do a simple subtraction.  


    (Note that before we actually do the subtractions, we calculated 
    a background level n the continuum images.   This is intended
    to leave whatever backgrouground level in the emission line images
    unaffected.)


    The files that are produced are nominally the following (for LMC_c42_T07):
    LMC_c42_T07.ha_sub_r.fits     - pure ha image using r-band for continuum
    LMC_c42_T07.s2_sub_r.fits     - pure s2 image using r-rand for continuum
    LMC_c42_T07.r_sub.fits        - pure r-band image after emission lines are subtracted       
    LMC_c42_T07.ha_s2_sub.fits    - the emission line portion of the r-band image 

    LMC_c42_T07.ha_sub_n708.fits  - pure ha image based on subtrcting the n708 image
    LMC_c42_T07.s2_sub_n708.fits  - pure ha image based on subtrcting the n708 image
    LMC_c42_T07.n708_sub.fits     - n708 image (after subtracting a biased median) 


    The outputs are written to a subdirectory of DECam_SUB


Primary routines:

    doit
    make_rband_subtractions
    make_n708_subtractions
    make_n540_subtractions

Notes:
                                       
History:

230717 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np


from astropy.io import fits
from astropy.table import Table
from glob import glob

import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

from log import *




def fits_deep_copy(orig):
    '''
    Make a deep copy of a fits image
    
    Astropy does not have a way, apparently to make 
    a deep copy of an image; instead it makes shallow
    copies, which means that often one does not know
    what one is dealing with. 
    
    The basic problem here is that if one does not
    make deep copies, one will not necessary know
    what image you are working with, one that is
    pristine or one that has been modified
    
    This little routine
    avoids this issue.
    
    
    '''
    xhead=orig[0].header.copy()
    xdata=orig[0].data.copy()
    hdu=fits.PrimaryHDU(xdata)
    x=fits.HDUList([hdu])
    x[0].header=xhead
    return x




def make_rband_subtractions(ha='data/LMC_c42_T07.N662.t800.fits',s2='data/LMC_c42_T07.N673.t800.fits',r='data/LMC_c42_T07.r.t060.fits',outroot='test'):
    
    ha_exists=False
    s2_exists=False
    r_exists=False


    
    try:
        zha=fits.open(ha)
        ha_exists=True
    except:
        print('CleanStars: The   Ha file could not be openened: %s' % ha)
    
    try:
        zs2=fits.open(s2)
        s2_exists=True
    except:
        print('CleanStars: The  SII file could not be openened: %s' % s2)
         
    try:
        zr=fits.open(r)
        r_exists=True
    except:
        print('CleanStars: The    R file could not be openened: %s' % r)    
        
    if ha_exists==False or s2_exists==False or r_exists==False:
        return
    

    if ha_exists:
        mean,median,std=sigma_clipped_stats(zha[0].data,sigma_lower=2,sigma_upper=1,grow=3)
        print('CleanStars: ha:  ',mean,median,std)
        zha[0].data-=median

    
    if s2_exists:
        mean,median,std=sigma_clipped_stats(zs2[0].data,sigma_lower=2,sigma_upper=1,grow=3)
        print('CleanStars: s2:  ',mean,median,std)
        zs2[0].data-=median

    if outroot=='':
        word=ha.split('.')
        outroot=word[0]
        # Use the directory path/first_word as the root, which is the usual case

    narrow=fits_deep_copy(zha)
    narrow[0].data=zha[0].data+zs2[0].data-2*zr[0].data
    narrow.writeto(outroot+'.ha_s2_sub.fits',overwrite=True)

    r_pure=fits_deep_copy(zr)
    r_pure[0].data-= (130/1500)* narrow[0].data
    
    mean,median,std=sigma_clipped_stats(r_pure[0].data,sigma_lower=2,sigma_upper=1,grow=3)
    print(mean,median,std)
    r_pure[0].data-=median
    
    r_pure[0].header['PROCTYPE']='LineSubtracted'
    r_pure[0].header['SFILTER']=('Ha + SII','Filter of image used for subtracttion')
    r_pure.writeto(outroot+'.r_sub.fits',overwrite=True) 
    
    
    
    x=fits_deep_copy(zha)
    x[0].data-=r_pure[0].data
    x[0].header['PROCTYPE']='StarSubtracted'
    x[0].header['SFILTER']=(zr[0].header['FILTER'],'Filter of image used for subtracttion')
    x.writeto(outroot+'.ha_sub_r.fits',overwrite=True)        

    x=fits_deep_copy(zs2)
    x[0].data-=r_pure[0].data
    x[0].header['PROCTYPE']='StarSubtracted'
    x[0].header['SFILTER']=(zr[0].header['FILTER'],'Filter of image used for subtracttion')
    x.writeto(outroot+'.s2_sub_r.fits',overwrite=True)      
    
    return

def make_n540_subtractions(o3='data/LMC_c42_T07.N501.t800.fits', n540='data/LMC_c42_T07.N540.t300.fits',outroot='test3'):
    '''
    For subtraction of N540 from N501 and assume no emission line 
    contamination, and since everything is scaled to the same level
    we just do a straight subtraction.
    '''
    o3_exists=False
    n540_exists=False

    try:
        zo3=fits.open(o3)
        o3_exists=True
    except:
        print('CleanStars: The OIII file could not be opened: %s' % o3)

    try:
        zn540=fits.open(n540)
        n540_exists=True
    except:
        print('CleanStars: The N540 file could not be opened: %s' % n540)
        return

    if o3_exists:
        mean,median,std=sigma_clipped_stats(zo3[0].data,sigma_lower=2,sigma_upper=1,grow=3)
        print('CleanStars: oiii:  ',mean,median,std)
        zo3[0].data-=median

    mean,median,std=sigma_clipped_stats(zn540[0].data,sigma_lower=2,sigma_upper=1,grow=3)
    print('CleanStars: n540: ',mean,median,std)
    zn540[0].data-=median

    zn540[0].header['PROCTYPE']='LineSubtracted'
    zn540.writeto(outroot+'.n540_sub.fits',overwrite=True)

    if outroot=='':
        word=o3.split('.')
        outroot=word[0]
        # Use the directory path/first_word as the root, which is the usual case


    if o3_exists and n540_exists:
        # zha[0].data-=zn708[0].data
        zo3[0].header['PROCTYPE']='StarSubtracted'
        zo3[0].header['SFILTER']=(zn540[0].header['FILTER'],'Filter of image used for star subtraction')
        zo3[0].data-=zn540[0].data
        zo3.writeto(outroot+'.o3_sub_n540.fits',overwrite=True)



def make_n708_subtractions(ha='data/LMC_c42_T07.N662.t800.fits',s2='data/LMC_c42_T07.N673.t800.fits',n708='data/LMC_c42_T07.N708.t400.fits',outroot='test2'):
    '''
    For narrow band subtraction we assume there is no emission line
    contamination, and since everything is scaled to the same level
    we just do a straight subtractin
    '''
    ha_exists=False
    s2_exists=False

    n708_exists=False
    
    try:
        zha=fits.open(ha)
        ha_exists=True
    except:
        print('CleanStars: The   Ha file could not be opened: %s' % ha)
    
    try:
        zs2=fits.open(s2)
        s2_exists=True
    except:
        print('CleanStars: The  SII file could not be opened: %s' % s2)
      
      
    try:
        zn708=fits.open(n708)
        n708_exists=True
    except:
        print('CleanStars: The N708 file could not be opened: %s' % n708)
        return

    if ha_exists:
        mean,median,std=sigma_clipped_stats(zha[0].data,sigma_lower=2,sigma_upper=1,grow=3)
        print('CleanStars: ha:  ',mean,median,std)
        zha[0].data-=median

    
    if s2_exists:
        mean,median,std=sigma_clipped_stats(zs2[0].data,sigma_lower=2,sigma_upper=1,grow=3)
        print('CleanStars: s2:  ',mean,median,std)
        zs2[0].data-=median

    mean,median,std=sigma_clipped_stats(zn708[0].data,sigma_lower=2,sigma_upper=1,grow=3)
    print('CleanStars: n708: ',mean,median,std)
    zn708[0].data-=median

    zn708[0].header['PROCTYPE']='LineSubtracted'
    zn708.writeto(outroot+'.n708_sub.fits',overwrite=True) 
    
    if outroot=='':
        word=ha.split('.')
        outroot=word[0]
        # Use the directory path/first_word as the root, which is the usual case

        
    if ha_exists and n708_exists:
        # zha[0].data-=zn708[0].data
        zha[0].header['PROCTYPE']='StarSubtracted'
        zha[0].header['SFILTER']=(zn708[0].header['FILTER'],'Filter of image used for star subtraction')
        zha[0].data-=zn708[0].data
        zha.writeto(outroot+'.ha_sub_n708.fits',overwrite=True)

    if s2_exists and n708_exists:
        # zs2[0].data-=zn708[0].data
        zs2[0].header['PROCTYPE']='StarSubtracted'
        zs2[0].header['SFILTER']=(zn708[0].header['FILTER'],'Filter of image used for star subtraction')
        zs2[0].data-=zn708[0].data
        zs2.writeto(outroot+'.s2_sub_n708.fits',overwrite=True)
        

def doit(xdir='data',outdir='data'):
    '''
    Routine to produce star subtracted images
    '''

    files=glob('%s/*.fits' % xdir)
    if len(files)==0:
        print("CleanStars: No files found for %s" % xdir)
        return

    # Now the files are of the following form: LMC_c30_T01.N673.fits
    print(files)
    records=[]
    qfile=[]
    qroot=[]
    qfilt=[]
    for one in files:
        word=one.split('/')
        filename=word[-1]
        word=filename.split('.')
        root=word[0]
        xfilt=word[1]
        record=[one,root,xfilt]
        records.append(record)
        qfile.append(one)
        qroot.append(root)
        qfilt.append(xfilt)



    r=ha=s2=n708=n540=o3='none'
    xtab=Table([qfile,qroot,qfilt],names=['Filename','Root','Image'])
    xtab.sort('Image')
    xtab.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)
    print(xtab)
    zroot='test'
    for one in xtab:
        if one['Image']=='N662':
            ha=one['Filename']
            zroot=one['Root']
        elif one['Image']=='N673':
            s2=one['Filename']
            zroot=one['Root']
        elif one['Image']=='N501':
            o3=one['Filename']
            zroot=one['Root']
        elif one['Image']=='r':
            r=one['Filename']
        elif one['Image']=='N708':
            n708=one['Filename']
        elif one['Image']=='N540':
            n540=one['Filename']
        else:
            print('CleanStars: Unknown fits file: %s '% one['Filename'])

    print(ha,s2,r,n708,o3,n540)


    if r!='none':
        make_rband_subtractions(ha,s2,r,outroot='%s/%s' % (outdir,zroot))
    else:
        print('CleanStars: No r band image found')

    if n708!='none':
        make_n708_subtractions(ha,s2,r,outroot='%s/%s' % (outdir,zroot))
    else:
        print('CleanStars: No N708 image found')

    if n540!='none':
        make_n540_subtractions(o3,n540,outroot='%s/%s' % (outdir,zroot))
    else:
        print('CleanStars: No N540 image found')



def steer(argv):
    '''
    This is just a steering routine for running swarp on one or more
    tiles from the command line

    '''
    field=''
    tiles=[]
    xall=False
    bsub=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-bsub':
            bsub=True
        elif field=='':
            field=argv[i]
        else:
            tiles.append(argv[i])
        i+=1

    if xall:
        tiles=[]
        i=1
        while i<17:
            tiles.append('T%02d' % i)
            i+=1

    open_log('%s.log' % field,reinitialize=False)

    if bsub:
        xindir='DECam_SWARP2/'
        xoutdir='DECam_SUB2/'
    else:
        xindir='DECam_SWARP1/'
        xoutdir='DECam_SUB1/'

    for tile in tiles:
        indir='%s/%s/%s/' % (xindir,field,tile)
        outdir='%s/%s/%s' % (xoutdir,field,tile)

        if os.path.isdir(indir)==False:
            log_message('CleanStars: Cannot find input directory: %s ' % indir)
        else:
            if os.path.isdir(outdir)==False:
                os.makedirs(outdir)

            doit(indir,outdir)

            log_message('CleanStars: Finished  %s %s ' % (field,tile))

    close_log()


    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
