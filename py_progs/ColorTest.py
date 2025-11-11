#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

See how well images have been subtracted


Command line usage (if any):

    Usage ColorTest.py -h -dir dirname file1..

    where:
        -h prints out this help and exits
        -dir searches the directory dirname and
            all subdirectories for continuum subtracted
            emission line files, processes these images
        file1 ... if -xdir is not specified carries out
            the analysis on specific files

Description:  

    The routine access the GAIA database to find stars in the
    field.  It then carries out forced photometry on emission
    line image from which the stars were subtracted and on
    the continuum image that was used for subtraction.  It
    produces a plot that shows the amount of over-subtaction
    or under-subtraction

Primary routines:

    doit

Notes:

    The routine looks for files with specific names to decide
    what images are are those that can be analyzed.  If new
    filters are involved some changes are requried.
                                       
History:

251108 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt



import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import join
import matplotlib.cm as cm
import PhotCompare
from kred import ImageSum



def color_compare(cont_image,subtracted_image,forced=False):

    print('Making new GaiCat file')
    ra,dec,size_deg=PhotCompare.get_size(cont_image)
    gaia_file=PhotCompare.get_gaia(ra, dec, size_deg,outroot='',nmax=-1)
    # So at this point I have the Gaifile

    cont_phot=PhotCompare.do_forced_photometry(cont_image,gaia_file,'')
    print(cont_phot)
    sub_phot=PhotCompare.do_forced_photometry(subtracted_image,gaia_file,'')
    print(sub_phot)
    return cont_phot,sub_phot,gaia_file




def plot_both(gaia,final,title=''):
    plt.figure(3,(12,6))
    plt.clf()
    plt.subplot(1,2,1)
    # plt.plot(gaia['G']-gaia['R'],gaia['R'],'.',alpha=0.01)
    sc=plt.scatter(gaia['G']-gaia['R'],gaia['R'],c=gaia['G']-gaia['R'],marker='.',cmap='plasma',vmin=-1,vmax=1,alpha=0.01)
    plt.xlabel('G-R')
    plt.ylabel('R')
    plt.xlim(-2,3)
    plt.ylim(12,23)

    cbar = plt.colorbar(sc)         # use the scatter as the mappable
    cbar.set_label('G - R')

    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)


    plt.subplot(1,2,2)

    sc = plt.scatter(
        np.array(final['Net_1']),
        np.array(final['Net_2']),
        marker='.',
        c=final['G-R'],
        cmap='plasma',
        vmin=-1,
        vmax=1,
        alpha=0.01
    )

    plt.ylim(-1000, 1000)
    plt.xlim(10, 1e4)
    plt.xlabel('Continuum Flux')
    plt.ylabel('Measured Flux in Subtracted Image')

    cbar = plt.colorbar(sc)         # use the scatter as the mappable
    cbar.set_label('G - R')

    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0)

    if title=='':
        plt.suptitle("Comparison of Gaia Colors and Image Fluxes", fontsize=16)
    else:
        plt.suptitle(title, fontsize=16)

    plt.tight_layout()

def get_base_name(filename='DECam_SWARP2/LMC_c35/T06/LMC_c35_T06.r.fits'):
    word=filename.split('/')
    base=word[-1].replace('.fits','')
    return base


def find_cont(filenames):
    '''
    Locate the companion files for subtracted files
    '''
    names=[]
    for name in filenames:
        name=name.replace('SUB','SWARP')
        name=name.replace('ha_sub_r','r')
        name=name.replace('s2_sub_r','r')
        names.append(name)
    return names

def doit(continuum_file='DECam_SWARP2/LMC_c35/T06/LMC_c35_T06.r.fits',subtracted_file='DECam_SUB2/LMC_c35/T06/LMC_c35_T06.ha_sub_r.fits'):

    try:
        cont_phot,sub_phot,gaia_file=color_compare(continuum_file,subtracted_file)
    except ValueError:
        print('Unsucessful for theis combination: %s %s' % (continuum_file,subtracted_file))
        return
    cont=ascii.read(cont_phot)
    sub=ascii.read(sub_phot)
    gaia=ascii.read(gaia_file)
    gaia['id']=np.arange(len(gaia))+1
    gaia['G-R']=gaia['G']-gaia['R']

    cont=join(cont,gaia,keys=['id'])
    sub=join(sub,gaia,keys=['id'])
    final=join(cont['id','Net','G-R'],sub['id','Net'],keys=['id'])

    base_cont=get_base_name(continuum_file)
    base_sub=get_base_name(subtracted_file)
    title='GAIA Color/Flux Comparison: %s and %s' % (base_sub,base_cont)
    plot_both(gaia,final,title)
    os.makedirs('Figs_color',exist_ok=True)
    plt.savefig('Figs_color/%s.png' % base_sub)
    return


def do_dir(xdir='DECam_SUB2/LMC_c37/T07',nrow_max=-1):
    '''
    Sort out how to call doit based on the directories
    '''

    xtab=ImageSum.table_create(xdir,outname=None)
    # The two or ultimated 3 types of images we will have will all be indentied becasue the image bytpe will have a name including sub_
    select=[]
    for i in range(len(xtab)):
        if xtab['Image_type'][i].count('sub_'):
            select.append(i)

    xsub=xtab[select]
    if len(xsub)==0:
        print('Error: There were no recognized continuum-subtracted images in %s' % xdir)
        return
    names=[]
    for name in xsub['filename']:
        name=name.replace('SUB','SWARP')
        name=name.replace('ha_sub_r','r')
        name=name.replace('s2_sub_r','r')
        name=name.replace('ha_sub_n708','n708')
        name=name.replace('s2_sub_n708','n708')
        names.append(name)
    xsub['cont_file']=names
    print(xsub)
    xsub.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)
    for one_row in xsub:
        doit(one_row['cont_file'],one_row['filename'])
    return




def steer(argv):
    '''
    Run the script given choices from the command line

    Usage ColorTest.py -h -xdir whatever file1..

    '''

    nrow_max=30000
    files=[]
    xdir=''

    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(__doc__)
            return
        elif argv[i]=='-dir':
            i+=1
            xdir=argv[i]
        elif argv[i]=='-nmax':
            i+=1
            nrow_max=int(argv[i])
        elif argv[i][0]=='-':
            print('Unknown switch ',argv)
            return
        else:
            files.append(argv[i])
        i+=1


    if xdir!='':
        do_dir(xdir=xdir,nrow_max=nrow_max)
        return

    cont_files=find_cont(filenames)
    xtab=Table([files,cont_files],names=['filename','cont_file'])
    for one_row in xtab:
        doit(continuum_file=one_row['cont_file'],subtracted_file=one_row['filename'])


    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)  
    else:
        print (__doc__)
