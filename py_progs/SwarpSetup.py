#!/usr/bin/env python
# coding: utf-8

'''

Create Combined images of the different filters in the field using swarp

This can only be run after PrepFiles and SumFiles have benn run.    The routines here generate 
inputs to run swarp which combines the individual CCD images into tile images.   This is normally
carried out inside a Jupyter scripts.

It is also possible to run swarp once the inputs have been generated.   This can be done
either in a Jupyter script or from the command line.

To run swarp from the command line (one must be in the normal run directory) since
a particular directory structure is assumed here)

Usage:   SwarpSetup.py [-all] [-bsub] field [tiles]

where -all will cause swarp to be run on all 16 tiles.  With these inputs, the routine will
use the files ending in _sw.tab to set up run files

and -bsub directs the routine to use data for which an addtioal backgound subtraction
algorithm has been used.  In this case the Swarp commmands are written and to
the DECam_SWARP/field/tile_b directory, so that data results which are
background sutbracted and those that are not can be compared.

If one wants to run only 1 or a few tiles then the command will be something like

Swarp.py LMC_c42  T01 T02 T03




'''


import os, stat
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from glob import glob
import timeit
import subprocess
from log import *

exec_dir=os.getcwd()


SWARPDIR='DECam_SWARP'
PREPDIR=os.path.abspath('DECam_PREP')


def get_sum_tab(field='LMC_c42',tile='T07'):
    '''
    Read a table that summarizes information about
    the various files athat are used for a tile
    
    The location of the table files is currently hardocaded
    '''
    
    xname='Summary/%s_%s.txt' % (field,tile)
    try:
        xtab=ascii.read(xname)
    except:
        print('Error: Could not find %s' % xname)
        raise IOError
    return xtab



def summarize(field='LMC_c42',tile='T07'):
    '''
    Summarize the various fits files that are relevant to
    a specific tile
    '''


    try:
        x=get_sum_tab(field,tile)
    except:
        return -1, -1

    ra=np.average(x['CENRA1'])
    dec=np.average(x['CENDEC1'])

    print('The center of this tile is %.5f  %.5f' % (ra,dec))


    ha=x[x['FILTER']=='N662']
    s2=x[x['FILTER']=='N673']
    r=x[x['FILTER']=='r']
    n708=x[x['FILTER']=='N708']    
    print('Ha  images  : %3d' % len(ha))
    print('SII images  : %3d' % len(s2))
    print('R   images  : %3d' % len(r))
    print('N708 images : %3d' % len(n708))
    
    print('\nHa')
    if len(ha) >0:
        times,counts=np.unique(ha['EXPTIME'],return_counts=True)

        records=[times,counts]
        xtab=Table(records,names=['EXPTIME','Number'])
        xtab['FILTER']='N662'
        ztab=xtab.copy()

        i=0
        while i<len(times):
            print('   exp %8s  no %4d' % (times[i],counts[i]))
            i+=1
    else:
        print('This tile contains no Ha images')

    print('\nSII')

    if len(s2)>0:
        times,counts=np.unique(s2['EXPTIME'],return_counts=True)

        records=[times,counts]
        xtab=Table(records,names=['EXPTIME','Number'])
        xtab['FILTER']='N673'
        try:
            ztab=vstack([ztab,xtab])
        except:
            ztab=xtab.copy()

        i=0
        while i<len(times):
            print('   exp %8s  no %4d' % (times[i],counts[i]))
            i+=1
    else:
        print('This tile contains no SII images')
        
    print('\nR')
    if len(r)>0:
        times,counts=np.unique(r['EXPTIME'],return_counts=True)

        records=[times,counts]
        xtab=Table(records,names=['EXPTIME','Number'])
        xtab['FILTER']='r'
        try:
            ztab=vstack([ztab,xtab])
        except:
            ztab=xtab.copy()

        i=0
        while i<len(times):
            print('   exp %8s  no %4d' % (times[i],counts[i]))
            i+=1    
    else:
        print('This tile contains no R band images')
        
    print('\nN708')

    if len(n708)>0:
        times,counts=np.unique(n708['EXPTIME'],return_counts=True)

        records=[times,counts]
        xtab=Table(records,names=['EXPTIME','Number'])
        if len(xtab)>0:
            xtab['FILTER']='N708'
        try:
            ztab=vstack([ztab,xtab])
        except:
            ztab=xtab.copy()

        i=0
        while i<len(times):
            print('   exp %8s  no %4d' % (times[i],counts[i]))
            i+=1
    else:
        print('This tile contains no N708 images')

    return ztab  



def create_swarp_dir(field='LMC_c42',tile='T07',bsub=False):
    '''
    Create a diretory for the swarp outputs if it does not exist

    Note that this routine does not require the field name or
    tile to be anything specific.  It just creates a subdirecroy
    of SWARPDIR
    '''
    if bsub==False:
        outdir='%s1/%s/%s/' % (SWARPDIR,field,tile)
    else:
        outdir='%s2/%s/%s/' % (SWARPDIR,field,tile)

    if os.path.isdir(outdir)==False:
        os.makedirs(outdir)
    return outdir



xdefault='''
# Default configuration file for SWarp 2.38.0
# EB 2019-07-10
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME          %s.fits      # Output filename
WEIGHTOUT_NAME       %s.weight.fits # Output weight-map filename

HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
CENTER                %f, %f  # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.27            # Pixel scale
IMAGE_SIZE          8700,8700          # Image size (0 = AUTOMATIC)

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          N               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE             2048             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)

#------------------------------ Miscellaneous ---------------------------------
DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
'''





def create_swarp_command(field='LMC_c42',tile='T07',filt='Ha',exp=[800],defaults=xdefault,bsub=False):
    '''
    Generate the inputs necessary to run swarp where exp corresponds to one
    or more exposure times for a particular filter and tile.  

    230619 - Added code to allow commands to be created in the _b directory if bsub=True
    '''
    print('\n### Creating swarp inputs for %s tile %s and filter %s' % (field,tile,filt))
    x=get_sum_tab(field,tile)
    xx=x[x['FILTER']==filt]
    if len(xx)==0:
        print('There are no observations with filter %s')
        return
    
    if type (exp)== int:
        exp=[exp]
    
    ra=np.average(x['CENRA1'])
    dec=np.average(x['CENDEC1'])

    i=0
    while i<len(exp):
        xxx=xx[xx['EXPTIME']==exp[i]]
        if len(xxx)==0:
            print('There are no observations with exposure %s' % exp[i])
        else:
            print('There are %d observations with exposure %d' % (len(xxx),exp[i]))
        if i==0:
            xxxx=xxx.copy()

        else:
            xxxx=vstack([xxxx,xxx])
            
        i+=1
    
    if len(xxxx)==0:
        print('There were no files with exposure times given by ',exp)
        return
    else:
        print('There will be %d files summed' % len(xxxx))
        

    xtile=tile

    # if bsub==True:
    #     xtile='%s_b' % tile


    xdir=create_swarp_dir(field,xtile,bsub)
    
    

    i=0
    exp_string=''
    while i<len(exp):
        xtime='.t%03d' %exp[i]   
        exp_string+=xtime
        i+=1
    
    root='%s_%s.%s' % (field,tile,filt)  
    root+=exp_string
    # print('root ',root)    
    
    filelist='%s.in' % root
    name='%s/%s' % (xdir,filelist)
    
    # print('Filelist ',name)
    
    f=open(name,'w')
    for one in xxxx:
        xname='%s/%s/%s/%s'% (PREPDIR,field,xtile,one['Filename'])
        f.write('%s\n' % xname)
    f.close()
    
    default_name='%s.default' % (root)
    # print('swarp ',default_name)
    f=open('%s/%s' % (xdir,default_name),'w')
    
    # print('Centering on ',ra,dec)
    zdefault=defaults % (root,root,ra,dec)
    f.write(zdefault)
    f.close()
    
    f=open('%s/%s.run' % (xdir,root),'w')
    f.write('swarp @%s -WEIGHTOUT_NAME  /dev/null -c %s \n' % (filelist,default_name))
    f.close()
    
    os.chmod('%s/%s.run' % (xdir,root),stat.S_IRWXU)
    print('### Finished swarp inputs for %s tile %s and filter %s\n' % (field,tile,filt))
    return   
    
    
def create_commmands_for_one_tile(field='LMC_c42',tile='T07',defaults=xdefault,bsub=False):
    '''
    Generate the standard set of commands for on tile, based on all of the separate
    exposure times that exist for the cell

    230619 - Added code to try to handle the situation where inputs are 
    somewhat inconsistent, assuming that if the tile is given as _b, one 
    actually wants to use the background subtracted data

    '''

    xtab=summarize(field,tile)

    # Homoogenize the inputs in a situation where _b has been added to the file name
    if tile.count('_b'):
        tile=tile.replace('_b','')
        bsub=True

    # tabfile='Summary/%s_%s.txt' % (field,tile)
    # try:
    #     xtab=ascii.read(tabfile)
    # except:
    #     print('Error: Could not read %s' % tabfile)
    #     raise IOError

    for one in xtab:
        create_swarp_command(field=field,tile=tile,filt=one['FILTER'] ,exp=[one['EXPTIME']],defaults=defaults,bsub=bsub)
    return




def steer(argv):
    '''
    This is just a steering routine for running swarp on one or more
    tiles from the command line

    '''
    field=''
    tiles=[]
    xall=False
    redo=True
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
        elif argv[i][0]=='-':
            print('Error: Unknwon switch  %s' % argv[i])
            return
        elif field=='':
            field=argv[i]
        else:
            tiles.append(argv[i])
        i+=1

    if xall:
        # Assumme all directories with fits files should be searched
        # Assume we could have both background subtracted and non
        # background subtracted data to deal with
        xfiles=glob('%s/%s/*/*.fits' % (PREPDIR,field))
        xdirs=[]
        for one in xfiles:
            words=one.split('/')
            xdirs.append(words[-2].replace('_b',''))

        tiles=np.unique(xdirs)


    open_log('%s.log' % field,reinitialize=False)
    for one in tiles:
        create_commmands_for_one_tile(field,one,bsub=bsub)
        if bsub:
            log_message('SwarpSetup for %s %s for modified background images' % (field,one))
        else:
            log_message('SwarpSetup for %s %s for original images' % (field,one))
    close_log()



    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


