#!/usr/bin/env python
# coding: utf-8

'''
Create inputs for combining images of the different filters in the field with Swarp

This rotutined can only be run after PrepFiles and SumFiles have benn run.    The routines here generate 
inputs to run swarp which combines the individual CCD images into tile images.   
To run swarp from the command line (one must be in the normal run directory) since
a particular directory structure is assumed here)

Usage:   SwarpSetup.py [-all] [-bsub] [-ave_pos field [tiles]

where 

    -all will cause swarp to be run on all 16 tiles.  With these inputs, the routine will
        use the files ending in _sw.tab to set up run files
    -bsub directs the routine to use data for which an addtioal backgound subtraction
        algorithm has been used.  In this case the Swarp commmands are written and to
        the DECam_SWARP2/field/tile directory and the data are taken from the DECAM_PREP2/field/tile
        directory
    -ave_pos causes the output fits file to be centered on the average position of all of the
        files which will be swarped


If one wants to run only 1 or a few tiles then the command will be something like

SwarpSetup.py LMC_c42  T01 T02 T03

Note: 

    The current default is to create input files for Swarp that are centered on
    the position designated in the configuration file, unlike what was done
    initially which was to center on the average position of all of the files
    to be swarped. To restore the old behavior use the switch -ave_pos

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
    the various files that are used for a tile
    
    The location of the table files is currently hardocaded
    '''
    
    xname='Summary/%s_%s.txt' % (field,tile)
    try:
        xtab=ascii.read(xname)
        print('SwarpSetup: Read %s' % xname)
    except:
        print('SwarpSetup: Error: Could not find %s' % xname)
        raise IOError
    return xtab


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





def create_swarp_command(field='LMC_c42',tile='T07',image='N673',defaults=xdefault,bsub=False,use_config=True):
    '''
    Generate the inputs necessary to run swarp where exp corresponds to one
    or more exposure times for a particular filter and tile.  

    230619 - Added code to allow commands to be created in the _b directory if bsub=True
    241005 - Added code to center on the config position as the default
    '''
    print('\n###SwarpSetup:  Creating swarp inputs for %s tile %s and image %s' % (field,tile,image))
    x=get_sum_tab(field,tile)
    print('create swarp comman.info)d')
    xtab.info()

    xx=x[x['Image']==image]
    if len(xx)==0:
        print('SwarpSetup: There are no observations with filter %s')
        return
    
    
    header=x.meta['comments']

    ra=-99.
    dec=-99.
    if use_config==True:
        try:
            for one in header:
                word=one.split()
                if word[0]=='RA':
                    ra=eval(word[1])
                elif word[0]=='DEC':
                    dec=eval(word[1])
                else:
                    print('SwarpSetup: Unknown header line ',one)
        except:
            print('SwarpSetup: No valid header to summary file')

    if ra==-99. or dec==-99.:
        print('SwarpSetup: Using average position of images to set up SWARP ouput image')
        ra=np.average(x['CENRA1'])
        dec=np.average(x['CENDEC1'])
    else:
        print('SwarpSetup: Using configuration file RA and DEC to set up SWARP ouput image')


    print('SwarpSetup: The center of this tile is %.5f  %.5f' % (ra,dec))




    xdir=create_swarp_dir(field,tile,bsub)
    
    
    
    root='%s_%s.%s' % (field,tile,image)  
    
    filelist='%s.in' % root
    name='%s/%s' % (xdir,filelist)
    
    # print('Filelist ',name)
    
    f=open(name,'w')
    for one in xx:
        if bsub:
            xname='%s2/%s/%s/%s'% (PREPDIR,field,tile,one['Filename'])
        else:
            xname='%s/%s/%s/%s'% (PREPDIR,field,tile,one['Filename'])

        if os.path.isfile(xname):
            f.write('%s\n' % xname)
        else:
            print('SwarpSetup: Warning: %s expected but not found' % xname)
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
    print('###SwarpSetup:  Finished swarp inputs for %s tile %s and image %s\n' % (field,tile,image))
    return   
    
    
def create_commmands_for_one_tile(field='LMC_c42',tile='T07',defaults=xdefault,bsub=False,use_config=True):
    '''
    Generate the standard set of commands for on tile, based on all of the separate
    exposure times that exist for the cell

    230619 - Added code to try to handle the situation where inputs are 
    somewhat inconsistent, assuming that if the tile is given as _b, one 
    actually wants to use the background subtracted data

    '''
    xtab=get_sum_tab(field,tile)

    if len(xtab)==0:
        print('SwarpSetup: Error: SwarpSetup.py: Nothing to swarp')
        return

    # Homoogenize the inputs in a situation where _b has been added to the file name
    if tile.count('_b'):
        tile=tile.replace('_b','')
        bsub=True


    print('xtab.info')
    xtab.info()
    images=np.unique(xtab['Image'])

    for one_image in images:
        create_swarp_command(field=field,tile=tile,image=one_image,defaults=defaults,bsub=bsub,use_config=use_config)
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
    use_config=True 

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-bsub':
            bsub=True
        elif argv[i]=='-ave_pos':
            use_config=False
        elif argv[i][0]=='-':
            print('SwarpSetup: Error: Unknwon switch  %s' % argv[i])
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
        print('Knox ',xdirs)





    open_log('%s.log' % field,reinitialize=False)
    for one in tiles:
        create_commmands_for_one_tile(field,one,bsub=bsub,use_config=use_config)
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


