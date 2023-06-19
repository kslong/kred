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

Usage:   Swarp.py [-all] [-bsub] field [tiles]

where -all will cause swarp to be run on all 16 tiles.

and   -bsub will cause swarp to be run on the files in background subtracted _b directories

If one wants to run only 1 or a few tiles then the command will be something like

Swarp.py LMC_c42  1 3 5


'''


import os, stat
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from glob import glob
import timeit
import subprocess

exec_dir=os.getcwd()


SWARPDIR='DECam_SWARP'
PREPDIR=os.path.abspath('DECam_PREP')


def get_sum_tab(field='LMC_c42',tile='T07'):
    '''
    Read a table that summarizes information about
    the various files athat are used for a tile
    
    The location of the table files is currently hardocaded
    '''
    
    xname='Summary/%s_%s_imsum.txt' % (field,tile)
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

    ra=np.average(x['RA'])
    dec=np.average(x['Dec'])

    print('The center of this tile is %.5f  %.5f' % (ra,dec))


    ha=x[x['Filter']=='Ha']
    s2=x[x['Filter']=='SII']
    r=x[x['Filter']=='R']
    n708=x[x['Filter']=='N708']    
    print('Ha  images  : %3d' % len(ha))
    print('SII images  : %3d' % len(s2))
    print('R   images  : %3d' % len(r))
    print('N708 images : %3d' % len(n708))
    
    print('\nHa')
    times,counts=np.unique(ha['Exptime'],return_counts=True)

    records=[times,counts]
    xtab=Table(records,names=['Exptime','Number'])
    xtab['Filter']='Ha'
    ztab=xtab.copy()

    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1

    print('\nSII')
    times,counts=np.unique(s2['Exptime'],return_counts=True)

    records=[times,counts]
    xtab=Table(records,names=['Exptime','Number'])
    xtab['Filter']='SII'
    ztab=vstack([ztab,xtab])

    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1
        
    print('\nR')
    times,counts=np.unique(r['Exptime'],return_counts=True)

    records=[times,counts]
    xtab=Table(records,names=['Exptime','Number'])
    xtab['Filter']='R'
    ztab=vstack([ztab,xtab])

    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1    
        
    print('\nN708')
    times,counts=np.unique(n708['Exptime'],return_counts=True)

    records=[times,counts]
    xtab=Table(records,names=['Exptime','Number'])
    xtab['Filter']='N708'
    ztab=vstack([ztab,xtab])


    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1

    print(ztab)
    return ra,dec



def create_swarp_dir(field='LMC_c42',tile=7):
    '''
    Create a diretory for the swarp outputs if it does not exist
    '''
    outdir='%s/%s/%s/' % (SWARPDIR,field,tile)
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





def create_swarp_command(field='LMC_c42',tile='T07',filt='Ha',exp=[800],defaults=xdefault):
    '''
    Generate the inputs necessary to run swarp where exp corresponds to one
    or more exposure times for a particular filter and tile.  
    '''
    print('\n### Creating swarp inputs for %s tile %s and filter %s' % (field,tile,filt))
    x=get_sum_tab(field,tile)
    xx=x[x['Filter']==filt]
    if len(xx)==0:
        print('There are no observations with filter %s')
        return
    
    if type (exp)== int:
        exp=[exp]
    
    ra=np.average(x['RA'])
    dec=np.average(x['Dec'])

    i=0
    while i<len(exp):
        xxx=xx[xx['Exptime']==exp[i]]
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
        
    xdir=create_swarp_dir(field,tile)
    
    # print('xdir ',xdir)
    

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
        xname='%s/%s/%s/%s'% (PREPDIR,field,tile,one['Filename'])
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
    f.write('swarp @%s -c %s\n' % (filelist,default_name))
    f.close()
    
    os.chmod('%s/%s.run' % (xdir,root),stat.S_IRWXU)
    print('### Finished swarp inputs for %s tile %s and filter %s\n' % (field,tile,filt))
    return   
    
    



def run_swarp(field='LMC_c42',tile='T07'):
    '''
    run_all of the swarp commands in a particular field and tile

    The normal outputs from swarp are writeen to a .txt file,
    a few are written to the screen here so that one can see 
    that the routines have run correctly.  
    '''
    xstart=qstart=start_time=timeit.default_timer()
    run_dir='%s/%s/%s/' % (SWARPDIR,field,tile)

    try:
        os.chdir(run_dir)
    except:
        print('Error: Could not cd to %s' % run_dir)
        os.chdir(exec_dir)

    run_files=glob('*.run')
    # print(run_files)
    nfiles=len(run_files)
    n=1
    for one in run_files:
        print('\n***Beginning %s (%d of %d)' % (one,n,nfiles))
        outfile=one.replace('.run','.txt')
        command=one
        xout=subprocess.run(command,shell=True,capture_output=True)
        current_time=timeit.default_timer()

        print('***Writing last portion of swarp outputs\n')

        z=open(outfile,'w')
        zz=xout.stderr.decode()
        lines=zz.split('\n')
        xlines=[]
        for line in lines:
            if line.count('\x1b[1M>')==0:
                xlines.append(line)
                z.write('%s\n' % line)
        z.close()


        for xline in xlines[-10:]:
            print(xline)

        print('***Finished writing end of outputs')

        print('***Finished %s (%d of %d) in %.1f s\n' % (one,n,nfiles,current_time-start_time))
        start_time=current_time


        n+=1


    xcurrent_time=timeit.default_timer()
    print('\n***Completely done for tile %s in %.1f s\n' % (tile,xcurrent_time-xstart))


    os.chdir(exec_dir)
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
        elif argv[i]=='-finish':
            redo=False
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
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

    for one in tiles:
        if bsub and one.count('_b')==0:
            one='%s_b' % one
        run_swarp(field,one)


    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


