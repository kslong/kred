#!/usr/bin/env python
# coding: utf-8

# # Create Combined images of the different filters in the field
# 
# This can only be run after PrepFiles and SumFiles have benn run.    The routines here generate 
# inputs to run swarp which combines the individual CCD images into tile images.  It does not
# actually run swarp at present



import os, stat
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import vstack


def get_sum_tab(field='LMC_c42',tile=7):
    
    xname='../workspace/Summary/%s_%d_imsum.txt' % (field,tile)
    xtab=ascii.read(xname)
    return xtab



def summarize(field='LMC_c42',tile=7):
    x=get_sum_tab(field,tile)

    ra=np.average(x['RA'])
    dec=np.average(x['Dec'])

    print('The center of thist tile is %.5f  %.5f' % (ra,dec))


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
    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1

    print('\nSII')
    times,counts=np.unique(s2['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1
        
    print('\nR')
    times,counts=np.unique(r['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1    
        
    print('\nN708')
    times,counts=np.unique(n708['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        print('   exp %8s  no %4d' % (times[i],counts[i]))
        i+=1
    return ra,dec



def create_swarp_dir(field='LMC_c42',tile=7):
    '''
    Create a diretory for the swarp outputs if it does not exist
    '''
    outdir='../workspace/%s_d/%d/' % (field,tile)
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

SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
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





def create_swarp_command(field='LMC_c42',tile=7,filt='Ha',exp=[30,800],defaults=xdefault):
    '''
    Generate the inputs necessary to run swarp
    '''
    print('\n### Creating swarp inputs for %s tile %d and filter %s' % (field,tile,filt))
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
    
    root='%s_%d.%s' % (field,tile,filt)  
    root+=exp_string
    # print('root ',root)    
    
    filelist='%s.in' % root
    name='%s/%s' % (xdir,filelist)
    
    # print('Filelist ',name)
    
    f=open(name,'w')
    for one in xxxx:
        xname='../../%s/%d/%s'% (field,tile,one['Filename'])
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
    print('### Finished swarp inputs for %s tile %d and filter %s\n' % (field,tile,filt))
    return   
    
    



def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
        elif field=='':
            field=argv[i]
        else:
            tiles.append(int(argv[i]))
        i+=1

    if xall:
        prep_all_tiles(field,redo)
    else:
        if len(tiles)==0:
            print('The tiles to be processed must be listed after the field, unless -all is invoked')
        for one in tiles:
            prep_one_tile(field,one,redo)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


