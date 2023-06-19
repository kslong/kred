#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create and run swarp on to create images
that can be used for x-matching back
background 


Command line usage (if any):

    usage: BackPrep.py [-all ] [-run] field tile

Description:  

    BackPrep.py without the -run simply sets up the Swarps

    With -run it will actaul run all the swarps

                                       
History:

230504 ksl Coding begun

'''




import os, stat

from glob import glob
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats

SWARPDIR='DECam_BACK'
PREPDIR=os.path.abspath('DECam_PREP')

from glob import glob
import timeit
import subprocess

exec_dir=os.getcwd()

def create_swarp_dir(field='LMC_c42',tile='T07'):
    '''
    Create a diretory for the swarp outputs if it does not exist
    '''
    outdir='%s/%s/%s/' % (SWARPDIR,field,tile)
    if os.path.isdir(outdir)==False:
        os.makedirs(outdir)
    return outdir



default='''
# Default configuration file for SWarp 2.38.0
# EB 2019-07-10
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME         foo.fits      # Output filename
WEIGHTOUT_NAME      foo.weight.fits # Output weight-map filename

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
CENTER                %.4f, %.4f       # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            2            # Pixel scale
IMAGE_SIZE          2000,2000          # Image size (0 = AUTOMATIC)

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







def prep_one_tile(field='LMC_c42',tile='T07'):
    try:
        sumfile='Summary/%s_%s.txt' % (field,tile)
        x=ascii.read(sumfile)
    except:
        print('Could not find %s' % sumfile)
        raise IOERROR


    
    indir='%s/%s/%s' % (PREPDIR,field,tile)
    if os.path.exists(indir)==False:
        print('Could not find the directory' % indir)
        raise IOERROR



    # Now create a directory if it does not exist

    outdir='%s/%s/%s' % (SWARPDIR,field,tile)
    if os.path.exists(outdir)==False:
        os.makedirs(outdir)
        print('The output directory did not exist, so it is being created')
    else:
        print('The output directory exists')

    
    ra=np.average(x['CENRA1'])
    dec=np.average(x['CENDEC1'])

    # This needs to be relative to the outdir

    indir='%s/%s/%s/' % (PREPDIR,field,tile)

    f=open('%s/swarp.default' % outdir,'w')
    f.write(default % (ra,dec))
    f.close()

    runfile='%s/RunSwarp' %outdir
    f=open(runfile,'w')
    for one in x:
        root=one['Filename'].replace('.fits','')
        f.write('swarp %s/%s.fits -IMAGEOUT_NAME %s.fits  -WEIGHTOUT_NAME  %s.weight.fits -c swarp.default\n' % (indir,root,root,root))
        # gen_one_swarp_command(one['Filename'],ra,dec,indir,outdir)
    f.close()
    os.chmod(runfile,stat.S_IRWXU)



    return 


def run_back(field='LMC_c42',tile='T07'):
    '''
    run swarp in the for the background fields
    '''

    xstart=qstart=start_time=timeit.default_timer()
    run_dir='%s/%s/%s/' % (SWARPDIR,field,tile)

    try:
        os.chdir(run_dir)
    except:
        print('Error: Could not cd to %s' % run_dir)
        os.chdir(exec_dir)

    command='RunSwarp'
    outfile='Swarp.out.txt'
    xout=subprocess.run(command,shell=True,capture_output=True)

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


    current_time=timeit.default_timer()
    print('\n***Completely done for tile %s in %.1f s\n' % (tile,current_time-xstart))


    os.chdir(exec_dir)
    return



                              
                              
def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    xrun=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-run':
            xrun=True
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s' % argv[i])
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

    if len(tiles)==0:
        print('The tiles to be processed must be listed after the field, unless -all is invoked')
    for one in tiles:
        prep_one_tile(field,one)

    if xrun:
        for one in tiles:
            print('Starting to run %s %s' % (field,one))
            run_back(field,one)


    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


