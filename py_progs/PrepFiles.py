#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:

Prepare files for Combining into tile images
 
This is where any processing of the individal files is done that cannot be done by swarp.

Prior to running PrepFiles one must have 

* Downladed the MEF an put them in the standard directory structure
* Run MefSum on the apprate fields to create tables that contain the RA's and DECs of the corrners of CCDs
* Run SetupTile to define which images and extensions are contained in specific tiles


 

Usage:  PrepFiles.py [-h] [-finish] [-all] [-sub_back] [-np 4] Field_name  [T01 T02 ... one or more tile numbers]

where -h   --- to print this documentation and exit
       -finish --> do not redo files that have already been processed
      -all    --> do all tiles for a fields
      -sub_back --> causes background to be estimated from individula images
      -np 4   --. To run with multiple threads

If -all is not provide then there must be one or more tile numbes, e. g. T01 T05 T07, listed


Description:

    If the preliminaries indicated above have taken place, PrepFiles will create directories if 
    necessary to store the results.  The program then reads tables files in the Summary directory
    where there is a table that identify the files and extension associated with a tile, and
    processes these files individually.

Primary routines:


Notes:

    As written one can only PrepFiles from a specific MEF directory


History:

230513 ksl Coding begun

'''




import os
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from glob import glob
import numpy as np
import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)


MEFDIR='DECam_MEF/'
PREPDIR='DECam_PREP/'




def setup_one_tile(field='LMC_c42',tile='T07'):
    # First check that we have data 
    
    mefdir='%s/%s/mef/' % (MEFDIR,field)

    if os.path.isdir(mefdir) == False:
        print('Error: %s does not appear to exist' % mefdir)
        raise IOError
    
    xfiles='%s/*ooi*fits.fz' % (mefdir)
    files=glob(xfiles)
    # print('There are %d files in the mef directory' % (len(files)))
    if len(files)==0:
        print('Could not locate any files in the mef directory')
        raise IOError
        
    # Now create diretories to receive the processed files if that is necessary
    
    prepdir='%s/%s/%s' % (PREPDIR,field,tile)
    
    if os.path.isdir(prepdir)==False:
        print('Creating the Prep directory as %s' % prepdir)
        os.makedirs(prepdir)
    else:
        print('The Prep directory %s already exists' % prepdir )
        
    return  mefdir,prepdir


def prep_one_file(filename='DECam_MEF/LMC_c42/mef/c4d_211111_024404_ooi_N673_v1.fits.fz',ext='S2',
                     prepdir='DECam_PREP/LMC_c42/T07/',subtract_back=False, redo=True):
    '''
    Prepare one file for combining with swarp.  This version scales by MAGZERO
    where 1DN will be a flux corresponding to mag 28
    History:

    230504 - Removed .fz from output file names.  Note that this does not mean that the
    data is not compressed.   Added extra keywords to primary header of output file
    that may be useful for swarp, including the rescaled saturation level

    :
    '''

 
    # find the extension we want given thee file name
    
    words=filename.split('/')
    new_file=words[-1]
    outfile='%s/%s' %(prepdir,new_file)
    outfile=outfile.replace('.fits.fz','_%s.fits' % (ext))
    print(outfile)
    
    if redo==False and os.path.isfile(outfile):
        return outfile
    
    f=fits.open(filename)    

    hdu0=f[0].copy()
    hdu1=f[ext].copy()
    fout=fits.HDUList([hdu0,hdu1])

    
    factor=10**(0.4*(28.0 - fout[0].header['MAGZERO']))
    fout[1].data*=factor

    mean,median,std=sigma_clipped_stats(f[1].data,sigma_lower=3,sigma_upper=2,grow=3)

    if subtract_back:
        fout[1].data-=median

    fout[0].header['XFACTOR']=factor

    if subtract_back:
        fout[0].header['SUB_BACK']='TRUE'
    else:
        fout[0].header['SUB_BACK']='FALSE'

    fout[0].header['XMED']=median

    xgain=0.5*(f[1].header['GAINA']+f[1].header['GAINB'])

    fout[0].header['GAIN']=xgain*factor

    xsat=np.minimum(f[1].header['SATURATA'],f[1].header['SATURATB'])

    fout[0].header['SATURAT']=xsat*factor


    xread=0.5*(fout[1].header['RDNOISEA']+fout[1].header['RDNOISEB'])

    fout[0].header['RDNOISE']=xread*factor

    
    fout.writeto(outfile,overwrite=True)
    return outfile
    


def prep_one_tile(field='LMC_c42',tile='T07',sub_back=False,redo=False):
    '''
    Prep a tile with a single processor
    '''

    tab_name='Summary/%s_%s.txt' % (field,tile)
    try:
        xtab=ascii.read(tab_name)
    except:
        print('Failed to read: %s' % tab_name)
        return


    try:
        mefdir,prepdir=setup_one_tile(field,tile)
    except:
        print('Failed to set up tiles')
        return

    
    xfiles=[]
    for one in xtab:
        xfiles.append('%s/%s.fits.fz' % (mefdir,one['Root']))

    xtab['Filename']=xfiles


    
    start_time = timeit.default_timer()
    
    nn=len(xtab)
    n=50
    if nn>200:
        n=100

    i=0
    for one in xtab:
        if i % n == 0:
            elapsed = timeit.default_timer() - start_time
            print('Completed %4d of %4d files in %f s' % (i,nn,elapsed))

        xfile=prep_one_file(one['Filename'],one['EXTNAME'],prepdir,sub_back,redo)
        i+=1
        
    elapsed = timeit.default_timer() - start_time
    print('All done in %f s' % elapsed)    
    return

def get_no_jobs(jobs):
    '''
    Check how many jobs are running
    '''
    njobs=0
    for one in jobs:
        if one.is_alive():
            njobs+=1
    return njobs


def xprep_one_tile(field='LMC_c42',tile='T07',sub_back=False, redo=False,nproc=4):
    '''
    Process tiles with using more than one core.

    This routine creates a bundle of files to process 
    at once. This is necessary because multiprocessing
    opens many files (note that in linux everything is
    a file, so this does not refer to fits files), and
    this can exceed the file limit.  The sympton of this
    is a message that says too many files are open.  If
    one sees this message, one may need to reduce n
    below.

    '''





    tab_name='Summary/%s_%s.txt' % (field,tile)
    try:
        xtab=ascii.read(tab_name)
    except:
        print('Failed to read: %s' % tab_name)
        return


    try:
        mefdir,prepdir=setup_one_tile(field,tile)
    except:
        print('Failed to set up tiles')
        return

    
    xfiles=[]
    for one in xtab:
        xfiles.append('%s/%s.fits.fz' % (mefdir,one['Root']))

    xtab['Filename']=xfiles


    
    start_time = timeit.default_timer()
    
    nn=len(xtab)
    nncount=0

    n=200  # This is the maximum number of jobs that will be run in a Cycle

    bounds=[]
    k=0
    while k<nn:
        bounds.append(k)
        k+=100
    if bounds[-1]<nn:
        bounds.append(nn)


    k=1
    while k<len(bounds):
        files=xtab['Filename'][bounds[k-1]:bounds[k]]
        exten=xtab['EXTNAME'][bounds[k-1]:bounds[k]]

        print(' Started processing of files %4d - %4d' % (bounds[k-1],bounds[k]))

        jobs=[]
        i=0
        for one in files:
            p=multiprocessing.Process(target=prep_one_file,args=(one,exten[i],prepdir,sub_back,redo))
            jobs.append(p)
            i+=1


        i=0
        while i<nproc and i<len(jobs):
            t = time.localtime()
            one=jobs[i]
            one.start()
            nncount+=1
            i+=1


        njobs=get_no_jobs(jobs)

        while i<len(jobs):
            time.sleep(2)
            njobs=get_no_jobs(jobs)

            while njobs<nproc and i<len(jobs):
                t = time.localtime()
                one=jobs[i]
                one.start()
                nncount+=1
                njobs+=1
                i+=1

        p.join()
        p.close()
        
        elapsed = timeit.default_timer() - start_time
        print('Finished processing of files %4d - %4d of %4d total at elapsed time %8.1f' % (bounds[k-1],bounds[k],nn, elapsed))
        k+=1

    print('Completed multiprocessing of field %s  tile %s ' % (field,tile))



        

    
def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True 
    nproc=1
    xback=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i]=='sub_back':
            xback=True
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
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

    if len(tiles)==0:
        print('The tiles to be processed must be listed after the field, unless -all is invoked')
        return

    for one in tiles:
        if nproc<2:
             prep_one_tile(field,one,xback,redo)
        else:
            print('Processing in parallel with %d cores' % (nproc))
            xprep_one_tile(field,one,xback,redo,nproc)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__) 




