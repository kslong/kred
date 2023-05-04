#!/usr/bin/env python
# coding: utf-8

'''
Prepare files for Combining into tile images

This extracts a relevent set of CCD images from the raw MEF files and writes them to a simpler version of the same files with be desired science data in the only image extension.  It also does background subraction etc if that is requred
# 
The names are the same as used by Photpipe
 
This the first step after the data has been prepared with Photopipe, and deposited in the raw directories
j 
This is where any processing of the individal files is done that cannot be done by swarp**
 
This versions uses zero points for rescaling, and does background subtrraction

Usage:  PrepFiles.py [-finish] [-all]  Field_name  [one or more tile numbers]

where -finish --> do not redo files that have already been processed
      -all    --> do all tiles for a fields

If -all is not provide then there must be one or more tile numbes, e. g. 1 5 7, listed
'''
 



import os
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from glob import glob
import numpy as np
import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)




def setup_one_tile(field='LMC_c42',tile=7):
    # First check that we have data 
    
    rawdir='../rawdata/%s/%d' % (field,tile)

    if os.path.isdir(rawdir) == False:
        print('Error: %s does not appear to exist' % rawdir)
        raise IOError
    
    xfiles='%s/%s*fits.fz' % (rawdir,field)
    files=glob(xfiles)
    print('There are %d files in the raw directory' % (len(files)))
    if len(files)==0:
        raise IOError
        
    # Now create diretories to receive the processed files if that is necessary
    
    workdir='../workspace/%s/%d/' % (field,tile)
    
    if os.path.isdir(workdir)==False:
        print('Creating the Working directory as %s' % workdir)
        os.makedirs(workdir)
    else:
        print('The working directory %s already exists' % workdir )
        
    return  rawdir,workdir,files


def prep_one_file(filename='../rawdata/LMC_c42/7/LMC_c42.211111.1050216_ooi_N673_v1_S17.fits.fz',
                     workdir='../workspace/LMC_c42/7/',redo=True):
    '''
    Prepare one file for combining with swarp.  This version scales by MAGZERO
    where 1DN will be a flux corresponding to mag 28
    '''

 
    # find the extensiton we want given thee file name
    
    words=filename.split('/')
    new_file=words[-1]
    outfile='%s/%s' %(workdir,new_file)
    
    if redo==False and os.path.isfile(outfile):
        return outfile
    
    f=fits.open(filename)    
    xname=words[-1]
    words=xname.split('.')
    words=words[2].split('_')
    ext=words[-1]
    
    names=[]
    i=1
    while i<len(f):
        if f[i].name!=ext:
            names.append(f[i].name)
        i+=1
    
    for name in names:
        f.pop(name)
        
    
    # factor=1./f[0].header['EXPTIME']
    factor=10**(0.4*(28.0 - f[0].header['MAGZERO']))
    f[1].data*=factor

    mean,median,std=sigma_clipped_stats(f[1].data,sigma_lower=3,sigma_upper=2,grow=3)

    f[1].data-=median

    f[0].header['XFACTOR']=factor
    f[0].header['XMED']=median
    
    f.writeto(outfile,overwrite=True)
    return outfile
    


def prep_one_tile(field='LMC_c42',tile=7,redo=False):
    try:
        rawdir,workdir,files=setup_one_tile(field,tile)
    except:
        print('Failed')
        return
    
    print(rawdir,workdir,files[0])
    
    start_time = timeit.default_timer()
    
    nn=len(files)
    n=50
    if nn>200:
        n=100

    i=0
    for one in files:
        if i % n == 0:
            elapsed = timeit.default_timer() - start_time
            print('Completed %4d of %4d files in %f s' % (i,nn,elapsed))

        xfile=prep_one_file(one,workdir,redo)
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


def xprep_one_tile(field='LMC_c42',tile=7,redo=False,np=4):
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

    try:
        rawdir,workdir,xfiles=setup_one_tile(field,tile)
    except:
        print('Failed')
        return
    
    
    start_time = timeit.default_timer()
    
    nn=len(xfiles)
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
        files=xfiles[bounds[k-1]:bounds[k]]

        print(' Started processing of files %4d - %4d' % (bounds[k-1],bounds[k]))

        jobs=[]
        for one in files:
            p=multiprocessing.Process(target=prep_one_file,args=(one,workdir,redo))
            jobs.append(p)


        i=0
        while i<np and i<len(jobs):
            t = time.localtime()
            one=jobs[i]
            one.start()
            nncount+=1
            i+=1


        njobs=get_no_jobs(jobs)

        while i<len(jobs):
            time.sleep(2)
            njobs=get_no_jobs(jobs)

            while njobs<np and i<len(jobs):
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

    print('Completed multiprocessing of field %s  tile %d ' % (field,tile))



        

def prep_all_tiles(field='LMC_c42',redo=False):
    '''
    Prepare all of the tiles in one set of data
    '''

    i=1
    while i<17:
        prep_one_tile(field,tile=i,redo=redo)
        i+=1
    
def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True 
    np=1

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i]=='-np':
            i+=1
            np=int(argv[i])
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
            if np<2:
                prep_one_tile(field,one,redo)
            else:
                print('Processing in parallel with %d cores' % (np))
                xprep_one_tile(field,one,redo,np)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__) 




