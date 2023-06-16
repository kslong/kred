#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:

Prepare files for Combining with Swpart
 
This is where any processing of the individal files is done that cannot be done by swarp.

Prior to running PrepFiles one must have 

* Downladed the MEF an put them in the standard directory structure
* Run MefSum on the apprate fields to create tables that contain the RA's and DECs of the corrners of CCDs

Run this routine to rescale the images to the same maagnitude scale and optionally subtract
background

 

Usage:  PrepMef.py [-h] [-finish] [-sub_back] [-np 4] Field_name or names

where -h   --- to print this documentation and exit
       -finish --> do not redo files that have already been processed
      -sub_back --> causes background to be estimated from individula images
      -np 4   --. To run with multiple threads



Description:

    If the preliminaries indicated above have taken place, PrepMef   will create directories if 
    necessary to store the results.  The program then reads tables files in the Summary directory
    where there is a table that identify the Mef files and extenstions associated with them.
    it then processes the each Mefile individual.

    The outputfiles are stored in directories named data, that have specific place in the
    file structure. 

    So all of the processed files for LMC_c42 will be stored in DECam_PREP/LMC_c42/data

    With this approach, one must create links to the apprtate files next

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
SUMDIR='Summary/'



def prep_one_det(filename='DECam_MEF/LMC_c42/mef/c4d_211111_024404_ooi_N673_v1.fits.fz',ext='S2',
                     prepdir='DECam_PREP/LMC_c42/T07/',back=0, redo=True):
    '''
    Prepare one file for combining with swarp.  This version scales by MAGZERO
    where 1DN will be a flux corresponding to mag 28
    History:

    230504 - Removed .fz from output file names.  Note that this does not mean that the
    data is not compressed.   Added extra keywords to primary header of output file
    that may be useful for swarp, including the rescaled saturation level
    
    230611 - Addapted fro prep_one_file.  This version allows for several backgroud
    options, based on information assembled by MEFSum

    :
    '''


    # find the extension we want given thee file name

    words=filename.split('/')
    new_file=words[-1]
    outfile='%s/%s' %(prepdir,new_file)
    outfile=outfile.replace('.fits.fz','_%s.fits' % (ext))
    # print(outfile)

    if redo==False and os.path.isfile(outfile):
        return outfile

    f=fits.open(filename)

    hdu0=f[0].copy()
    hdu1=f[ext].copy()
    fout=fits.HDUList([hdu0,hdu1])


    factor=10**(0.4*(28.0 - fout[0].header['MAGZERO']))
    fout[1].data*=factor

    # mean,median,std=sigma_clipped_stats(f[1].data,sigma_lower=3,sigma_upper=2,grow=3)


    fout[0].header['XFACTOR']=factor

    if back!=0:
        fout[0].header['SUB_BACK']='TRUE'
        fout[1].data-=back
    else:
        fout[0].header['SUB_BACK']='FALSE'

    # fout[0].header['XMED']=median

    xgain=0.5*(f[1].header['GAINA']+f[1].header['GAINB'])

    fout[0].header['GAIN']=xgain*factor

    xsat=np.minimum(f[1].header['SATURATA'],f[1].header['SATURATB'])

    fout[0].header['SATURAT']=xsat*factor


    xread=0.5*(fout[1].header['RDNOISEA']+fout[1].header['RDNOISEB'])

    fout[0].header['RDNOISE']=xread*factor


    fout.writeto(outfile,overwrite=True)
    return outfile



def prep_one_mef(field='LMC_c42',root='c4d_190109_061931_ooi_N662_v1',back_type='min',redo=False,outdir=''):
    '''
    Split the mef files into their individual extenstions and prepare them for the downstream parts of the processing
    
    This version stores all of the actual output files in a single directory for each field. Setup_Tiles (will be modified
    for the new file structure
    
    '''



    if redo==False:
        print('Beginning    %s in field %s with back_type %s (No reprocessing)' % (root,field,back_type))
    else:
        print('Beginning    %s in field %s with back_type %s (Reprocessins is on)' % (root,field,back_type))
    
    mef_file='%s/%s_mef.tab' % (SUMDIR,field)
    try:
        
        mef=ascii.read(mef_file)
        mef=mef[mef['Root']==root]
    except:
        print('Error: Could not locate %s' % (mef_file))
        return
    
    det_file='%s/%s_det.tab' % (SUMDIR,field)
    try:
        det=ascii.read(det_file)
        det=det[det['Root']==root]
    except:
        print('Error: Could not locate %s' % (det_file))
        return
    
    time_start=timeit.default_timer()
    
    # print(len(det),len(mef))
    
    if back_type=='none':
        back=0
    elif back_type=='med':
        back=np.median(det['Med'])
    elif back_type=='min':
        back=np.min(det['Med'])
    else:
        print('Error: Unknown option for background subtraction %s' % (back_type))
        return
    
    # print(back)
    
    if outdir=='':
        outdir='%s/%s/data/' % (PREPDIR,field)
    if os.path.isdir(outdir)==False:
        print('Creating Prep Dir as :',outdir)
        os.makedirs(outdir)
    
    mef_fits_file='DECam_MEF/%s/mef/%s.fits.fz' % (field,root)
    
    ndone=0
    for one in det:
        outfile='%s/%s_%s.fits' % (outdir,root,one['EXTNAME'])
        if os.path.isfile(outfile) and redo==False :
            pass
        else:
            prep_one_det(mef_fits_file,one['EXTNAME'],outdir)
            ndone+=1
    
    elapsed = timeit.default_timer() - time_start
    print('Finished mef %s in %s - processed %d ext in %.1f s'  % (root,field,ndone,elapsed))
    
    return

    
    


def prep_one_field(field='LMC_c42',back_type='none',redo=False,outdir=''):
    '''
    Prep a field with a single processor

    This routine identifies the mefs to be processed when one
    is in single processor mode
    '''

    tab_name='Summary/%s_mef.tab' % (field)
    try:
        xtab=ascii.read(tab_name)
    except:
        print('Failed to read: %s' % tab_name)
        return

    xtab=xtab[xtab['Field']==field]

    
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

        prep_one_mef(field=field,root=one['Root'],back_type=back_type,redo=redo,outdir=outdir)
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


def xprep_one_field(field='LMC_c42',back_type='none', redo=False,outdir='',nproc=4):
    '''
    Process the mef files in a field usingwith using more than one core.

    This routine creates a bundle of files to process 
    at once. This is necessary because multiprocessing
    opens many files (note that in linux everything is
    a file, so this does not refer to fits files), and
    this can exceed the file limit.  The sympton of this
    is a message that says too many files are open.  If
    one sees this message, one may need to reduce n
    below.

    '''



    tab_name='Summary/%s_mef.tab' % (field)
    try:
        xtab=ascii.read(tab_name)
    except:
        print('Failed to read: %s' % tab_name)
        return

    xtab=xtab[xtab['Field']==field]

    
    start_time = timeit.default_timer()

    jobs=[]
    for one in xtab:
        p=multiprocessing.Process(target=prep_one_mef,args=[field,one['Root'],back_type,redo,outdir])
        jobs.append(p)


    i=0
    while i<nproc and i<len(jobs):
        t = time.localtime()
        one=jobs[i]
        one.start()
        i+=1


    njobs=get_no_jobs(jobs)

    while i<len(jobs):
        time.sleep(2)
        njobs=get_no_jobs(jobs)

        while njobs<nproc and i<len(jobs):
            t = time.localtime()
            one=jobs[i]
            one.start()
            njobs+=1
            i+=1

    p.join()
    p.close()
        
    elapsed = timeit.default_timer() - start_time
    print('Completed multiprocessing of field %s  ' % (field))

    return



        

    
def steer(argv):
    '''
    This is just a steering routine
    '''
    fields=[]
    xall=False
    redo=True 
    nproc=1
    xback='none' 
    outdir=''

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
            xback='med'
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
        else:
            fields.append(argv[i])
        i+=1




    for one_field in fields:
        if nproc<2:
             prep_one_field(one_field,xback,redo)
        else:
            print('Processing in parallel with %d cores' % (nproc))
            xprep_one_field(one_field,xback,redo,outdir,nproc)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__) 



