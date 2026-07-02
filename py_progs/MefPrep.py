#!/usr/bin/env python
# coding: utf-8


"""Prepare files for combining with Swarp

Space Telescope Science Institute

Command Line Usage
------------------

::

    MefPrep.py [-h] [-all] [-finish] [-back_none] [-back_min] [-np N]
               [-zp COLUMN] FIELD [FIELD ...]

Rescales individual DECam CCD images to a common flux scale (1 DN = mag 28)
and optionally subtracts an initial background estimate.  Output files are
written to ``DECam_PREP/{field}/data/``.

Prerequisites: MEF files downloaded into the standard directory structure and
``MefSum.py`` run to create the ``Summary/{field}_mef.tab`` inventory tables.

Optional Arguments
------------------

-h
    Print this help and exit.

-all
    Process all fields found under ``DECam_MEF/`` (interactive confirmation
    required).

-finish
    Skip output files that already exist (do not reprocess).

-back_none
    Do not subtract any background estimate.

-back_min
    Subtract the minimum CCD background across all detectors in each exposure
    (instead of the default mode estimate).

-np N
    Number of parallel worker processes (default: 1).

-zp COLUMN
    Use empirical zero points from the named column in
    ``Summary/{field}_mef.tab`` (e.g. ``-zp ZP_smash_r``) instead of the
    ``MAGZERO`` header keyword.  The value must satisfy ``20 < ZP < 35``; if
    the column is absent or the value is invalid, the script falls back to the
    header ``MAGZERO`` with a warning.  Generates ``ZP_USE`` and
    ``ZP_SRC`` keywords in every output header.

Output Header Keywords
----------------------

* ``ZP_USE`` -- zero point actually applied for flux scaling
* ``ZP_SRC`` -- source of that zero point (column name or ``MAGZERO``)
* ``XFACTOR``     -- multiplicative scale factor applied to pixel data
* ``SUB_BACK``    -- ``TRUE`` / ``FALSE`` depending on whether background was subtracted
* ``GAIN``, ``SATURAT``, ``RDNOISE`` -- rescaled detector properties

Examples
--------

Standard single-field run with background subtraction (mode)::

    MefPrep.py -np 8 LMC_c42

Re-run using empirical SMASH r-band zero points::

    MefPrep.py -np 8 -zp ZP_smash_r LMC_c42

Skip already-processed files, no background subtraction::

    MefPrep.py -finish -back_none LMC_c42

Notes
-----

Unlike ``MefSum.py``, individual files are processed in separate threads so
all requested cores are typically utilised.  If too many files are open
simultaneously, reduce ``-np``.

Version History
---------------

230513 ksl
    Coding begun

230621 ksl
    Revised so that the default is to subtract background

260608 ksl
    Add -zp COLUMN to use empirical zero points from Summary/_mef.tab instead
    of the header MAGZERO.  ZP_USE and ZP_SRC written to output headers.
    Fixed pre-existing bug where computed back value was not passed to prep_one_det.

"""





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
from log import *


MEFDIR='DECam_MEF/'
CCDDIR='DECam_CCD/'
SUMDIR='Summary/'



def prep_one_det(filename='DECam_MEF/LMC_c42/mef/c4d_211111_024404_ooi_N673_v1.fits.fz',ext='S2',
                     prepdir='DECam_PREP/LMC_c42/T07/',back=0, redo=True, magzero=None, zp_col=''):
    '''
    Prepare one file for combining with swarp.  This version scales by MAGZERO
    where 1DN will be a flux corresponding to mag 28
    History:

    230504 - Removed .fz from output file names.  Note that this does not mean that the
    data is not compressed.   Added extra keywords to primary header of output file
    that may be useful for swarp, including the rescaled saturation level

    230611 - Addapted from prep_one_file.  This version allows for several backgroud
    options, based on information assembled by MEFSum

    :
    '''


    # find the extension we want given thee file name

    words=filename.split('/')
    new_file=words[-1]
    outfile='%s/%s' %(prepdir,new_file)
    outfile=outfile.replace('.fits.fz','_%s.fits' % (ext))

    if redo==False and os.path.isfile(outfile):
        return outfile

    f=fits.open(filename)

    hdu0=f[0].copy()
    hdu1=f[ext].copy()
    fout=fits.HDUList([hdu0,hdu1])

    if magzero is not None:
        zp = magzero
        zp_src = zp_col if zp_col else 'table'
    else:
        zp = fout[0].header['MAGZERO']
        zp_src = 'MAGZERO'

    factor=10**(0.4*(28.0 - zp))
    fout[1].data*=factor

    fout[0].header['XFACTOR'] = factor
    fout[0].header['ZP_USE'] = (zp,     'Zero point used for flux scaling')
    fout[0].header['ZP_SRC'] = (zp_src, 'ZP source: column name or MAGZERO')

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



def prep_one_mef(field='LMC_c42',root='c4d_190109_061931_ooi_N662_v1',back_type='min',redo=False,outdir='',zp_col=''):
    '''
    Split the mef files into their individual extenstions and prepare them for the downstream parts of the processing
    
    This version stores all of the actual output files in a single directory for each field. Setup_Tiles (will be modified
    for the new file structure


    Note that background subtraction here, if is carriout out, subtractis the same value from all of 
    the CCDs in a single image.  

    What is done is slightly confusing.  MEFSum.py has estimated a different background in each of
    the CCDs in the exposure and stored this in the det file.  Here we take all of those measurements
    for the mef image and either take the median of those values or the minimum.

    
    '''



    if redo==False:
        print('MefPrep: Beginning    %s in field %s with back_type %s (No reprocessing)' % (root,field,back_type))
    else:
        print('MefPrep: Beginning    %s in field %s with back_type %s (Reprocessins is on)' % (root,field,back_type))
    
    mef_file='%s/%s_mef.tab' % (SUMDIR,field)
    try:
        mef=ascii.read(mef_file)
        mef=mef[mef['Root']==root]
    except:
        print('MefPrep: Error: Could not locate %s' % (mef_file))
        return

    det_file='%s/%s_det.tab' % (SUMDIR,field)
    try:
        det=ascii.read(det_file)
        det=det[det['Root']==root]
    except:
        print('MefPrep: Error: Could not locate %s' % (det_file))
        return

    time_start=timeit.default_timer()

    if back_type=='none':
        back=0
    elif back_type=='mode':
        back=np.median(det['Mode'])
    elif back_type=='med':
        back=np.median(det['Med'])
    elif back_type=='min':
        back=np.min(det['Med'])
    else:
        print('MefPrep: Error: Unknown option for background subtraction %s' % (back_type))
        return

    # resolve zero point: use named column from _mef.tab if valid, else fall back to header
    magzero = None
    if zp_col:
        if zp_col in mef.colnames:
            val = float(mef[zp_col][0])
            if np.isfinite(val) and 20.0 < val < 35.0:
                magzero = val
            else:
                print('MefPrep: Warning: %s=%s for %s is invalid, falling back to header MAGZERO' % (zp_col, val, root))
        else:
            print('MefPrep: Warning: column %s not in %s, falling back to header MAGZERO' % (zp_col, mef_file))

    if outdir=='':
        outdir='%s/%s/data/' % (CCDDIR,field)
    if os.path.isdir(outdir)==False:
        print('MefPrep: Creating Prep Dir as :',outdir)
        os.makedirs(outdir,exist_ok=True)

    mef_fits_file='DECam_MEF/%s/mef/%s.fits.fz' % (field,root)

    ndone=0
    for one in det:
        outfile='%s/%s_%s.fits' % (outdir,root,one['EXTNAME'])
        if os.path.isfile(outfile) and redo==False :
            pass
        else:
            prep_one_det(mef_fits_file, one['EXTNAME'], outdir, back=back, redo=redo,
                         magzero=magzero, zp_col=zp_col)
            ndone+=1
    
    elapsed = timeit.default_timer() - time_start
    print('MefPrep: Finished mef %s in %s - processed %d ext in %.1f s'  % (root,field,ndone,elapsed))
    
    return

    
    


def prep_one_field(field='LMC_c42',back_type='none',redo=False,outdir='',zp_col=''):
    '''
    Prep a field with a single processor

    This routine identifies the mefs to be processed when one
    is in single processor mode
    '''

    tab_name='Summary/%s_mef.tab' % (field)
    try:
        xtab=ascii.read(tab_name)
    except:
        print('MefPrep: Failed to read: %s' % tab_name)
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
            print('MefPrep: Completed %4d of %4d files in %f s' % (i,nn,elapsed))

        prep_one_mef(field=field,root=one['Root'],back_type=back_type,redo=redo,outdir=outdir,zp_col=zp_col)
        i+=1
        
    elapsed = timeit.default_timer() - start_time
    print('MefPrep: All done in %f s' % elapsed)    
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


def xprep_one_field(field='LMC_c42',back_type='none', redo=False,outdir='',nproc=4,zp_col=''):
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
        p=multiprocessing.Process(target=prep_one_mef,args=[field,one['Root'],back_type,redo,outdir],
                                   kwargs={'zp_col': zp_col})
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
    xback='mode'
    outdir=''
    zp_col=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i]=='-back_none':
            xback='none'
        elif argv[i]=='-back_min':
            xback='min'
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i]=='-zp':
            i+=1
            zp_col=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown switch  %s' % argv[i])
            return
        else:
            fields.append(argv[i])
        i+=1

    if xall==True:
        xfields=glob('%s/*' % MEFDIR)

        response = input("Did you really want to run MefSum on all (%d) fields (yes/no): " % (len(xfields)))

        # Convert the response to lowercase for case-insensitive comparison
        response = response.lower()

        if response.count('y'): 
            print("OK, beginning MefSum for all of the fields ")
        else:
            print("OK, it is easy to get confused on the inputs")
            return

        for one in xfields:
            words=one.split('/')
            fields.append(words[-1])


    xtime_start=timeit.default_timer()

    for one_field in fields:
        open_log('%s.log' % one_field,reinitialize=False)
        log_message('Starting MefPrep on %s with back=%s zp_col=%s' %(one_field,xback,zp_col if zp_col else 'MAGZERO'))

        if nproc<2:
             prep_one_field(one_field,xback,redo,outdir,zp_col)
        else:
            print('Processing in parallel with %d cores' % (nproc))
            xprep_one_field(one_field,xback,redo,outdir,nproc,zp_col)
        log_message('Finished MefPrep on %s with back=%s zp_col=%s' %(one_field,xback,zp_col if zp_col else 'MAGZERO'))
        close_log()

    elapsed = timeit.default_timer() - xtime_start
    print('This entire job ran to completion in %d s' % elapsed)

    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__) 




