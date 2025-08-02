#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create a table (or tables) that contains information about 
how overlap regions of various images compare


Command line usage (if any):

    usage: BackStats.py  [-np 8] [-all] [-rm] [-sigma_clipped] field T01 T02 ...

    where field is a field name, and the following numbers are one 
    or more tiles

    -all implies to do all tiles

    -rm causes the files that were created by BackPrep to be deleted
    if the stats were carried out successfully.  This only occurs
    at the end of the program so that if something fails along
    the way, none of the backround files will have been deleted

    -np 8  implies to perform calculations in parallel, with in
    this case 8 threads

    -sigma_clipped is a diagnostic mode (which greatly increases
    the run time (but which calculates the mean, median, and std
    of the differences with sigma_clipped versions of these
    routines.)

Description:  



Primary routines:


Notes:

    This routine can take a while to run since individual tiles
    involve a large number of files and one must search for
    essentially  all x-matches.  If there are 100 files
    one needs to look at potentially 100 x 99 choices, and 
    in some cases we have close to 1000 files

    Running it with as many threads as possile is helpful.

                                       
History:

230505 ksl Coding begun and routine parallelized
230617 ksl Adapted to new version of routines
230704 ksl The mode and various additional stats have
            been added. The mode is now used for
            the purpose of setting the background 
            levels

'''

import sys
import os
import psutil
from glob import glob
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
from astropy.stats import sigma_clipped_stats
from multiprocessing import Pool
from log import *
import shutil
from medianrange import *
import timeit
import gc
import time

ierror=False

# import warnings
# warnings.filterwarnings('error', category=SyntaxWarning)


def halfsamplemode(inputData, axis=None):

    """Compute mode of an array using the half-sample algorithm"""

    if axis is not None:
        return np.apply_along_axis(halfsamplemode, axis, inputData)

    vals = np.sort(inputData.ravel())
    n = len(vals)
    while n > 2:
        nhalf = (n+1)//2
        isub = (vals[n-nhalf:]-vals[:nhalf]).argmin()
        vals = vals[isub:isub+nhalf]
        n = nhalf
    if n == 2:
        return 0.5*(vals[0]+vals[1])
    else:
        return vals[0]

def quartile_stats(xdata):
    '''
    Compute skew, etc, of a distribution with outlier 
    rejection

    See: https://towardsdatascience.com/skewness-and-kurtosis-with-outliers-f43167532c69
    '''

    vals=np.sort(xdata.ravel())

    # eliminate outliers
    n=len(vals)
    q3=vals[int(0.75*n)]
    q1=vals[int(0.25*n)]
    iqr=q3-q1
    xmin=q1-1.5*iqr
    xmax=q3+1.5*iqr

    xvals=vals[vals>xmin]
    xxvals=xvals[xvals<xmax]


    std=np.std(xxvals)
    mean=np.average(xxvals)
    med=np.median(xxvals)

    skew=(xxvals-mean)**3/std**3
    skew=np.average(skew)

    kurt=(xxvals-mean)**4/std**4
    kurt=np.average(kurt)

    # print('%10.3f %10.3f %10.3f %10.3f %10.3f' % (med,mean,std,skew,kurt))

    return skew,kurt







BACKDIR='DECam_BACK'

def calculate_one(file,xmatch_files,indir,calc_sigma_clipped=False,npix_min=100):
    '''
    Get statistics and make an estimate in
    the difference backgroound from the overlap region of
    images all of which must have the same 
    WCS (or actually image size)  

    one file to a number of file

    where:
    file is the name of a fits files files for
        which serves as a base file
    xmatch_files is a list of fits 
        for which offsets will be calculated
    imdir is the directory where the fits
        files are locatd
    calc_sigma_clipped means that in addtition to
        calculationg the standard set of statistics, sigma
        cliped values of the mean,median, and std are
        calculated. This increases the program run
        time significantly and so is used here as
        an option.
    npix_min is the minimum numver of pixels
        that must overlap the image

    The routine returns a list of records which 
    gives various statistics of the overlap regions,
    including the mode, mean,median and std 
    in the overlap region
    '''
    # print('Starting %s  with %d matches' % (file,len(xmatch_files)),flush=True)
    global ierror

    one=fits.open('%s/%s' % (indir,file))
    records=[]
        
    j=0
    while j<len(xmatch_files):
        one_record=[]
        # print('Diag: starting %d %s' % (j,xmatch_files[j]),flush=True)
        two=fits.open('%s/%s' % (indir,xmatch_files[j]))
        try:
            overlap=np.select([one[0].data*two[0].data!=0],[1],default=0)
            nonzero=np.count_nonzero(overlap)
        except:
            log_message('Error: Problem comparing %s and %s' % (file,xmatch_files[j]))
            nonzero=0
            ierror=True

        if nonzero>npix_min:
            ok=True
            one_record=[file,xmatch_files[j],nonzero]

            # We now need a mask, but we now need the good values 0 and the masked
            # values to be something else, so we subtract 1
            xmask=overlap-1
            one_masked=np.ma.array(one[0].data,mask=xmask)
            two_masked=np.ma.array(two[0].data,mask=xmask)
            xxdelta=two_masked-one_masked
            xdelta=xxdelta.compressed()
            # med_true,med_low,med_hi=medianrange(xdelta)
            mean=np.ma.average(xdelta)
            med=np.ma.median(xdelta)
            std=np.ma.std(xdelta)
            if np.isnan(med):
                ok=False

            if ok:
                # add the unbiased values
                one_record=one_record+[mean,med,std]
                # add the clipped values
                if calc_sigma_clipped: 
                    mean1,median1,std1=sigma_clipped_stats(xxdelta,sigma_lower=2,sigma_upper=2,grow=3)
                    if np.isnan(median1):
                        ok=False
                else:
                    mean1=-99.
                    median1=-99.
                    std1=-99.

            # Add the skew and kurtosis
            if ok:
                one_record=one_record+[mean1,median1,std1]
                skew,kurt=quartile_stats(xdelta)
                if np.isnan(skew):
                    ok=False

            if ok:
                one_record=one_record+[skew,kurt]
                # add the mode
                xmode=halfsamplemode(xdelta)
                if np.isnan(xmode):
                    ok=False

            if ok:
                one_record=one_record+[xmode]
                # print('Diag: Succeeded for files %s and %s' % (file,xmatch_files[j]),flush=True)
                records.append(one_record)
            else:
                print('BackStats: Failed with (NaNs for) files %s and %s' % (file,xmatch_files[j]),flush=True)
        # else:
        #     print('BackStats: Failed with %d non_zero < %s pixels for files %s and %s' % (nonzero,npix_min,file,xmatch_files[j]))

        two.close()
        del two
        gc.collect()
        j+=1
    one.close()
    print('BackStats: Finished %3d x-matches for  %s/%s ' % (len(records),indir,file),flush=True)
    return records


        
        
def do_one_tile(field='LMC_c42',tile='T07',nproc=1,calc_sigma_clipped=False):

    sumfile='Summary/%s_%s_overlap.txt' % (field,tile)
    try:
        x=ascii.read(sumfile)
    except:
        print('BackStats: Could not find overlap file %s ' % sumfile)
        raise IOError


    # Read the summary file so we can attached the exposure time 

    imsumfile='Summary/%s_%s.txt' % (field,tile)
    try:
        imsum=ascii.read(imsumfile)
        imsum=join(imsum,x['Filename','XEXPTIME'],join_type='left')
    except:
        print('BackStats: Could not read imsum file: %s' % imsum)
        raise IOError


    indir='%s/%s/%s' % (BACKDIR,field,tile)
    if os.path.exists(indir)==False:
        print('BackStats: Could not find the directory %s with the data' % indir)
        raise IOError

    
    qfiles=np.unique(x['Filename'])

    records=[]

    i=0
    all_inputs=[]
    while i<len(qfiles):
        files=x[x['Filename']==qfiles[i]]
        all_inputs.append([qfiles[i],files['XFilename'],indir,calc_sigma_clipped])
        i+=1
    
    log_message('BackStats:  %s %s has %d files to process' % (field,tile,len(all_inputs)))

    if nproc<2:
        i=0
        while i<len(all_inputs):
            xrecords=calculate_one(all_inputs[i][0],all_inputs[i][1],all_inputs[i][2],all_inputs[i][3])
            records=records+xrecords
            print('Debug: calculate_one complete (%d) for %s len %d' % (i, all_inputs[i][0],len(xrecords)))
            i+=1
    else:
        with Pool(nproc) as p:
            zrecords=p.starmap(calculate_one,all_inputs)
        records=[]
        print('Debug: finished starmap',flush=True)
        for one in zrecords:
            records=records+one

    log_message('Backstats: calculate_one is complete for all files in %s %s' % (field,tile))
    xrecords=np.array(records)
    
    xtab=Table(xrecords,names=['file1','file2','npix','mean','med','std','mean_clipped','med_clipped','std_clipped','skew','kurt','Delta'])

    os.sync()

    # Next lines are a cheat to get the formats to be sensible
    print('Debug: xtab with %d lines is complete' % len(xtab),flush=True)
    xtab.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)
    time.sleep(1.0)
    xtab=ascii.read('foo.txt')
    print('Debug: read back foo.txt with %d lines' % (len(xtab)),flush=True)

    # Now rescale to account for pixel size
    scale_factor=(0.2631*0.2631)/(2.*2.)

    xtab['mean']*=scale_factor
    xtab['med']*=scale_factor
    xtab['std']*=scale_factor
    if calc_sigma_clipped:
        xtab['mean_clipped']*=scale_factor
        xtab['med_clipped']*=scale_factor
        xtab['std_clipped']*=scale_factor
    xtab['Delta']*=scale_factor

    print('Debug: rescaled' ,flush=True)

    # colnames=xtab.colnames
    # i=0
    # while i<len(colnames):
    #     z=xtab[colnames[i]]
    #     xxx=[]
    #     xfloat=True
    #     for one_value in z:
    #         try:
    #             xxx.append(eval(one_value))
    #         except:
    #             xfloat=False
    #     if xfloat==True:
    #         xtab.remove_column(colnames[i])
    #        xtab[colnames[i]]=xxx
    #    i+=1

    xtab.info()

    colnames = xtab.colnames.copy()

    for colname in colnames:
        col = xtab[colname]
    
        # Skip if already numeric
        if col.dtype.kind in 'biufc':
            continue
    
        print(f'Debug: Processing column {colname}...', flush=True)
    
        try:
            # Use astropy's built-in conversion - much safer than eval()
        
            # Try to convert to float
            numeric_col = col.astype(float)
        
            # Replace the column
            xtab.remove_column(colname)
            xtab[colname] = numeric_col
            print(f'Debug:   Converted {colname} to float', flush=True)
        
        except (ValueError, TypeError):
            # Keep as string if conversion fails
            print(f'Debug:   Kept {colname} as string', flush=True)
            continue

    print('Debug: column conversion complete', flush=True)        

    print('Debug: removed columns' ,flush=True)

    ztab=xtab[colnames]

    print('Debug: created ztab' ,flush=True)


    ztab['mean'].format='.3f'
    ztab['mean'].format='.3f'
    ztab['med'].format='.3f'
    ztab['std'].format='.3f'
    ztab['med_clipped'].format='.3f'
    ztab['mean_clipped'].format='.3f'
    ztab['std_clipped'].format='.3f'
    ztab['skew'].format='.3f'
    ztab['kurt'].format='.3f'
    ztab['Delta'].format='.3f'



    print('Debug: assembled ztab',flush=True)

    # Now attach the filter and exposure time to this
    # The code below is to try to deal with a memory leak problme in astropy join

    imsum.rename_column('Filename','file1')
    
    gc.collect()
    print('Debug: analyzing tables before join', flush=True)

    # Check table sizes and memory usage
    print(f'Debug: ztab: {len(ztab)} rows, {len(ztab.colnames)} columns')
    print(f'Debug: imsum subset: {len(imsum)} rows, 5 columns')

    #Check for indices (still the #1 crash cause)
    print(f'Debug: ztab indices: {len(ztab.indices) if hasattr(ztab, "indices") else 0}')
    print(f'Debug: imsum indices: {len(imsum.indices) if hasattr(imsum, "indices") else 0}')

    # Check column names and types
    print(f'Debug: ztab columns: {ztab.colnames}')
    imsum_subset = imsum['file1','Image','FILTER','EXPTIME','XEXPTIME']
    print(f'Debug: imsum subset columns: {imsum_subset.colnames}')

    # CRITICAL: Check for common columns (join keys)
    common_cols = set(ztab.colnames) & set(imsum_subset.colnames)
    print(f'Debug: Common columns for join: {common_cols}')

    if not common_cols:
        print('ERROR: No common columns found - join will fail!')
    else:
        # Check data types of join columns
        for col in common_cols:
            print(f'Debug:   {col}: ztab dtype={ztab[col].dtype}, imsum dtype={imsum_subset[col].dtype}')

    # Check memory usage
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024**2
    print(f'Debug: Current memory usage: {mem_mb:.1f} MB')

    gc.collect()

    ztab.info()
    imsum.info()
    print('Debug: before join', flush=True)

    # print('Debug: before join',flush=True)

    ztab_result=join(ztab,imsum['file1','Image','FILTER','EXPTIME','XEXPTIME'],join_type='left')
    del ztab
    ztab=ztab_result
    del ztab_result
    gc.collect()

    print('Debug: joined ztab',flush=True)

    out_name='Summary/%s_%s_xxx.txt' % (field,tile)
    ztab.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    print('Debug: Wrote %s with %d lines' % (out_name,len(ztab)))
    return ztab

def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    nproc=0
    xremove=False
    calc_sigma_clipped=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-rm':
            xremove=True
        elif argv[i]=='-sigma_clipped':
            calc_sigma_clipped=True
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i][0]=='-':
            print('BackStats: Error: Unknown switch %s' % argv[i])
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
        print('BackStats: The tiles to be processed must be listed after the field, unless -all is invoked')
    open_log('%s.log' % field,reinitialize=False)
    for one in tiles:
        if calc_sigma_clipped:
            log_message('BackStats: Starting %s %s with %d processors (w/ sigma_clipped stats)' % (field,one,nproc))
        else:
            log_message('BackStats: Starting %s %s with %d processors (w/o sigma_clipped_stats)' % (field,one,nproc))


        start_time=timeit.default_timer()
        try:
            do_one_tile(field,one,nproc,calc_sigma_clipped)
        except IOError:
            log_message('BackStats: Failed on %s %s' % (field,one))
            return

        current_time=timeit.default_timer()
        log_message('BackStats: Finished %s %s in %.1f s' % (field,one,current_time-start_time))

    if ierror:
        log_message('Error: This run of BaskStats.py had errors')
        log_message('Error: BackPrep files not removed')
    elif xremove:
        for one in tiles:
            xdir='%s/%s/%s' % (BACKDIR,field,one)
            print('BackStats: There were no errors reported, removing %s' % (xdir))
            shutil.rmtree(xdir)


    close_log()
    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


