#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create a table (or tables) that contains information about 
how overlap regions of various images compare


Command line usage (if any):

    usage: BackStats.py  [-np 8] [-all] [-rm] field T01 T02 ...

    where field is a field name, and the following numbers are one 
    or more tiles

    -all implies to do all tiles
    -rm causes the files that were created by BackPrep to be deleted
    if the stats were carried out successfully.  This only occurs
    at the end of the program so that if something fails along
    the way, none of the backround files will have been deleted

    -np 8  implies to perform calculations in parallel, with in
    this case 8 threads

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

'''

import sys
import os
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

def calculate_one(file,xmatch_files,indir,npix_min=100):
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
    npix_min is the minimum numver of pixels
        that must overlap the image

    The routimee returns a list of records which 
    give the mean,median and std of each file 
    in the overlap region
    '''
    # print('Starting %s  with %d matches' % (file,len(xmatch_files)))

    one=fits.open('%s/%s' % (indir,file))
    records=[]
        
    j=0
    while j<len(xmatch_files):
        one_record=[]
        two=fits.open('%s/%s' % (indir,xmatch_files[j]))
        overlap=np.select([one[0].data*two[0].data!=0],[1],default=0)
        nonzero=np.count_nonzero(overlap)

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
            med_true,med_low,med_hi=medianrange(xdelta)
            mean=np.ma.average(xdelta)
            med=np.ma.median(xdelta)
            std=np.ma.std(xdelta)
            if np.isnan(med):
                ok=False
                continue

            # add the unbiased values
            one_record=one_record+[mean,med,std]

            # add the clipped values
            mean1,median1,std1=sigma_clipped_stats(xxdelta,sigma_lower=2,sigma_upper=2,grow=3)
            if np.isnan(median1):
                ok=False
                continue

            # Add the skew and kurtosis
            one_record=one_record+[mean1,median1,std1]
            skew,kurt=quartile_stats(xdelta)
            if np.isnan(skew):
                ok=False
                continue
            one_record=one_record+[skew,kurt]

            # add the mode
            xmode=halfsamplemode(xdelta)
            one_record=one_record+[med]

            if ok==True:
                records.append(one_record)

        two.close()
        j+=1
    one.close()
    print('Finished %3d x-matches for  %s/%s ' % (len(records),indir,file))
    return records


        
        
def do_one_tile(field='LMC_c42',tile='T07',nproc=1):

    sumfile='Summary/%s_%s_overlap.txt' % (field,tile)
    try:
        x=ascii.read(sumfile)
    except:
        print('Could not find overlap file %s ' % sumfile)
        raise IOError


    # Read the summary file so we can attached the exposure time 

    imsumfile='Summary/%s_%s.txt' % (field,tile)
    try:
        imsum=ascii.read(imsumfile)
    except:
        print('Could not read imsum file: %s' % imsum)
        raise IOError


    indir='%s/%s/%s' % (BACKDIR,field,tile)
    if os.path.exists(indir)==False:
        print('Could not find the directory %s with the data' % indir)
        raise IOError

    
    qfiles=np.unique(x['Filename'])

    records=[]

    i=0
    all_inputs=[]
    while i<len(qfiles):
        files=x[x['Filename']==qfiles[i]]
        all_inputs.append([qfiles[i],files['XFilename'],indir])
        i+=1
    
    log_message('BackStats:  %s %s has %d files to process' % (field,tile,len(all_inputs)))

    if nproc<2:
        i=0
        while i<len(all_inputs):
            xrecords=calculate_one(all_inputs[i][0],all_inputs[i][1],all_inputs[i][2])
            records=records+xrecords
            i+=1
    else:
        with Pool(nproc) as p:
            zrecords=p.starmap(calculate_one,all_inputs)
        records=[]
        for one in zrecords:
            records=records+one

    xrecords=np.array(records)
    
    xtab=Table(xrecords,names=['file1','file2','npix','mean1','med1','std1','mean2','med2','std2','skew','kurt','Delta'])
    # Next lines are a cheat to get the formats to be sensible
    xtab.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)
    xtab=ascii.read('foo.txt')

    # Now rescale to account for pixel size
    scale_factor=(0.27*0.27)/(2.*2.)

    xtab['mean1']*=scale_factor
    xtab['med1']*=scale_factor
    xtab['std1']*=scale_factor
    xtab['mean2']*=scale_factor
    xtab['med2']*=scale_factor
    xtab['std2']*=scale_factor
    xtab['Delta']*=scale_factor

    colnames=xtab.colnames
    i=0
    while i<len(colnames):
        z=xtab[colnames[i]]
        xxx=[]
        xfloat=True
        for one_value in z:
            try:
                xxx.append(eval(one_value))
            except:
                xfloat=False
        if xfloat==True:
            xtab.remove_column(colnames[i])
            xtab[colnames[i]]=xxx
        i+=1
        

    ztab=xtab[colnames]
    ztab['mean1'].format='.3f'
    ztab['mean2'].format='.3f'
    ztab['med1'].format='.3f'
    ztab['med2'].format='.3f'
    ztab['std1'].format='.3f'
    ztab['std2'].format='.3f'
    ztab['skew'].format='.3f'
    ztab['kurt'].format='.3f'
    ztab['Delta'].format='.3f'




    # Now attach the filter and exposure time to this

    imsum.rename_column('Filename','file1')

    ztab=join(ztab,imsum['file1','FILTER','EXPTIME'],join_type='left')

    out_name='Summary/%s_%s_xxx.txt' % (field,tile)
    ztab.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d lines' % (out_name,len(ztab)))
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

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-rm':
            xremove=True
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
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
    open_log('%s.log' % field,reinitialize=False)
    for one in tiles:
        log_message('BackStats: Starting %s %s with %d processors' % (field,one,nproc))
        start_time=timeit.default_timer()
        try:
            do_one_tile(field,one,nproc)
        except IOError:
            log_message('BackStats: Failed on %s %s' % (field,one))
            return

        current_time=timeit.default_timer()
        log_message('BackStats: Finished %s %s in %.1f s' % (field,one,current_time-start_time))
    close_log()

    if xremove:
        for one in tiles:
            xdir='%s/%s/%s' % (BACKDIR,field,one)
            print('Removing %s' % (xdir))
            shutil.rmtree(xdir)


    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


