#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Sumarize mef files that constitute the raw data from
DECam observations of the LMC and SMC

This routines assume a standard setup with the a
local link to the raw DECam mef files (or part of
them with the same format as on Box


Command line usage (if any):

    usage: MefSum.py [-h] [-all] [-np 3]  [-det] [-mef] field1 field2 ...

    where 
        -h prints the documentation
        -np 3  causes the processing to be carried out with a a given no of threads
        -all causes all files in the MEF directories to be processed. This should
        only be sued with cauthion since it will take considerable time, and so
        the user is asked to confirm this option.
        -det just runs the individual ccd portion of the process (this is diagnostic)
        -mef just runs the overall mef portions (this is diagnostic)
        -sigma_clipped (Causes some statitics to be calculaed using sigma clipping).
            This is a diagnostic mode and is not normally used.  The estimates
            are now based on the mode

    and fields comprise one or more fields e.g LMC_c42  to be processed

Description:  

    The routine reads the MEF files and summarizes information from
    the headers in two tables, a mef.tab file and a det.tab file
    The routine also calculates some statstics associated with the
    counts in each detector, that are used to carry out an
    intial background subtraction.

Primary routines:


Notes:

    At present this does not check that a field has been processed

    Here multiple processors are used only to handle multiple fields,
    and so if one as for 4 threads, but only provides 3 fields
    only 3 threads will be used.  (This is something that should
    be changed, but requires a slightly different set of 
    multiprocessing commands)


                                       
History:

230513 ksl Coding begun

'''

import sys
import os

from astropy.io import fits, ascii
from glob import glob
from astropy.stats import sigma_clipped_stats
from astropy.table import Table,vstack
from astropy.wcs import WCS
import numpy as np
import timeit
import time
import multiprocessing
from multiprocessing import Pool
from log import *
from medianrange import *
multiprocessing.set_start_method("spawn",force=True)




# # Summarize a set of MEF files used for DECam data
# 
# The purpose of this routine is to created a summary of the files that appear in a particular collection (usually
# and exposure sequence of MEF files).  This is an overall summary of the MEF files that are available
# and does not itself have information about the individual CCD detectors (which is dealt with else where)

# MEFDIR is the location of the the root directory for all of the directories that conatina subsets of MEF files (e.g. LMC_c42) that are going to be used. 
# One shouuld be wary of writing extraneous files to the MEDDIR
# 
# 
# Note that the MEF dir can actually just contain symbolic links to the real locations for the files, which might be spread over several disks, but
# this will be the "master" directory that kred uses to find data.
# 
# These routines create a directory called Summary in the MEDDIR that will be used for storing various
# tables that describe the data.



MEFDIR='DECam_MEF/'

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



def get_keyword(key,ext):
    '''
    where key is just the name of a keyword
    and ext points to a specific extension
    of an already opened fits file, e.g x[0]
    
    Note that one does not include .header in
    the cal

    History:

    230520 - Fix to handle the fact that there
    are some files for which the Band is not filled
    in.
    
    '''
    if key=='FILTER':
        try:
            value=ext.header['FILTER']
            if value.count('N662'):
                answer='N662'
            elif value.count('N673'):
                answer='N673'
            elif value.count('N708'):
                answer='N708'
            elif value.count('r DECam'):
                answer='r'
            else:
                print('Could not identify filter from ',answer)
                answer='Unknown'
            return answer
        except:
            answer='Unknown'
        return answer

    try:
        answer=ext.header[key]
    except:
        answer='Unknown'
    return answer
                          

def get_mef_overview(field='LMC_c45'):
    

    files=glob('%s/%s/mef/*ooi*fz' % (MEFDIR,field))
    if len(files)==0:
        print('Error: could not locate any files with %s/%s/mef/*ooi*fz' % (MEFDIR,field))
        return []

    records=[]
    
    keys=['OBJECT','FILTER','EXPTIME','MAGZERO','SEEING','OBSID']
    
    for one_file in files:
        x=fits.open(one_file)
        record=[]
        name=one_file.split('/')
        name=name[-1]
        name=name.replace('.fits.fz','')  
        record=[name]  
        for one_key in keys:
            record.append(get_keyword(one_key,x[0]))
        
        records.append(record)
        x.close()
        
    records=np.array(records)   
    
    
        
    tab_names=['Root']+keys
    
    xtab=Table(records,names=tab_names)
    xtab['Field']=field
    
    # Next step is a cheat to improve translate numbers to numbers rather than strings
    
    
    xtab.write('goo_%s.txt' % field ,format='ascii.fixed_width_two_line',overwrite=True)
    ztab=ascii.read('goo_%s.txt' %field)
    os.remove('goo_%s.txt' %field)
    
    ztab.sort(['FILTER','EXPTIME'])


    ztab.write('Summary/%s_mef.tab' % field ,format='ascii.fixed_width_two_line',overwrite=True)
    
    
    return 
    
                          

def get_det_overview(field='LMC_c45',calc_sigma_clipped=False):
    '''
    Get information about the individual CCDs, including 
    some statistics associated with the data that is 
    useful for estimating background
    '''
    
    files=glob('%s/%s/mef/*ooi*fz' % (MEFDIR,field))
    if len(files)==0:
        print('Error: could not locate any files with %s/%s/mef/*ooi*fz' % (MEFDIR,field))
        return []


    
    records=[]
    
    keys=['EXTNAME','CENRA1','CENDEC1','COR1RA1','COR1DEC1','COR2RA1','COR2DEC1','COR3RA1','COR3DEC1','COR4RA1','COR4DEC1']
    
    for one_file in files:
        start_time = timeit.default_timer()
        log_message('Starting %s ' % one_file)
        x=fits.open(one_file)

        name=one_file.split('/')
        name=name[-1]
        name=name.replace('.fits.fz','')  
 
        i=1
        while i<len(x):
            record=[name] 

            for one_key in keys:
                try:
                    record.append(get_keyword(one_key,x[i]))
                except:
                    log_message('Error: Problem with keyword %s for %s extension %d'% (one_key,one_file,i))
                    record.append(-99.9999)
            try:
                mode=halfsamplemode(x[i].data)
                if calc_sigma_clipped:
                    # print('sigma clipping',i)
                    mean,median,std=sigma_clipped_stats(x[i].data,sigma_lower=3,sigma_upper=2,grow=3)
                else:
                    # print('normal sigma ',i)
                    mean=np.average(x[i].data)
                    median=np.median(x[i].data)
                    std=np.std(x[i].data)
                    # mean=median=std=-99.
            except:
                log_message('Error: unable to get mean etc for %s extension %d'% (one_file,i))
                mode=mean=median=std=-999.

            record.append(mean)
            record.append(median)
            record.append(std)
            record.append(mode)
        
            records.append(record)
            i+=1
        log_message('finished %s in %.1f s' % (one_file,timeit.default_timer()-start_time))
        x.close()
        
    records=np.array(records)   
    
    
        
    tab_names=['Root']+keys+['Mean','Med','STD','Mode']
    
    xtab=Table(records,names=tab_names)
    xtab['Field']=field
    
    # Next step is a cheat to improve translate numbers to numbers rather than strings
    
    xtab.write('foo_%s.txt' % field ,format='ascii.fixed_width_two_line',overwrite=True)
    ztab=ascii.read('foo_%s.txt' %field)
    os.remove('foo_%s.txt' %field)


    bad=[]
    i=0
    while i<len(ztab):
        if ztab['COR1DEC1'][i]=='Unknown' or ztab['Mode'][i]==-999.:
            bad.append(i)
        i+=1


    if len(bad)>0:
        xbad=ztab[bad]
        print('Warning: There were %d CCD images of %d total without RAs and DECs in header' % (len(bad),len(ztab)))
        print('Writing the list of bad images to Summary/%s_det_bad.tab' % field )
        xbad.write('Summary/%s_det_bad.tab' % field ,format='ascii.fixed_width_two_line',overwrite=True)

        ztab.remove_rows(bad)
    
    
    print('Writing Summary/%s_det.tab with %d rows' % (field,len(ztab))) 
    ztab.write('Summary/%s_det.tab' % field ,format='ascii.fixed_width_two_line',overwrite=True)

    return 




def do_one(field='LMC_c45',do_mef=True,do_det=True,calc_sigma_clipped=False):
        '''
        Simply get all the information one wants about the imaages associated
        with a single field
        '''
        open_log('%s.log' % field,reinitialize=True)
        log_message('Starting MefSum on %s'% field)
        if do_mef:
            get_mef_overview(field)
        if do_det:
            get_det_overview(field,calc_sigma_clipped)
        close_log()
        return


def steer(argv):
    '''
    This is just a sterring routine, 
    and at present is not very useful,
    '''

    if os.path.isdir(MEFDIR)==False:
        print('Error: These routines assume a directory DECam_MEF exists, and contains the raw mef files')
        return

    if os.path.isdir('./Summary')==False:
        os.mkdir('./Summary')
        print('Created a new Summary directory')
    

    nproc=1
    fields=[]

    xall=False
    xdet=True
    xmef=True
    calc_sigma_clipped=False


    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-det':
            xdet=True
            xmef=False
        elif argv[i]=='-mef':
            xmef=True
            xdet=False
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i]=='-sigma_clipped':
            calc_sigma_clipped=True
        elif argv[i][0]=='-':
            print('Error: Unknown switch  %s' % argv[i])
            return
        else:
            fields.append(argv[i])
        i+=1

    if xall==True:
        xfields=glob('%s/*/mef/' % MEFDIR)

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
            fields.append(words[-3])
        fields=np.sort(fields)

    start_time = timeit.default_timer()

    if len(fields)<nproc:
        nproc=len(fields)
        print('Reducing the number of proocess threads to %d, which is the number of fields' % nproc)

    print('These are the fields that will be processed')
    for one in fields:
        print(one)

    if nproc==1:
        for one in fields:
            do_one(one,xmef,xdet,calc_sigma_clipped)
    else:
        xinputs=[]
        for one in fields:
            xinputs.append([one,xmef,xdet,calc_sigma_clipped])
        with Pool(nproc) as p:
            zrecords=p.starmap(do_one,xinputs)

    elapsed = timeit.default_timer() - start_time

    print('Processed %d fields with %d threads in %.1f s or %.1f s per field' % (len(fields), nproc,elapsed,elapsed/len(fields)))



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys

    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)


