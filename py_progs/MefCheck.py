#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Check whether one or more fits files have headers
that conform to the fits standard, by looking at
the warnings that astropy generaates


Command line usage (if any):

    usage: check.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

230706 ksl Coding begun

'''

import sys
import warnings
import re
import os

from astropy.io import fits, ascii
from glob import glob
from astropy.table import Table,vstack
import numpy as np
import timeit
import time
import multiprocessing
from multiprocessing import Pool
from log import *
multiprocessing.set_start_method("spawn",force=True)


MEFDIR='DECam_MEF/'



# Function to raise a test warning
def test(filename='c4d_190108_063145_ooi_N662_v1.fits.fz',outputfile='out.txt'):
    x=fits.open(filename)
    try:
        #original_stdout=sys.stdout
        #sys.stdout = open('/dev/null', 'w')
        i=0
        while i<len(x):
            if i==0:
                foo=x[i].header['EXPTIME']
            else:
                foo=x[i].header['EXTNAME']
            i+=1
        # x.info()
        # print(x[0].header['EXPTIME'])
        #sys.stdout=original_stdout
    except AttributeError:
        print('Got Attibutue Error')


# Function to catch and count warnings
def count_warnings(filename):
    warning_counts = {}
    original_showwarning = warnings.showwarning

    def catch_warning(message, category, filename, lineno, file=None, line=None):
        warning_key = re.split(r':\s*', str(message), 1)[0]
        warning_counts[warning_key] = warning_counts.get(warning_key, 0) + 1

    warnings.showwarning = catch_warning

    # Call the function that raises the warnings
    test(filename)

    warnings.showwarning = original_showwarning

    return warning_counts



def do_one(filename):

    warning_counts = count_warnings(filename)

    # Print the warning counts
    for warning_message, count in warning_counts.items():
        print(f"Warning: {warning_message}, Count: {count}")
    if len(warning_counts.items())==0:
        print('file %s is ok' % filename)
        return 'Good'
    else:
        print('file %s is NOK' % filename)
        return 'Bad'




def check_mefs(field='LMC_c45'):


    files=glob('%s/%s/mef/*ooi*fz' % (MEFDIR,field))
    if len(files)==0:
        print('Error: could not locate any files with %s/%s/mef/*ooi*fz' % (MEFDIR,field))
        return []

    xstat=[]
    for one in files:
        xx=do_one(one)
        xstat.append(xx)

    xfiles=[]
    for one in files:
        words=one.split('/')
        xfiles.append(words[-1])


    qtab=Table([xfiles,xstat],names=['Filename','Header_Quality'])
    qtab.write('Summary/%s_mef_qual.tab' % field,format='ascii.fixed_width_two_line',overwrite=True)



    return




def do_one_field(field='LMC_c45'):
        '''
        Simply get all the information one wants about the imaages associated
        with a single field
        '''
        # open_log('%s.log' % field,reinitialize=True)
        # log_message('Starting MefSum on %s'% field)
        check_mefs(field)
        # close_log()
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


    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch  %s' % argv[i])
            return
        else:
            fields.append(argv[i])
        i+=1

    if os.path.isdir('./Summary')==False:
        os.mkdir('./Summary')
        print('Created a new Summary directory')

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

    start_time = timeit.default_timer()

    if len(fields)<nproc:
        nproc=len(fields)
        print('Reducing the number of proocess threads to %d, which is the number of fields' % nproc)

    print('These are the fields that will be processed')
    for one in fields:
        print(one)

    if nproc==1:
        for one in fields:
            do_one_field(one)
    else:
        xinputs=[]
        for one in fields:
            xinputs.append([one])
        with Pool(nproc) as p:
            zrecords=p.starmap(do_one_field,xinputs)

    elapsed = timeit.default_timer() - start_time

    print('Processed %d fields with %d threads in %.1f s or %.1f s per field' % (len(fields), nproc,elapsed,elapsed/len(fields)))



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys

    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
