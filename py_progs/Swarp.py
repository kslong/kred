#!/usr/bin/env python
# coding: utf-8

'''

Create Combined images of the different filters in the field using swarp

This can only be run after PrepFiles and SumFiles have benn run.    The routines here generate 
inputs to run swarp which combines the individual CCD images into tile images.   This is normally
carried out inside a Jupyter scripts.


To run swarp from the command line (one must be in the normal run directory) since
a particular directory structure is assumed here)

Usage:   Swarp.py [-all] [-bsub] field [tiles]

where -all will cause swarp to be run on all 16 tiles.

and   -bsub will cause swarp to be run on the files in background subtracted _b directories

If one wants to run only 1 or a few tiles then the command will be something like

Swarp.py LMC_c42  T01 T03 T05


'''


import os, stat
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from glob import glob
import timeit
import subprocess
from log import *

exec_dir=os.getcwd()


SWARPDIR='DECam_SWARP'
PREPDIR=os.path.abspath('DECam_PREP')


def run_swarp(field='LMC_c42',tile='T07',bsub=False):
    '''
    run_all of the swarp commands in a particular field and tile

    The normal outputs from swarp are writeen to a .txt file,
    a few are written to the screen here so that one can see 
    that the routines have run correctly.  
    '''
    xstart=qstart=start_time=timeit.default_timer()

    if bsub==False:
        run_dir='%s1/%s/%s/' % (SWARPDIR,field,tile)
    else:
        run_dir='%s1/%s/%s/' % (SWARPDIR,field,tile)


    try:
        os.chdir(run_dir)
    except:
        print('Error: Could not cd to %s' % run_dir)
        os.chdir(exec_dir)

    run_files=glob('*.run')
    # print(run_files)
    nfiles=len(run_files)
    n=1
    for one in run_files:
        print('\n***Beginning %s (%d of %d)' % (one,n,nfiles))
        outfile=one.replace('.run','.txt')
        command=one
        xout=subprocess.run(command,shell=True,capture_output=True)
        current_time=timeit.default_timer()

        print('***Writing last portion of swarp outputs\n')

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

        print('***Finished writing end of outputs')

        print('***Finished %s (%d of %d) in %.1f s\n' % (one,n,nfiles,current_time-start_time))
        start_time=current_time


        n+=1


    xcurrent_time=timeit.default_timer()
    print('\n***Completely done for tile %s in %.1f s\n' % (tile,xcurrent_time-xstart))


    os.chdir(exec_dir)
    return 



def steer(argv):
    '''
    This is just a steering routine for running swarp on one or more
    tiles from the command line

    '''
    field=''
    tiles=[]
    xall=False
    redo=True
    bsub=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-bsub':
            bsub=True
        elif argv[i]=='-finish':
            redo=False
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

    open_log('%s.log' % field,reinitialize=False)

    for one in tiles:
        if bsub and one.count('_b')==0:
            log_message('Swarp: Start run on %s %s with modified background' % (field,one))
            one='%s_b' % one
        else:
            log_message('Swarp: Start run on %s %s ' % (field,one))

        run_swarp(field,one,bsub)

        log_message('Swarp: Finished run n %s %s ' % (field,one))

    close_log()


    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


