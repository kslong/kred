#!/usr/bin/env python
# coding: utf-8

'''

Create Combined images of the different filters in the field using swarp

This can only be run after PrepFiles, SumFiles and SwarpSetup have benn run.    
SwarpSetup will have prepared the run directories and the inpus.

To run swarp from the command line (one must be in the normal run directory) since
a particular directory structure is assumed here)

Usage:   Swarp.py [-all] [-bsub] field [tiles]

where -all will cause swarp to be run on all 16 tiles.

and   -bsub will cause swarp to be run in the DECAM_SWARP2 directories, and should
be operating on files that are in background subtracted DECAM_PREP2 directories

If one wants to run only 1 or a few tiles then the command will be something like

Swarp.py LMC_c42  T01 T03 T05

Notes:
    It would be sensible to combine SetupSwarp and Swarp
'''


import os, stat
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from glob import glob
import timeit
import subprocess
from log import *
from kred_repo_info import get_kred_git_commit_info
from astropy.wcs import WCS

exec_dir=os.getcwd()


SWARPDIR='DECam_SWARP'
PREPDIR=os.path.abspath('DECam_PREP')





def add_header_keywords(root,run_dir='./',field='',tile=''):
    '''
    Add header keywords to an the results of Swarp
    '''
    if run_dir=='':
        run_dir='./'
    swarp_file='%s/%s.fits' % (run_dir,root)
    if os.path.isfile(swarp_file)==False:
        print('Error: swarp: Could not find %s to add header info' % swarp_file)
        return
    inputs='%s/%s.in' % (run_dir,root)
    if os.path.isfile(inputs)==False:
        print('Error: swarp: Could not find file inputs %s to add header info' % inputs)
        return

    f=open(inputs)
    lines=f.readlines()
    hdr_file=lines[0].strip()
    try:
        xhdr=fits.open(hdr_file)
        hdr=xhdr['PRIMARY'].header
    except:
        print('Failed to open %s' % hdr_file)
        return

    try:
        swarp=fits.open(swarp_file,mode='update')
        shdr=swarp['PRIMARY'].header
    except:
        print('Failed to open %s' % swarp_file)
        return

    xobject=field
    if len(tile)>0:
        xobject='%s_%s' % (field,tile)
    shdr['PROPID']='DeMCELS'
    shdr['OBJECT']=xobject
    shdr['EXPTIME']=hdr['EXPTIME']
    shdr['TELESCOP']=hdr['TELESCOP']
    shdr['OBSERVAT']=hdr['OBSERVAT']
    shdr['INSTRUME']=hdr['INSTRUME']
    shdr['FILTER']=hdr['FILTER']
    shdr['N_EXP']=len(lines)
    shdr['PROCTYPE']='Stacked'

    commit,commit_date=get_kred_git_commit_info()

    shdr['PIPELINE'] = 'kred'
    shdr['COMMIT'] = commit
    shdr['PDATE']=commit_date

    naxis1=shdr['NAXIS1']
    naxis2=shdr['NAXIS2']

    corners_pixels = np.array([[0, 0], [0, naxis2-1], [naxis1-1, 0], [naxis1-1, naxis2-1]])
    wcs=WCS(shdr)

    corners_ra_dec = wcs.all_pix2world(corners_pixels, 0)

    print(corners_ra_dec)
    for i in range(4):
        shdr[f'COR{i+1}RA1'] = corners_ra_dec[i, 0]
        shdr[f'COR{i+1}DEC1'] = corners_ra_dec[i, 1]

    swarp.flush()
    swarp.close()


    return



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
        run_dir='%s2/%s/%s/' % (SWARPDIR,field,tile)


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

        # Note that when this routine is run we are in the actual directory where Swarp is tun
        xroot=one.replace('.run','')
        add_header_keywords(xroot,'',field,tile)

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
        if bsub :
            log_message('Swarp: Start run on %s %s with modified background' % (field,one))
            one='%s' % one
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


