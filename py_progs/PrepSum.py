#!/usr/bin/env python
# coding: utf-8

'''

Synopsis:  

This produces summary tables of all of the data processed with PrepFiles


Command line usage (if any):

    usage: PremSum.py [-all]  Field_name [one or more tiles]

where
      -all    --> do all tiles for a fields

If -all is not provide then there must be one or more tile numbes, e. g. 1 5 7, listed


Description:  

Prepare a Summary of the files that have been to be used in prepped  for the Creation of combined images
with swarp.


* This asumes that PrepFiles has been run for a tile or tiles.  
* The next step is combine files


Primary routines:


Notes:

                                       
History:
    230504 - Modified so files gets fits files with and without .fz extension
    230516 - Replaced SumFile.py with this routine that follows new directory structure


'''

import sys
from astropy.io import ascii
import timeit
import os
from astropy.io import fits
from astropy.table import Table
from glob import glob
import numpy as np

PREPDIR='DECam_PREP'



def get_tile_files(field='LMC_c42',tile='T07'):
    '''
    Verify that the directory exists where the individual observations
    files are stored for a tile, and return the names of the tile files
    '''
    
    tile_dir='%s/%s/%s/' % (PREPDIR,field,tile)
    print('The tile directory is ',tile_dir)
    if os.path.isdir(tile_dir)==False:
        print('Error: the tile directory %s does not exist' % tile_dir)
        return tile_dir,[]
    
    xfiles='%s/c*fits*' % tile_dir
    files=glob(xfiles)

    print('get_tile_files found %d files' % len(files))

    if len(files)==0:
        print('Error no files found for %s' % xfiles)
    
    i=0
    while i<len(files):
        words=files[i].split('/')
        files[i]=words[-1]
        i+=1
    
        
    
    return tile_dir,files



def get_one_keyword(keyword='object',ext=0,filenames='',work_dir=''):
    
    value=[]
    nerr=0
    for file in filenames:
        x=fits.open('%s/%s' % (work_dir,file))
        try:
            z=x[ext].header[keyword]
            value.append(z)
        except:
            nerr+=1
            value.append('unknown')
    
    if nerr>0:
        print('Failed to find %s in extension %d in %d of %d attempts' % (keyword,ext,nerr,len(filenames)))
    
    return value



def make_summary(filenames='',work_dir='',outroot=''):

    start_time = timeit.default_timer()
    objects=get_one_keyword('object',0,filenames,work_dir)
    xtab=Table([filenames,objects],names=['Filename','Object'])

    filter=get_one_keyword('filter',0, filenames,work_dir)
    xfilter=[]
    for one in filter:
        zfilt='Unknown'
        if one.count('N708'):
            zfilt='N708'
        elif one.count('N673'):
            zfilt='SII'
        elif one.count('N662'):
            zfilt='Ha'
        elif one[0]=='r':
            zfilt='R'

        xfilter.append(zfilt)
    xtab['Filter']=xfilter
    
    xtab['Detector']=get_one_keyword('DETPOS',1,filenames,work_dir)
    
    xtab['Exptime']=get_one_keyword('exptime',0,filenames,work_dir)
    
    gaina=np.array(get_one_keyword('gaina',1,filenames,work_dir))
    gainb=np.array(get_one_keyword('gainb',1,filenames,work_dir))
    xtab['Gain']=0.5*(gaina+gainb)
    xtab['Gain'].format='.2f'
    xtab['MAGZERO']=get_one_keyword('MAGZERO',0,filenames,work_dir)
    xtab['MAGZPT']=get_one_keyword('MAGZPT',0,filenames,work_dir)    
    xtab['SEEING']=get_one_keyword('SEEING',0,filenames,work_dir)   
    xtab['Sky']=get_one_keyword('AVSKY',1,filenames,work_dir)
    xtab['RA']=get_one_keyword('CENRA1',1,filenames,work_dir)
    xtab['Dec']=get_one_keyword('CENDEC1',1,filenames,work_dir)
    xtab.sort(['Filter','Exptime','Object','Detector','Filename'])
    
    if outroot=='':
        foo=work_dir.replace('/',' ')
        foo=foo.strip()
        words=foo.split()
        print(work_dir)
        print(words)
        outroot='%s_%s_imsum' % (words[-2],words[-1])
        print ('Outroot ',outroot)
        
    xtab.write('Summary/%s.txt' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    elapsed = timeit.default_timer() - start_time
    print('All done in %f s' % elapsed)   
    return xtab



def do_one_tile(field='LMC_c42',tile='T07'):
    tile_dir,files=get_tile_files(field,tile)
    if len(files)==0:
        print('Houston: There is a problem for field %s and tile %s.  No files were found' % (field,tile))
        return
    make_summary(filenames=files,work_dir=tile_dir,outroot='')

    summarize(field,tile)
    return


def do_all_tiles(field='LMC_c42'):
    '''
    Prepare all of the tiles in one set of data
    '''

    i=1
    while i<17:
        do_one_tile(field,tile='T%02d' % i)
        i+=1

def get_sum_tab(field='LMC_c42',tile='T07'):
    '''
    Read a table that summarizes information about
    the various files athat are used for a tile

    The location of the table files is currently hardocaded
    '''

    xname='Summary/%s_%s_imsum.txt' % (field,tile)
    xtab=ascii.read(xname)
    return xtab

def summarize(field='LMC_c42',tile='T07'):
    '''
    Summarize the various fits files that are relevant to
    a specific tile

    History

    230504 - Moved this functionality from Swarp.py to SumFiles.py, and 
    created an summary in the Summary directory)
    '''

    f=open('Summary/%s_%s.exsum.txt' % (field,tile),'w')

    string=('Summary for field %s tile %s' % (field,tile))

    print(string)
    f.write('%s/n' % string)

    x=get_sum_tab(field,tile)

    ra=np.average(x['RA'])
    dec=np.average(x['Dec'])

    string='The center of this tile is %.5f  %.5f' % (ra,dec)

    print(string)
    f.write('%s\n' % string)


    ha=x[x['Filter']=='Ha']
    s2=x[x['Filter']=='SII']
    r=x[x['Filter']=='R']
    n708=x[x['Filter']=='N708']

    string=('Ha  images  : %3d' % len(ha))
    print(string)
    f.write('%s\n' % string)

    string=('SII images  : %3d' % len(s2))
    print(string)
    f.write('%s\n' % string)

    string=('R   images  : %3d' % len(r))
    print(string)
    f.write('%s\n' % string)

    string=('N708 images : %3d' % len(n708))
    print(string)
    f.write('%s\n' % string)

    string=('\nHa')
    print(string)
    f.write('%s\n' % string)

    times,counts=np.unique(ha['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        string=('   exp %8s  no %4d' % (times[i],counts[i]))
        print(string)
        f.write('%s\n' % string)
        i+=1

    string=('\nSII')
    print(string)
    f.write('%s\n' % string)
    times,counts=np.unique(s2['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        string=('   exp %8s  no %4d' % (times[i],counts[i]))
        print(string)
        f.write('%s\n' % string)
        i+=1

    string=('\nR')
    print(string)
    f.write('%s\n' % string)

    times,counts=np.unique(r['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        string=('   exp %8s  no %4d' % (times[i],counts[i]))
        print(string)
        f.write('%s\n' % string)
        i+=1

    string=('\nN708')
    print(string)
    f.write('%s\n' % string)

    times,counts=np.unique(n708['Exptime'],return_counts=True)
    i=0
    while i<len(times):
        string=('   exp %8s  no %4d' % (times[i],counts[i]))
        print(string)
        f.write('%s\n' % string)
        i+=1

    print('\n')
    f.close()
    return ra,dec


    
def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
        elif field=='':
            field=argv[i]
        else:
            tiles.append(argv[i])
        i+=1

    if xall:
        do_all_tiles(field)
    else:
        if len(tiles)==0:
            print('The tiles to be processed must be listed after the field, unless -all is invoked')
        for one in tiles:
            do_one_tile(field,one)



    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)

