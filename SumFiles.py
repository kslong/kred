#!/usr/bin/env python
# coding: utf-8

'''

Synopsis:  

This produces summary tables of all of the data processed with PrepFiles


Command line usage (if any):

    usage: SumFiles.py [-all]  Field_name [one or more tiles]

where
      -all    --> do all tiles for a fields

If -all is not provide then there must be one or more tile numbes, e. g. 1 5 7, listed


Description:  

Prepare a Summary of the files that have been to be used in prepped  for the Creation of combined images
with swarp.


* This asumes that PrepFiles has been run for a tile or tiles.  
* This is the second step after the data have been downloaded, and the zeroth stage of photpipe has been run.  
* The next step is combine files


Primary routines:


Notes:
                                       
History:


'''

import sys
from astropy.io import ascii
import timeit
import os
from astropy.io import fits
from astropy.table import Table
from glob import glob
import numpy as np



def get_tile_files(field='LMC_c42',tile=7):
    '''
    Verify that the directory exists where the individual observations
    files are stored for a tile, and return the names of the tile files
    '''
    
    tile_dir='../workspace/%s/%d/' % (field,tile)
    print('The tile directory is ',tile_dir)
    if os.path.isdir(tile_dir)==False:
        print('Error: the tile directory %s does not exist' % tile_dir)
        return []
    
    xfiles='%s/L*fits.fz' % tile_dir

    files=glob(xfiles)
    print('get_tile_files found %d files' % len(files))
    
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
        
    xtab.write('../workspace/Summary/%s.txt' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    elapsed = timeit.default_timer() - start_time
    print('All done in %f s' % elapsed)   
    return xtab



def do_one_tile(field='LMC_c45',tile=7):
    tile_dir,files=get_tile_files(field,tile)
    make_summary(filenames=files,work_dir=tile_dir,outroot='')


def do_all_tiles(field='LMC_c42'):
    '''
    Prepare all of the tiles in one set of data
    '''

    i=1
    while i<17:
        do_one_tile(field,tile=i)
        i+=1
    
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
            tiles.append(int(argv[i]))
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

