#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Calculate a which images obtained with a filter overlap
with which others


Command line usage (if any):

    usage: FindOverlaps.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

230524 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import join, vstack



def get_files(field,tile,band='N673'):
    '''
    Idenfify all of the files associated
    with a specific filter in a tile
    '''

    det_file='Summary/%s_det.tab' % (field)
    try:
        det=ascii.read(det_file)
    except:
        print('Error could not locate %s ' % det_file)
        raise IOError 

    mef_file='Summary/%s_mef.tab' % (field)  
    try:
        mef=ascii.read(mef_file)
    except:
        print('Erorr: could not locate %s' % mef_file)
        raise IOError

    xall=join(det,mef,join_type='left')

    xtab=xall[xall['BAND']==band]

    if len(xtab)==False:
        print('No approapriate files were found for %s' % band)
        return []

    return xtab




def get_overlap(row,xxtab):
    '''
    Get all of the images which overlap
    with a given rwo in the table
    '''
    
    one_row=xxtab[row]
    try: 
        xtab=xxtab[row+1:]
    except:
        print('Got to end')
        return []
    dec_max=np.max([one_row['COR1DEC1'],one_row['COR2DEC1'],one_row['COR3DEC1'],one_row['COR4DEC1']],axis=0)
    dec_min=np.min([one_row['COR1DEC1'],one_row['COR2DEC1'],one_row['COR3DEC1'],one_row['COR4DEC1']],axis=0)
    ra_max=np.max([one_row['COR1RA1'],one_row['COR2RA1'],one_row['COR3RA1'],one_row['COR4RA1']],axis=0)
    ra_min=np.min([one_row['COR1RA1'],one_row['COR2RA1'],one_row['COR3RA1'],one_row['COR4RA1']],axis=0)
    
    xtab['Dec_max'] = np.max([xtab['COR1DEC1'],xtab['COR2DEC1'],xtab['COR3DEC1'],xtab['COR4DEC1']],axis=0)
    xtab['Dec_min'] = np.min([xtab['COR1DEC1'],xtab['COR2DEC1'],xtab['COR3DEC1'],xtab['COR4DEC1']],axis=0)
    
    xtab=xtab[xtab['Dec_min']<dec_max]
    
    # print(len(xtab))

    if len(xtab)==0:
            print('Failed to find any images with Dec less than %.5f ' % dec_max)
            return []
    
    xtab=xtab[xtab['Dec_max']>dec_min]
    if len(xtab)==0:
            print('Failed to find any images with Dec greater than %.5f ' % dec_min)
            return []
    # print(len(xtab))
    
    xtab['RA_max'] = np.max([xtab['COR1RA1'],xtab['COR2RA1'],xtab['COR3RA1'],xtab['COR4RA1']],axis=0)
    xtab['RA_min'] = np.min([xtab['COR1RA1'],xtab['COR2RA1'],xtab['COR3RA1'],xtab['COR4RA1']],axis=0)

    xtab=xtab[xtab['RA_min']<ra_max]
    if len(xtab)==0:
            print('Failed to find any images with RA  less than %.5f ' % ra_max)
            return []
    # print(len(xtab))

    # print(len(x))
    xtab=xtab[xtab['RA_max']>ra_min]
    if len(xtab)==0:
            print('Failed to find any images with RA  greater than %.5f ' % ra_min)
            return []
    
    xtab['delta_dec']=np.fabs(xtab['CENDEC1']-one_row['CENDEC1'])
    xtab['delta_ra']=np.fabs(xtab['CENRA1']-one_row['CENRA1'])*np.cos(one_row['CENDEC1']/57.29578)
    
    xtab['delta']=np.sqrt(xtab['delta_dec']**2+xtab['delta_ra']**2)
    
    xtab.sort(['delta','delta_dec','delta_ra'])
    
    xtab['delta'].format=xtab['delta_dec'].format=xtab['delta_ra'].format='0.3f'
    
    # print(len(xtab))
    
    return(xtab)
    
    
def do_one_tile(field,tile,band):
    xtab= get_files(field,tile,band)

    i=0
    imax=len(xtab)-2
    while i<imax:
        z=get_overlap(i,xtab)
        if len(z)>0:
            xout=z['Root','BAND','EXPTIME','delta','delta_ra','delta_dec']
            xout['XRoot']=xtab['Root'][i]
            if i==0:
                xx=xout.copy()
            else:
                xx=vstack([xx,xout])
        # print('Finished row',i)
        i+=1
    xx.sort(['delta'])

    outfile='Summary/%s_%s_%s_overlap.txt' % (field,tile,band)
    xx.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return xx

            
                
def steer(argv):
    '''
    This is just a steering routine for running swarp on one or more
    tiles from the command line

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
        tiles=[]
        i=1
        while i<17:
            tiles.append('T%02d' % i)
            i+=1

    for one in tiles:
        do_one_tile(field,one,'N673')


    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)
