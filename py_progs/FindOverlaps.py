#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Calculate a which images obtained with a filter overlap
with which others


Command line usage (if any):

    usage: FindOverlaps.py [-all] Field Tile1  Tile2 .....

    where -all implies that the ovelaps for all of the 
    tiles in afiled will be dound

Description:  

    Find the images that ovelap that have the same
    filter and the same exposure time in one or more tiles

    The routine assumes (for now) that we will not try
    to combine images with different exposure times,
    though this condition could be dropped

    The overlaps are found based on the RA, DEC corrners 
    of the immages.
    
    The routine outputs a single file that indicates
    which images overlap subject to thise conditions

    The routine DOES NOT calculate fluxes in the
    overlap regions, just that the ovelaps exist


Primary routines:

    do_one_tile

Notes:

    It is important to recognize that care has
    been taken to avoid double listing, that is
    to say in separate lines of the output file
    that 
    
    file_a matches file_b
    file_b matches file_a

    There is no significance to the order of
    filer_a and file_b.  
                                       
History:

230524 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import join, vstack



def get_files(field,tile,xfilter='Ha'):
    '''
    Idenfify all of the files associated
    with a specific filter in a tile
    '''

    det_file='Summary/%s_det.tab' % (field)
    try:
        det=ascii.read(det_file)
        filenames=[]
        for one in det:
            filenames.append('%s_%s.fits' % (one['Root'],one['EXTNAME']))
        det['Filename']=filenames
    except:
        print('Error could not locate %s ' % det_file)
        raise IOError 

    imsum_file='Summary/%s_%s.txt' % (field,tile)  
    try:
        imsum=ascii.read(imsum_file)
    except:
        print('Error: could not locate %s' % imsum_file)
        raise IOError

    xall=join(imsum,det,join_type='left')

    xtab=xall[xall['FILTER']==xfilter]

    if len(xtab)==False:
        print('No appropriate files were found for %s' % xfilter)
        return []

    return xtab




def get_overlap(row,xxtab):
    '''
    Get all of the images which have the
    same exposure time as the given row
    and also overlap with a given row 
    in the table

    Note eliminating images that have
    different exposure times implies
    that we do not intend to match backgrounds
    from images that have differnt times
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

    xtab=xtab[xtab['EXPTIME']==one_row['EXPTIME']]
    if len(xtab)==0:
        # print('Failed to find any images with the same exposre time ' % one_row['EXPTIME'])
            return []
    
    
    xtab=xtab[xtab['Dec_min']<dec_max]
    
    # print(len(xtab))

    if len(xtab)==0:
            # print('Failed to find any images with Dec less than %.5f ' % dec_max)
            return []
    
    xtab=xtab[xtab['Dec_max']>dec_min]
    if len(xtab)==0:
            # print('Failed to find any images with Dec greater than %.5f ' % dec_min)
            return []
    # print(len(xtab))
    
    xtab['RA_max'] = np.max([xtab['COR1RA1'],xtab['COR2RA1'],xtab['COR3RA1'],xtab['COR4RA1']],axis=0)
    xtab['RA_min'] = np.min([xtab['COR1RA1'],xtab['COR2RA1'],xtab['COR3RA1'],xtab['COR4RA1']],axis=0)

    xtab=xtab[xtab['RA_min']<ra_max]
    if len(xtab)==0:
            # print('Failed to find any images with RA  less than %.5f ' % ra_max)
            return []
    # print(len(xtab))

    # print(len(x))
    xtab=xtab[xtab['RA_max']>ra_min]
    if len(xtab)==0:
            # print('Failed to find any images with RA  greater than %.5f ' % ra_min)
            return []
    
    xtab['delta_dec']=np.fabs(xtab['CENDEC1']-one_row['CENDEC1'])
    xtab['delta_ra']=np.fabs(xtab['CENRA1']-one_row['CENRA1'])*np.cos(one_row['CENDEC1']/57.29578)
    
    xtab['delta']=np.sqrt(xtab['delta_dec']**2+xtab['delta_ra']**2)
    
    xtab.sort(['delta','delta_dec','delta_ra'])
    
    xtab['delta'].format=xtab['delta_dec'].format=xtab['delta_ra'].format='0.3f'
    
    # print(len(xtab))
    
    return(xtab)
    
    
def do_one_filter(field='LMC_c45',tile='T07',xfilter='Ha'):
    '''
    Find the overlaps for one filter in one tile
    '''
    xtab= get_files(field,tile,xfilter)

    if len(xtab)==0:
        return []

    i=0
    imax=len(xtab)-2
    while i<imax:
        z=get_overlap(i,xtab)
        if len(z)>0:
            xout=z['Filename','FILTER','EXPTIME','delta','delta_ra','delta_dec']
            xout['XFilename']=xtab['Filename'][i]
            if i==0:
                xx=xout.copy()
            else:
                xx=vstack([xx,xout])
        # print('Finished row',i)
        i+=1
    xx.sort(['EXPTIME','delta'])

    # outfile='Summary/%s_%s_%s_overlap.txt' % (field,tile,xfilter)
    # xx.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return xx

def do_one_tile(field='LMC_c45',tile='T07'):
    '''
    Find the overlaps for all of the filters in one tile
    '''

    print('\nStarting Field %s Tile %s' % (field,tile))
    
    xfilters=['N662','N673','r','N708']

    z=[]
    xlen=[]
    for one in xfilters:
        x=do_one_filter(field,tile,xfilter=one)
        xlen.append(len(x))
        if len(x)>0:
            z.append(x)
    all_overlaps=vstack(z)

    i=0
    while i<len(xfilters):
        print('For %5s there are %5d overlaps' % (xfilters[i],xlen[i]))
        i+=1


    outfile='Summary/%s_%s_overlap.txt' % (field,tile)
    all_overlaps.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s\n' % outfile)



                
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
        do_one_tile(field,one)


    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)
