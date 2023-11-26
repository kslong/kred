#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a diagnositc routine intended to allow one to
make plots of rows and columns in regions of images
that show breaks in the apparent flux


Command line usage (if any):

    usage: MakeXcut.py  -filt N662 -exptime 400 -size 200 -ymin 30 -ymax 100 xdir ra dec

    where

        -xfilt designates a fileter to use
        -exptime designates and expousre time to chose
        -size indicates the size in pixels of the xcuts
        -ymin and -ymax set the min and maxim y leve

    and 
        xdir is the directory in which to look for fits files
        ra and dec are the ra and dec of a position to inspect


Description:  

Primary routines:

    doit

Notes:
                                       
History:

231018 ksl Coding begun

'''

import os
from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from glob import glob

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from decimal import Decimal, ROUND_HALF_EVEN

def round2integer(arr):
    rounded_arr = np.vectorize(lambda x: int(Decimal(x).quantize(Decimal('1'), rounding=ROUND_HALF_EVEN)))(arr)
    return rounded_arr



def convert_coordinates(ra, dec, is_degrees=True):
    '''
    Convert coordiantes to degrees if
    necessary
    '''

    if is_degrees:
        return ra, dec
    else:
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg



def radec2pix(filename='c4d_221221_004943_ooi_N673_req_N24.fits',ra='01:02:38.6',dec='-71:59:08.23'):
    '''
    Open a file and see if the postion is inside the image.  If so incdicate the location and
    flux nearby.  Also collect some other information about the object
    '''
    # Load the FITS file
    hdulist = fits.open(filename)
    # hdulist.info()
    wcs = WCS(hdulist[1].header)
    try:
        xfilt=hdulist[0].header['FILTER']
        words=xfilt.split()
        xfilt=words[0]
        xtime=hdulist[0].header['EXPTIME']
    except:
        xfilt='Unknown'
        xtime==-99.

    # print(ra,dec)
    try:
        ra=eval(ra)
        dec=eval(dec)
    except:
        # print('what')
        xra, xdec =  convert_coordinates(ra,dec,False)
        ra=xra
        dec=xdec


    # print(ra,dec)
    # Create a SkyCoord object with the target coordinates
    target_coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    

    # print(target_coord)
    # Convert RA and Dec to pixel coordinates
    pixel_x, pixel_y = wcs.world_to_pixel(target_coord)
    # print(pixel_x,pixel_y)
    pixel_x=round2integer(pixel_x)
    pixel_y=round2integer(pixel_y)

    
    # Check if the pixel coordinates are within the image bounds
    if not (0 <= pixel_x < hdulist[1].header['NAXIS1'] and 0 <= pixel_y < hdulist[1].header['NAXIS2']):
        pixel_x=pixel_y=value=-99
    else:
        half=5
        ymin=pixel_y-half-1
        ymax=ymin+2*half
        xmin=pixel_x-half-1
        xmax=xmin+2*half
        if ymin<0:
            ymin=0
        if xmin<0:
            xmin=0
        if ymax>=hdulist[1].header['NAXIS2']:
            ymax=hdulist[1].header['NAXIS2']-1
        if xmax>=hdulist[1].header['NAXIS1']:
            xmax=hdulist[1].header['NAXIS1']-1
        xarray=hdulist[1].data[ymin:ymax,xmin:xmax]
        value=np.median(xarray)
    # Close the FITS file
    hdulist.close()
    return pixel_x,pixel_y,value,xfilt,xtime



def find_overlaps(xdir='.',ra='01:02:38.6',dec='-71:59:08.23',size=1000):
    
    glob_string='%s/*fits' % (xdir)
    files=glob(glob_string)
    if len(files)==0:
        print ('Foound no overlaps at {ra} {dec}')
        return []
    
    xfile=[]
    x=[]
    y=[]
    xmed=[]
    xfilt=[]
    xtime=[]
    for one_file in files:
        xx,yy,value,xxfilt,xxtime=radec2pix(one_file,ra,dec)
        if xx>0:
            x.append(xx)
            y.append(yy)
            xmed.append(value)
            xfile.append(one_file)
            xfilt.append(xxfilt)
            xtime.append(xxtime)
    if len(xfile)==0:
        if isinstance(ra,float):
            print('Files, but none with valid pixels at %.6f %.6f' % (ra,dec))
        else:
            print('Files, but none with valid pixels at %s %s' % (ra,dec))
        return []
    else:
        print('Found %d files that overlap' % (len(xfile)))
        xtab=Table([xfile,xfilt,xtime,x,y,xmed],names=['Filename','Filter','Exptime','x','y','Counts'])
        return xtab
        

def do_xcut_plot(ztab,xfilt='N662',exptime=400.,size=100,ymin=30,ymax=100):
    
    if xfilt!='':
        zztab=ztab[ztab['Filter']==xfilt]
        if len(zztab)==0:
            print('Files with overlaps and filter %s not found' % xfilt)
            return 0
    else:
        zztab=ztab
    
    if exptime>0:
        qtab=zztab[zztab['Exptime']==exptime]
        if len(qtab)==0:
            print('Files with overlaps and exptime %s not found' % xfilt)
            return 0
    else:
        qtab=zztab
    
    base=np.arange(size)-int(0.5*size)
    row_sum=np.zeros_like(base)
    row_num=np.zeros_like(base)
    col_sum=np.zeros_like(base)
    col_num=np.zeros_like(base)    
    
    plt.figure(1,(12,6))
    plt.clf()
    for one_row in qtab:
        f=fits.open(one_row['Filename'])
        row=f[1].data[one_row['y'],]
        xcol=np.arange(len(row))-one_row['x']
        xmin=int(one_row['x']-0.5*size)
        xmax=int(one_row['x']+0.5*size)
        offset=0
        if xmin<0:
            offset=-xmin
            xmin=0
        if xmax>=f[1].header['NAXIS1']:
            xmax=f[1].header['NAXIS1']-1

        row=row[xmin:xmax]
        xcol=xcol[xmin:xmax]
        # xcole contains the array elements offset by half size
        
        istart=offset
        istop=istart+len(row)
        row_num[istart:istop]+=1
        
        
        row_sum[istart:istop]=np.add(row_sum[istart:istop],row)
        
        plt.subplot(2,1,1)
        plt.plot(xcol,row,'.')
        
        
        col=f[1].data[:,one_row['x']]
        xrow=np.arange(len(col))-one_row['y']
        
        xmin=int(one_row['y']-0.5*size)
        xmax=int(one_row['y']+0.5*size)
        offset=0
        if xmin<0:
            offset=-xmin
            xmin=0
        if xmax>=f[1].header['NAXIS2']:
            xmax=f[1].header['NAXIS2']-1
        
        
        xrow=xrow[xmin:xmax]
        col=col[xmin:xmax]
        
        istart=offset
        istop=istart+len(col)
        
        col_num[istart:istop]+=1
        col_sum[istart:istop]=np.add(col_sum[istart:istop],col)
        
        plt.subplot(2,1,2)
        plt.plot(xrow,col,'.')
        
    row_ave=row_sum/row_num
    plt.subplot(2,1,1)
    plt.plot(base,row_ave,'k.')
    plt.xlabel('Column')
    plt.ylabel('Counts')
    plt.xlim(-size/2,size/2)
    if ymin!=0 or ymax!=0:
        plt.ylim(ymin,ymax)
    plt.tight_layout()
    
    col_ave=col_sum/col_num
    plt.subplot(2,1,2)
    plt.plot(base,col_ave,'k.')
    plt.xlabel('Row')
    plt.ylabel('Counts')
    plt.xlim(-size/2,size/2)
    if ymin!=0 or ymax!=0:
        plt.ylim(ymin,ymax)
    plt.tight_layout()
    plt.tight_layout()

    return 1

def doit(ra='01:02:38.6',dec='-71:59:08.23',xdir='.',xfilt='N662',exptime=400,size=200,ymin=30,ymax=100):
    '''
    Run the script
    '''
    
    xtab=find_overlaps(xdir,ra,dec,size)

    if len(xtab)==0:
        print('No files with this position')
        return

    xtab.write('xcut_tab.txt',format='ascii.fixed_width_two_line',overwrite=True)
    

    ireturn=do_xcut_plot(xtab,xfilt,exptime,size,ymin,ymax)

    if ireturn:
        plt.savefig('xcut_%s_%s.png' % (xfilt,exptime))
    else:
        filters=np.unique(xtab['Filter'])
        for one_filter in filters:
            foo=xtab[xtab['Filter']==one_filter]
            exptimes,counts=np.unique(foo['Exptime'],return_counts=True)
            j=0
            while j<len(exptimes):
                print('%12s %8.1f %3d' % (one_filter,exptimes[j],counts[j]))

                j+=1






def steer(argv):
    '''
    MakeXcut.py  -filt N662 -exptime 400 -size 200 -ymin 30 -ymax 100 xdir ra dec
    '''

    xfilt='N662'
    exptime=400
    size=200
    ymin=0
    ymax=100
    xdir=''
    ra=''
    dec=''  

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-filt':
            i+=1
            xfilt=argv[i]
        elif argv[i]=='-exptime':
            i+=1
            exptime=eval(argv[i])
        elif argv[i]=='-size':
            i+=1
            size=eval(argv[i])
            size=int(size)
        elif argv[i]=='-ymin':
            i+=1
            ymin=eval(argv[i])
        elif argv[i]=='-ymax':
            i+=1
            ymax=eval(argv[i])
        elif argv[i][0]=='-' and xdir=='':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif xdir=='':
            xdir=argv[i]
        elif ra=='':
            ra=argv[i]
        elif dec=='':
            dec=argv[i]
        i+=1


    if dec=='':
        print(__doc__)
        print('Error: Not enough inputs consumed ',argv)
        return
    if os.path.isdir(xdir)==False:
        print('Error: directory (%s) not found' % xdir)
        return



    doit(ra,dec,xdir,xfilt,exptime,size,ymin,ymax)

    return

            








# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)


