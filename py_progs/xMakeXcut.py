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

from astropy.wcs import WCS

def get_pixel_scale(filename):
    '''
    Return the pixel scale from and image 
    that is presumed to be in extension 1
    '''

    # Open the FITS file
    hdulist = fits.open(filename)

    # Assuming the WCS information is in the header of the second HDU (extension 1)
    wcs = WCS(hdulist[1].header)

    # Get the center coordinates of the image
    image_center = (hdulist[1].data.shape[1] / 2, hdulist[1].data.shape[0] / 2)

    # Convert the center pixel coordinates to world coordinates
    world_center = wcs.pixel_to_world(*image_center)

    # Get the coordinates of an adjacent pixel
    adjacent_pixel = (image_center[0] + 1, image_center[1])

    # Convert the adjacent pixel coordinates to world coordinates
    world_adjacent = wcs.pixel_to_world(*adjacent_pixel)

    # Calculate the pixel scale
    delta_x = world_adjacent.ra.deg - world_center.ra.deg
    delta_y = world_adjacent.dec.deg - world_center.dec.deg

    pixel_scale = ((delta_x**2 + delta_y**2)**0.5) * 3600  # Convert to arcseconds

    # Close the FITS file
    hdulist.close()

    return pixel_scal

def round2integer(arr):
    rounded_arr = np.vectorize(lambda x: int(Decimal(x).quantize(Decimal('1'), rounding=ROUND_HALF_EVEN)))(arr)
    return rounded_arr



def radec2deg(ra, dec):
    '''
    Convert coordiantes to degrees if
    necessary
    '''

    try:
        ra=eval(ra)
        dec=eval(dec)
        return ra,deg
    except:
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg



def radec2pix(x,ra=15.,dec=-71):
    '''
    Open a file and see if the position is inside the image.  If so incdicate the location and
    flux nearby.  Also collect some other information about the object

    If the position is not in the image, then -99, -99  is returned for the pixel position
    '''


    # Load the FITS file
    hdulist = x
    wcs = WCS(hdulist[1].header)


    print(ra,dec)
    # Create a SkyCoord object with the target coordinates
    target_coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    

    pixel_x, pixel_y = wcs.world_to_pixel(target_coord)
    # print(pixel_x,pixel_y)
    pixel_x=round2integer(pixel_x)
    pixel_y=round2integer(pixel_y)

    
    # Check if the pixel coordinates are within the image bounds
    if not (0 <= pixel_x < hdulist[1].header['NAXIS1'] and 0 <= pixel_y < hdulist[1].header['NAXIS2']):
        pixel_x=pixel_y=value=-99
    return pixel_x,pixel_y


def  make_xcut_tab(xdir):

    glob_string='%s/*fits' % (xdir)
    files=glob(glob_string)
    if len(files)==0:
        print ('Found no fits files in %s' % xdir)
        return []
    
    xfile=[]
    xfilt=[]
    xtime=[]
    for one_file in files:
        x=fits.open(one_file)
        xfile.apend(one_file)
        filt=hdulist[0].header['FILTER']
        words=filt.split()
        xfilt.append(words[0])
        xtime.append(hdulist[0].header['EXPTIME'])

    xtab=Table([xfile,xfilt,xtime],names=['Filename','Filter','Exptime'])
    xtab_name='%s/xcut_tab.txt' % xdir
    xtab.write(xtab_name,forrmat='ascii.fixed_width_two_line',overwrite=True)
    return xtab
    



def do_xcut_plot(ztab_name,ra=10,dec=-70,xfilt='N662',exptime=400.,size=100,ymin=30,ymax=100):

    try:
        ztab=ascii.read(ztab_name)
    except:
        print('Could not read %s' % ztab_name)
        return 0
    
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
        x,y=radec2pix(f,ra,dec)
        row=f[1].data[y,]
        xcol=np.arange(len(row))-x
        xmin=int(x-0.5*size)
        xmax=int(y+0.5*size)
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
        
        
        col=f[1].data[:,x]
        xrow=np.arange(len(col))-one_row['y']
        
        xmin=int(y-0.5*size)
        xmax=int(y+0.5*size)
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


    xtab_name='%s/xcut_tab.txt' % xdir
    if os.path.isfile(xtab_name)==False:
        make_xcut_tab(xdir)

    ra,dec=radec2deg(ra,dec)

    ireturn=do_xcut_plot(xtab_name,ra,dec,xfilt,exptime,size,ymin,ymax)

    if ireturn:
        plt.savefig('xcut_%s_%s.png' % (xfilt,exptime))
    else:
        xtab=ascii.read(xtab_name)
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


