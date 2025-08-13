#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

This routine simple allows one to sumarize some information
from all of the fits files in a particular part of the
kred directory structure, such as DECam_SUB2. 

It is intended to facillitate locating images that 
include a given RA and Dec on the sky

Command line usage (if any):

    usage: ImageSum.py [-h] [-out whatever] dirname 

    where:

        -h prints this documentation and exits
        -out whatever changes the name of the output
            file that is created.  Without this
            the name of the output file is based
            on dirname
        dir the directory which will be searched


Description:  

    The routine simply searches for fits files in
    a directory and all of its subdirectories 
    and produces a listing of all of the files,
    along with certain information derived from
    the header


Primary routines:

   steer - directs the routine 
   table_create - the main routine 

Notes:
                                       
History:

250813 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii,fits
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.wcs import WCS

def get_image_center_and_size_from_header(header):
    """
    Calculate the center coordinates (RA, Dec) and size in degrees 
    of a FITS image from its header containing WCS information.
    
    Parameters:
    -----------
    header : astropy.io.fits.Header
        FITS header containing WCS keywords and NAXIS information
        
    Returns:
    --------
    dict : Dictionary containing:
        - 'center_ra': RA of image center in degrees
        - 'center_dec': Dec of image center in degrees  
        - 'width_deg': Width of image in degrees
        - 'height_deg': Height of image in degrees
        - 'pixel_scale_ra_deg_per_pix': Pixel scale in RA direction (deg/pixel)
        - 'pixel_scale_dec_deg_per_pix': Pixel scale in Dec direction (deg/pixel)
        - 'image_shape': Tuple of (ny, nx) image dimensions
        - 'wcs': WCS object for further coordinate transformations
    """
    
    # Get image dimensions from header
    # Handle both 2D and higher dimensional data
    if 'NAXIS1' in header and 'NAXIS2' in header:
        nx = header['NAXIS1']  # Width (columns)
        ny = header['NAXIS2']  # Height (rows)
    else:
        raise ValueError("Header must contain NAXIS1 and NAXIS2 keywords")
    
    # Create WCS object from header
    wcs = WCS(header)
    
    # If WCS has more than 2 dimensions, take celestial subset
    if wcs.naxis > 2:
        wcs = wcs.celestial
    
    # Calculate center pixel coordinates (0-indexed)
    center_x = (nx - 1) / 2.0
    center_y = (ny - 1) / 2.0
    
    # Convert center pixel to world coordinates
    center_ra, center_dec = wcs.pixel_to_world_values(center_x, center_y)
    
    # Calculate corner coordinates to determine image size
    corners_pixel = np.array([
        [0, 0],           # bottom-left
        [nx-1, 0],        # bottom-right  
        [nx-1, ny-1],     # top-right
        [0, ny-1]         # top-left
    ])
    
    # Convert corners to world coordinates
    corners_world = wcs.pixel_to_world_values(corners_pixel[:, 0], corners_pixel[:, 1])
    corner_ra = corners_world[0]
    corner_dec = corners_world[1]
    
    # Calculate image size in degrees
    ra_min, ra_max = np.min(corner_ra), np.max(corner_ra)
    dec_min, dec_max = np.min(corner_dec), np.max(corner_dec)
    
    # Handle RA wrapping around 0/360
    if ra_max - ra_min > 180:
        # Likely wrapped around, recalculate
        corner_ra_wrapped = corner_ra.copy()
        corner_ra_wrapped[corner_ra_wrapped > 180] -= 360
        ra_min, ra_max = np.min(corner_ra_wrapped), np.max(corner_ra_wrapped)
        width_deg = ra_max - ra_min
    else:
        width_deg = ra_max - ra_min
        
    height_deg = dec_max - dec_min
    
    # Calculate approximate pixel scale at center
    # Small offset to calculate derivative
    offset = 0.5
    ra1, dec1 = wcs.pixel_to_world_values(center_x - offset, center_y)
    ra2, dec2 = wcs.pixel_to_world_values(center_x + offset, center_y)
    ra3, dec3 = wcs.pixel_to_world_values(center_x, center_y - offset)
    ra4, dec4 = wcs.pixel_to_world_values(center_x, center_y + offset)
    
    # Pixel scale in degrees per pixel
    pixel_scale_ra = abs(ra2 - ra1) / (2 * offset)
    pixel_scale_dec = abs(dec4 - dec3) / (2 * offset)
    
    return {
        'center_ra': center_ra,
        'center_dec': center_dec,
        'width_deg': abs(width_deg),
        'height_deg': abs(height_deg),
        'pixel_scale_ra_deg_per_pix': pixel_scale_ra,
        'pixel_scale_dec_deg_per_pix': pixel_scale_dec,
        'image_shape': (ny, nx),
        'wcs': wcs
    }

def print_header_image_info(header, description="FITS Header"):
    """
    Convenience function to print image information from header in a readable format.
    
    Parameters:
    -----------
    header : astropy.io.fits.Header
        FITS header containing WCS and dimension information
    description : str
        Optional description for the output
    """
    info = get_image_center_and_size_from_header(header)
    
    print(f"Image Information: {description}")
    print("=" * 50)
    print(f"Center RA:  {info['center_ra']:.6f}°")
    print(f"Center Dec: {info['center_dec']:.6f}°")
    print(f"Image width:  {info['width_deg']:.6f}° ({info['width_deg']*3600:.2f} arcsec)")
    print(f"Image height: {info['height_deg']:.6f}° ({info['height_deg']*3600:.2f} arcsec)")
    print(f"Image shape: {info['image_shape']} pixels")
    print(f"Pixel scale RA:  {info['pixel_scale_ra_deg_per_pix']:.8f} deg/pix ({info['pixel_scale_ra_deg_per_pix']*3600:.3f} arcsec/pix)")
    print(f"Pixel scale Dec: {info['pixel_scale_dec_deg_per_pix']:.8f} deg/pix ({info['pixel_scale_dec_deg_per_pix']*3600:.3f} arcsec/pix)")



def table_create(xdir='DECam_SUB2',outname=''):
    files=glob('%s/**/*.fits' % xdir,recursive=True)
    source=[]
    xtype=[]
    ra=[]
    dec=[]
    height=[]
    width=[]
    for one_file in files:
        name=one_file.split('/')[-1]
        word=name.split('.')
        xtype.append(word[-2])
        header = fits.getheader(one_file, ext=0)
        source.append(header['Object'])
        info=get_image_center_and_size_from_header(header)
        ra.append(info['center_ra'])
        dec.append(info['center_dec'])
        width.append(info['width_deg'])
        height.append(info['height_deg'])

    dec=np.array(dec)
    height=np.array(height)
    width=np.array(width)
    width*=np.cos(dec/(180./np.pi))
    xtab=Table([source,xtype,ra,dec,width,height,files],names=['Source_name','Image_type','RA','Dec','width','height','filename'])
    xtab['RA'].format='.5f'
    xtab['Dec'].format='.5f'
    xtab['width'].format='.2f'
    xtab['height'].format='.2f'
    if outname=='':
        qdir=xdir.replace('/','-')
        outname='Image_Sum_%s.txt' % qdir
    xtab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote summary with %d lines to %s' % (len(xtab),outname))

    types,counts=np.unique(xtab['Image_type'],return_counts=True)
    i=0
    while i < len(types):
        print('%20s  %5d' % (types[i],counts[i]))
        i+=1
    
    return xtab
            

def steer(argv):
    xdir=''
    xout=''
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            xout=argv[i]
        elif argv[i][0]=='-':
            print('Error: Cannot parse command line: ',argv)
            return
        else:
            xdir=argv[i]
        i+=1

    if xdir=='':
        print('Error: now direcory provided: ',argv)

    if os.path.isdir(xdir)==False:
        print('Could not find directory:',xdir)

    table_create(xdir,xout)
    return












# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
