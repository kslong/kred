#!/usr/bin/env python
# coding: utf-8

"""StarFind - Star Identification for PSF Creation

Space Telescope Science Institute

Synopsis
--------

Identifies stars suitable for creating PSF functions from astronomical images.

Command Line Usage
------------------

::

    Usage: StarFind.py -h -out root file1 file2 etc.

Description
-----------

This routine identifies stars suitable for PSF (Point Spread Function) creation.

Version History
---------------
251210 ksl
    Coding begun

251222 ksl
    Added to kred

"""


import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt

# # Idenify stars for creatio of psfs
# 
# 251209 - This is a continuation of what I was dooing previousl



import os
import numpy as np
from astropy.io import fits,ascii
from astropy.stats import mad_std
from astropy.stats import SigmaClip
from astropy.wcs import WCS
from astropy.wcs import NoConvergence
from astropy.wcs._wcs import InvalidCoordinateError
from astropy.table import Table,join,hstack
from astropy.coordinates import SkyCoord
import astropy.units as u

from photutils.detection import DAOStarFinder
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats

import matplotlib.pyplot as plt


from scipy.spatial import KDTree

import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)


import pathlib
import os.path as path
import requests


import time
from http.client import IncompleteRead



XDIR=''  # Part of a directory name; used to isolate different runs of PhotCompare


def read_table(filename):
    """
    This is a generic routine to try to read a table
    in fits or ascii format.  It is intended to accommodate 
    several different types of formats.
    """

    print('XXXX - filename ',filename)

    if not os.path.isfile(filename):
        raise IOError ('read_table: %s does not appear to exist' % filename)

    try:
        xtable=Table.read(filename)
    except:
        try:
            xtable=ascii.read(filename)
        except:
            raise IOError('read_table: %s exist, but could not be read' % filename)
    return xtable

def get_objects_from_image(filename='LMC_c48_T08.r.t060.fits',outroot=''):
    
    print ('starting')
    try:
        x=fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

    print('get_photometry: Beginning photometry of %s' % filename)

    xexptime=x['PRIMARY'].header['EXPTIME']
    try:
        xfilter=x['PRIMARY'].header['FILTER']
    except:
        words=filename.split('.')
        xfilter=words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter,filename))

    tab_dir='./TabPhot%s' % XDIR

    os.makedirs(tab_dir,exist_ok=True)

    
    if outroot=='':
        words=filename.split('/')
        outroot=words[-1].replace('.fits','')
        
    # Allow for the data to be in the first of second image
    if x[0].data is not None:
        image_wcs=WCS(x[0].header)
        image=x[0].data
    elif x[1].data is not None:
        image_wcs=WCS(x[1].header)
        image=x[1].data

    # print(np.median(image))

    image-=np.median(image)
    
    bkg_sigma = mad_std(image)  

    daofind = DAOStarFinder(fwhm=4.0, threshold=3.0 * bkg_sigma)  

    sources = daofind(image)  
    # sources.info()

    pos=image_wcs.pixel_to_world(sources['xcentroid'],sources['ycentroid'])
    sources['RA']=pos.ra.degree
    sources['Dec']=pos.dec.degree



    for col in sources.colnames:  
        if sources[col].dtype.kind=='f':
            sources[col].info.format = '%.8g'  # for consistent table output

    sources['Source_name']='Unknown'
    i=0
    while i<len(sources):
        sources['Source_name'][i]='X%06d' % (i+1)
        i+=1

    

    outname='%s/%s_sources.fits' % (tab_dir,outroot)
    sources.write(outname,format='fits',overwrite=True)

    return outname



def do_forced_photometry(filename='LMC_c48_T08.r.t060.fits',image_ext=0,object_file='objects.txt',nrows_max=-1,rstar=6,b_in=8,b_out=12,add_psf_metrics=False):
    """Do forced photometry based on RA/Dec positions from object file.

    Parameters
    ----------
    filename : str, optional
        FITS file with one or more image extensions. Default is 'LMC_c48_T08.r.t060.fits'.
    image_ext : int, optional
        Extension to be analyzed. Default is 0.
    object_file : str, optional
        File containing source positions. Default is 'objects.txt'.
    nrows_max : int, optional
        Maximum number of objects for forced photometry. If -1, processes all. Default is -1.
    rstar : float, optional
        Aperture radius in pixels for source extraction. Default is 6.
    b_in : float, optional
        Inner radius of background annulus in pixels. Default is 8.
    b_out : float, optional
        Outer radius of background annulus in pixels. Default is 12.
    add_psf_metrics : bool, optional
        If True, adds columns useful for PSF star selection. Default is False.

    Returns
    -------
    phot_table : astropy.table.Table
        Photometry results table.

    Notes
    -----
    This function returns the table but does not write it to disk.
    
    """


    try:
        x=fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

    image_wcs=WCS(x[image_ext].header)
    image=x[image_ext].data

    # Don't subtract median here - let local background handle it
    image_mask = (image == 0) | ~np.isfinite(image)
    NAXIS1=x[image_ext].header['NAXIS1']
    NAXIS2=x[image_ext].header['NAXIS2']

    sources=read_table(object_file)
    if 'G' in sources.colnames:
        good = ~sources['R'].mask
        sources=sources[good]

    print('XXX starting forced photometery with %d sources' % (len(sources)))

    coords = SkyCoord(ra=sources['RA']*u.deg, dec=sources['Dec']*u.deg)

    # Initialize with NaNs
    sources['xcentroid'] = np.nan
    sources['ycentroid'] = np.nan

    try:
        x_pix, y_pix = image_wcs.world_to_pixel(coords)
        sources['xcentroid'] = x_pix
        sources['ycentroid'] = y_pix
    except (NoConvergence, InvalidCoordinateError) as e:
        if isinstance(e, NoConvergence):
            sources['xcentroid'] = e.best_solution[0]
            sources['ycentroid'] = e.best_solution[1]
            print(f"Warning: {len(e.divergent)} coordinates failed to converge")
        else:
            print(f"Warning: Severe coordinate transformation error - skipping bad coordinates")

    margin=2*b_out

    mask = (
        np.isfinite(sources['xcentroid']) &
        np.isfinite(sources['ycentroid']) &
        (sources['xcentroid'] >= margin) &
        (sources['xcentroid'] < NAXIS1 - margin) &
        (sources['ycentroid'] >= margin) &
        (sources['ycentroid'] < NAXIS2 - margin)
    )

    sources = sources[mask]

    center_coord = image_wcs.pixel_to_world(NAXIS1/2, NAXIS2/2)

    npossible=len(sources)

    if nrows_max>0 and len(sources)>nrows_max:
        sources=random_rows(sources, nrows=nrows_max, seed=None)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  

    apertures = CircularAperture(positions, r=rstar)  
    annulus_apertures = CircularAnnulus(positions, r_in=b_in, r_out=b_out)

    # Calculate background statistics with sigma clipping
    sigclip = SigmaClip(sigma=3, maxiters=10)
    bkg_stats = ApertureStats(image, annulus_apertures, sigma_clip=sigclip, mask=image_mask)

    # Get aperture statistics on ORIGINAL image (not background-subtracted)
    aper_stats = ApertureStats(image, apertures, sigma_clip=None, mask=image_mask)

    # Calculate background-subtracted flux
    bkg_mean = bkg_stats.mean  # median would be more robust
    n_aper_pixels = aper_stats.sum_aper_area.value
    total_background = bkg_mean * n_aper_pixels
    net = aper_stats.sum - total_background

    # Error estimation
    bkg_std_per_pixel = bkg_stats.std
    n_bkg_pixels = bkg_stats.sum_aper_area.value

    error = np.sqrt(
        np.abs(net) +  # Poisson from source
        n_aper_pixels * bkg_std_per_pixel**2 +  # Background noise in aperture
        n_aper_pixels**2 * bkg_std_per_pixel**2 / n_bkg_pixels  # Background estimation error
    )

    # Create output table
    phot_table = Table()
    phot_table['id'] = np.arange(len(positions))
    phot_table['xcenter'] = positions[:, 0]
    phot_table['ycenter'] = positions[:, 1]
    phot_table['Raw'] = aper_stats.sum 
    phot_table['Bkg'] = total_background
    phot_table['BkgMean'] = bkg_mean
    phot_table['BkgStd'] = bkg_std_per_pixel
    phot_table['Net'] = net
    phot_table['ErrNet'] = error

    # For FWHM, create a local background-subtracted version
    # This is more reliable than using the global image

    # For FWHM, create a local background-subtracted version
    phot_table['FWHM'] = np.nan
    phot_table['Eccentricity'] = np.nan

    for i in range(len(positions)):
        # Extract cutout around source
        x_int, y_int = int(positions[i, 0]), int(positions[i, 1])
        cutout_size = int(2 * b_out) + 10
        y_min = max(0, y_int - cutout_size)
        y_max = min(NAXIS2, y_int + cutout_size)
        x_min = max(0, x_int - cutout_size)
        x_max = min(NAXIS1, x_int + cutout_size)

        if y_max > y_min and x_max > x_min:
            cutout = image[y_min:y_max, x_min:x_max] - bkg_mean[i]
            cutout_pos = [(positions[i, 0] - x_min, positions[i, 1] - y_min)]
            cutout_aper = CircularAperture(cutout_pos, r=rstar)
            cutout_stats = ApertureStats(cutout, cutout_aper, sigma_clip=None)

            # Extract value without units
            if np.isfinite(cutout_stats.fwhm.value):
                phot_table['FWHM'][i] = cutout_stats.fwhm.value
            if np.isfinite(cutout_stats.eccentricity):
                phot_table['Eccentricity'][i] = cutout_stats.eccentricity

    phot_table['Max'] = aper_stats.max
    phot_table['Min'] = aper_stats.min

    # Add minimal PSF quality metrics if requested
    if add_psf_metrics:
        # Signal-to-noise ratio
        phot_table['SNR'] = np.abs(phot_table['Net']) / phot_table['ErrNet']

        # Light concentration (peak/mean flux density)
        mean_flux_density = phot_table['Net'] / n_aper_pixels
        phot_table['Concentration'] = phot_table['Max'] / mean_flux_density

        # Background contamination flag (normalized to median)
        median_bkg_std = np.nanmedian(phot_table['BkgStd'])
        phot_table['BkgContam'] = phot_table['BkgStd'] / median_bkg_std

    # Magnitude calculation
    phot_table['phot_mag'] = 28 - 2.5 * np.log10(np.abs(phot_table['Net']))
    phot_table['phot_mag_raw'] = 28 - 2.5 * np.log10(np.abs(phot_table['Raw']))

    phot_table['phot_mag'] = np.select([phot_table['Net'] > 0], 
                                       [phot_table['phot_mag']], 
                                       default=-phot_table['phot_mag'])
    phot_table['phot_mag_raw'] = np.select([phot_table['Net'] > 0], 
                                           [phot_table['phot_mag_raw']], 
                                           default=-phot_table['phot_mag_raw'])

    for col in phot_table.colnames:  
        if phot_table[col].dtype.kind=='f':
            phot_table[col].info.format = '%.8g'

    pos = image_wcs.pixel_to_world(phot_table['xcenter'], phot_table['ycenter'])

    phot_table = hstack([phot_table, sources])

    if 'Source_name' not in phot_table.colnames:
        names = []
        for one in phot_table:
            names.append('x%05d' % one['id'])
        phot_table['Source_name'] = names    

    return phot_table


def do_one(filename, root='test'):
    """
    Find stars and perform photometry on an image.

    Parameters
    ----------
    filename : str
        Path to FITS image file.
    root : str, optional
        Output filename root. If empty, derived from filename.

    Notes
    -----
    Writes {root}_all_stars.fits containing all detected stars with photometry.
    Use PsfBuild to select PSF stars and build PSF models.
    """

    if root == '':
        root = filename.split('/')[-1]
        root = root.replace('.fits', '')
        root = root.replace('.fz', '')

    source_file = get_objects_from_image(filename=filename, outroot=root)
    phot = do_forced_photometry(filename, object_file=source_file, add_psf_metrics=True)

    all_out = '%s_all_stars.fits' % root
    print('Writing all stars: %s' % all_out)
    phot.write(all_out, format='fits', overwrite=True)




def steer(argv):
    '''
    This is generally just a steering routine

    Usage: StarFind -h -out root file1 file2 etc.
    '''

    filenames=[]
    root=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            root=argv[i]
        elif argv[i][0]=='-':
            print('Error: Could not intepret commands: ',argv)
            return
        else:
            filenames.append(argv[i])
        i+=1


    
    i=0
    for one_file in filenames:
        xroot=root
        if len(filenames)>1 and root != '':
            xroot='%s%02d' % (root,i+1)
        do_one(one_file,xroot)
        i+=1



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
