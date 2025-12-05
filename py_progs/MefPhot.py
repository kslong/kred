#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Carry out forced photomentry on (multi-extension)
files and produce tables of the outputs


Command line usage (if any):

    usage MefPhot.py [-h] [-np 8] [-r 6] [-b 8 12][-out root] file1 file2 ...

Description:  

    where
        -h prints out this documentaiton and exits
        -np 8 sets the number of processors to use (default 8)
        -r 6  the radius in pixels to be used for extraction of the stellar flux
        -b 8 12 the inner and outer radius of the background to be used
        -out root changes the root name of the output table

    As normally  run, the routine reads each image extension in one or
    more fits files, and then finds stars in the Gaia catalog
    of the region, and does carries out forced photometry of the
    gaia stars


    A table containing the results from the forced photmetry is produced
    and normally placed in the directory TabPhot


Primary routines:

    do_one where the work is actually done
    do_many when multiprocessing is used, this distributes the work

Notes:

    The Gaia data should be contained in a file Gaia_MagClouds.fits, located either locally
    or in kred/xdata
                                       
History:

251128 ksl Coding begun
251205 ksl This has been tested with version 2.3.0 of photutils; earlier versions seemed to "hang".  
    A typical MEF files takes or order 8 minutes on an M1 Mac.

'''

import sys
import os
import numpy as np
from astropy.io import fits,ascii
from photutils.detection import DAOStarFinder

from astropy.stats import mad_std
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from astropy.stats import SigmaClip

import matplotlib.pyplot as plt
from astropy.wcs import WCS

import matplotlib.pyplot as plt
from astropy.table import Table,join,hstack

from astropy.coordinates import SkyCoord
import astropy.units as u
import timeit
import time

from astropy.time import Time

import pathlib
import os.path as path
import requests

from scipy.spatial import KDTree
import numpy as np
from astropy.table import Table
from astropy.wcs import NoConvergence
from astropy.wcs._wcs import InvalidCoordinateError

import time


# import GaiaCat
from GaiaCat import get_gaia
import ImageSum
import numpy as np
import matplotlib as plt
from astropy.io import fits
from astropy.table import Table,join, vstack


def random_rows(tab, nrows, seed=None):
    """
    Randomly select rows from an Astropy Table without duplicates.

    Parameters
    ----------
    tab : astropy.table.Table
        Input table.
    nrows : int
        Number of rows to randomly select (must be ≤ len(tab)).
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    subtab : astropy.table.Table
        Table containing the randomly selected rows.
    """
    if nrows > len(tab):
        print("Requested more rows than available in table")
        return tab

    rng = np.random.default_rng(seed)
    indices = rng.choice(len(tab), size=nrows, replace=False)
    return tab[indices]


def read_table(filename):
    '''
    This is a generic routine to try to read a table
    in fits or ascii format.  It is intended to accommodate
    several different types of formats.
    '''

    #print('XXXX - filename ',filename)

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



def do_forced_photometry(filename='LMC_c48_T08.r.t060.fits',image_ext=1,object_file='objects.txt',nrows_max=-1,rstar=4,b_in=4,b_out=8):
    '''
    Do forced photometry based on ra and decs, where the object file contains a set of source positions

    where:
        filename is a file with one or more image extensions.
        image_ext is the extension to be analyised
        nrows_rows limits the total number of object for forced photometry to a value, all
            if -1
        rstar is the aperure readius in pixles used for source extraction
        b_in and b_out define the size of the background annulus

    This is set up primarily for working with the files created from the Gaia catalog
    n IIxxxxxxxeveeeeee'''
    
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
        phot_table[col].info.format = '%.8g'
        
    pos = image_wcs.pixel_to_world(phot_table['xcenter'], phot_table['ycenter'])

    phot_table = hstack([phot_table, sources])

    if 'Source_name' not in phot_table.colnames:
        names = []
        for one in phot_table:
            names.append('x%05d' % one['id'])
        phot_table['Source_name'] = names    

    return phot_table


def do_one(filename='foo.fits',outroot='',nrows_max=-1,rstar=4,b_in=4,b_out=8):
    '''
    The driving routine for forced photmetry of 
    a single file.  

    where filename is the name of the file and
        outroot causes the name of the
        output table  to be given by

        TapPhot/outroot.fits

        Otherwise portions of the file name are
        used.

    This calls the forced
    photmetry routine for each image extension
    in the file.
    '''

    try:
        x=fits.open(filename)
    except:
        print('Could not locate %s' % filename)
        raise IOError

    print('do_one: Starting  %s with radius %.1f and annulus %.1f %.1f '  % (filename,rstar,b_in,b_out))

    xexptime=x['PRIMARY'].header['EXPTIME']
    try:
        xfilter=x['PRIMARY'].header['FILTER']
    except:
        words=filename.split('.')
        xfilter=words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter,filename))
        
    image_extensions=ImageSum.list_image_extensions(filename)

    phot_tables=[]

    i=0
    for one_extension in np.array(image_extensions['EXT']):
        info=ImageSum.get_image_center_and_size_from_header(x[one_extension].header)
        ra=info['center_ra']
        dec=info['center_dec']
        width=info['width_deg']
        height=info['height_deg']
        # print(ra,dec,width,height)
        size=np.sqrt(width*width+height*height)/2.
        # gaia_file=GaiaCat.get_gaia(ra,dec,size)
        gaia_file=get_gaia(ra,dec,size)
        # print(gaia_file)
        phot_table=do_forced_photometry(filename,one_extension,gaia_file,nrows_max,rstar,b_in,b_out)
        phot_table['EXT']=one_extension
        phot_table['CCD']=image_extensions['NAME'][i]
        phot_tables.append(phot_table)
        i+=1

    phot=vstack(phot_tables,metadata_conflicts='silent')
    phot['Filter']=xfilter
    phot['Exptime']=xexptime
    phot['Filename']=filename

        

    os.makedirs('TabPhot',exist_ok=True)
    if outroot=='':
        outroot=filename.split('/')[-1]
        outroot=outroot.replace('.fz','')
        outroot=outroot.replace('.fits','')
    outfile='TabPhot/%s.fits' % outroot

    # Before writing store the metadata
    now = Time.now()               # Current time in UTCi  
    timestamp = now.isot 
    phot.meta['DATE'] = timestamp
    phot.meta['FILE']=filename
    phot.meta['RADIUS']=rstar
    phot.meta['B_IN']=b_in
    phot.meta['B_OUT']=b_out

    phot.write(outfile,format='fits',overwrite=True)
    print('do_one: Completed %s and written to %s' % (filename,outfile))
    return phot

from multiprocessing import Pool
from tqdm import tqdm
import multiprocessing as mp
import traceback


def _safe_do_one_with_index(args):
    """Wrapper that unpacks (index, filename) and calls do_one with numbered outroot"""
    index, filename, outroot, nrows_max, rstar, b_in, b_out = args
    try:
        # Create numbered outroot if outroot is specified
        if outroot != '':
            numbered_outroot = f"{outroot}_{index:03d}"
        else:
            numbered_outroot = ''

        do_one(filename, outroot=numbered_outroot, nrows_max=nrows_max,
               rstar=rstar, b_in=b_in, b_out=b_out)
        return (filename, True, None)
    except Exception as e:
        # Capture the full error message and traceback
        error_msg = f"{type(e).__name__}: {str(e)}"
        tb = traceback.format_exc()
        return (filename, False, error_msg, tb)


def do_many(filenames, outroot='', nrows_max=-1, rstar=4, b_in=4, b_out=8,
            n_processes=None, logfile=None, verbose_errors=False):
    """
    Parallelize aperture photometry across multiple files.

    Parameters:
    -----------
    filenames : list
        List of image files to process
    outroot : str
        Output root directory/prefix. If not empty, each file gets outroot_NNN
    nrows_max : int
        Maximum number of rows to process per extension (-1 for all)
    rstar : float
        Aperture radius in pixels
    b_in : float
        Inner radius of background annulus in pixels
    b_out : float
        Outer radius of background annulus in pixels
    n_processes : int, optional
        Number of processes to use. If None, uses CPU count - 1
    logfile : str, optional
        Path to write failed filenames and errors. If None, only prints to screen.
    verbose_errors : bool
        If True, print full tracebacks for errors

    Returns:
    --------
    failed_files : list of tuples
        List of (filename, error_message) for files that failed
    """
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)

    # Create list of (index, filename, params...) tuples
    args_list = [(i, fname, outroot, nrows_max, rstar, b_in, b_out)
                 for i, fname in enumerate(filenames)]

    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap(_safe_do_one_with_index, args_list),
                           total=len(filenames),
                           desc="Processing images"))

    # Collect failed files with error messages
    failed_files = []
    for result in results:
        if len(result) == 3:  # Success case
            fname, success, _ = result
        else:  # Failure case
            fname, success, error_msg, tb = result
            if not success:
                failed_files.append((fname, error_msg, tb))

    if failed_files:
        print(f"\n{len(failed_files)}/{len(filenames)} files failed to process")
        print("\nErrors:")

        for fname, error_msg, tb in failed_files:
            print(f"\n  {fname}:")
            print(f"    {error_msg}")
            if verbose_errors:
                print("    Full traceback:")
                for line in tb.split('\n'):
                    print(f"      {line}")

        if logfile:
            with open(logfile, 'w') as f:
                f.write(f"{len(failed_files)}/{len(filenames)} files failed\n\n")
                for fname, error_msg, tb in failed_files:
                    f.write(f"{fname}\n")
                    f.write(f"  Error: {error_msg}\n")
                    if verbose_errors:
                        f.write(f"  Traceback:\n")
                        for line in tb.split('\n'):
                            f.write(f"    {line}\n")
                    f.write("\n")
            print(f"\nFailed filenames and errors written to {logfile}")
    else:
        print(f"\nSuccessfully processed all {len(filenames)} files")

    return failed_files

def steer(argv):
    '''
    usage MefPhot.py [-h] [-np 8] [-r 6] [-b 8 12][-out root] file1 file2 ...
    '''

    filenames=[]
    np=8
    root=''
    rstar=6

    nrows_max=-1

    b_in=8
    b_out=12

    
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:3]=='-np':
            i+=1
            np=int(argv[i])
        elif argv[i][:4]=='-out':
            i+=1
            root=argv[i]
        elif argv[i][:2]=='-r':
            i+=1
            rstar=float(argv[i])
            b_in=rstar+1
            b_out=rstar+4
        elif argv[i][:2]=='-b':
            i+=1
            b_in=float(argv[i])
            i+=1
            b_out=float(argv[i])
        elif argv[i][0]=='-':
            print('Error: unknow switch:',argv)
            return
        elif argv[i].count('.txt') or argv[i].count('tab'):
            xtab=ascii.read(argv[i])
            filenames=xtab['Filename']
        else:
            filenames.append(argv[i])
        i+=1

    print('Starting with %d filenames and rstar of %.1f and background annnulus of %.1f %.1f' % (len(filenames),rstar,b_in,b_out))

    if rstar > b_in or b_in > b_out:
        print('UNPHYSICAL limits for photometry')
        return

    if len(filenames)==1 or np<2:
        for one_file in filenames:
            do_one(filename=one_file,outroot=root,nrows_max=nrows_max,rstar=rstar,b_in=b_in,b_out=b_out)
        return

    do_many(filenames, outroot=root, nrows_max=nrows_max, rstar=rstar, b_in=b_in, b_out=b_out, 
            n_processes=np, logfile=None)
    return






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)

