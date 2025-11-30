#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Compare the brightness es of objects found in kred-processed images to Gaia
assuming our standard rescaling


Command line usage (if any):

    usage: PhotCompare.py -h -dir DECamSWARP2/SMC_c01 -nmax 30000 -forced -unforced file1 file2 ...

Description:  

    The routines processes one or more files comparing Gaia photometry
    and images produced with kred, and produces a figure which is stored 
    in Figs_phot.  (The xmatch between GAIA and the image is stored
    in TabPhot)

    There are two basic modes, one which is invoked with -dir, and one if that argument is
    not present

    If a directory is given, then all of the fits files in that or any subdirecotry are processed.  If
    this is the case then any specific files are ignored

    If one or more files are given then only those files are processed.

    the various switches are as follows:

    -h prints out this help and quites
    -dir causes all files in the directory named and any subdirectory to be processed.  The is 
        a basic assumption made that these images are swarped versions of the original data
    -nmax places a limit on the number of positions that will be used for forced photometry in the 
        GAIA catalog.  If nmax<0 all positions are processed
    -forced causes the progrm to used forced photometry (this is the default)
    -unforced in this case the routine searches for sources in the image, and then x-matches the
        postions to GAIA. This is largely a diagnostic mode which might become necessary if there
        are concerns about the relative astrometry between GAIA and our images.  The results of
        the seach of the image are stored in TabPhot




Primary routines:

    do_many

Notes:

    The routine retrieves if necessary Gaia catalog information for
    an image (or group of images), carrieds out aperture photometry
    on the images, and then x-correlates the results.  

    The most time-consuming part of the process is Gaia catalog
    retrieval, so this is only done once, if all of the files
    have the same centers and sizes. The GAIA catalogs are stored
    in a subdirectory GAIA.  



    (There are some of functions that are not in the end used.)

                                       
History:

240318 ksl Coding begun
240527 ksl Speed up the catalog matching.
251105 ksl Split finding sources in an image from doing photometry
            on the sources

'''


# # Compare Photometry from image to image and from image to Gaia



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
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)



from astroquery.gaia import Gaia

import pathlib
import os.path as path
import requests
from gaiaxpy import calibrate

from scipy.spatial import KDTree
import numpy as np
from astropy.table import Table
from astropy.wcs import NoConvergence
from astropy.wcs._wcs import InvalidCoordinateError

import time
from http.client import IncompleteRead

from kred import ImageSum
from kred import GaiaCat


XDIR=''  # Part of a directory name; used to isolate different runs of PhotCompare


def read_table(filename):
    '''
    This is a generic routine to try to read a table
    in fits or ascii format.  It is intended to accommodate 
    several different types of formats.
    '''

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

def old_unique_rows_within_tol(tab, tol=0.01):
    """
    Return unique rows from an Astropy table based on approximate
    equality of RA, Dec, and Size within a given tolerance (in degrees).

    Parameters
    ----------
    tab : astropy.table.Table
        Table containing columns 'RA', 'Dec', and 'Size' (in degrees).
    tol : float, optional
        Matching tolerance in degrees. Default is 0.01°.

    Returns
    -------
    unique_tab : astropy.table.Table
        New table containing one representative row per unique group.
    """
    # Stack RA, Dec, Size into a NumPy array
    data = np.vstack([tab['RA'], tab['Dec'], tab['Size']]).T

    # Initialize list of unique rows
    unique_indices = []
    used = np.zeros(len(data), dtype=bool)

    for i in range(len(data)):
        if used[i]:
            continue
        diff = np.abs(data - data[i])
        mask = np.all(diff < tol, axis=1)
        used[mask] = True
        unique_indices.append(i)

    return tab[unique_indices]



def unique_rows_within_tol(tab, tol=0.01):
    """
    Return unique rows from an Astropy table based on approximate
    equality of RA, Dec, and Size within a given tolerance (in degrees).

    Parameters
    ----------
    tab : astropy.table.Table
        Table containing columns 'RA', 'Dec', and 'Size' (in degrees).
    tol : float, optional
        Matching tolerance in degrees. Default is 0.01°.

    Returns
    -------
    unique_tab : astropy.table.Table
        New table containing one representative row per unique group.
    mapping : np.ndarray
        Array of length len(tab) where mapping[i] gives the index in
        unique_tab that row i of the original table maps to.
    """
    # Stack RA, Dec, Size into a NumPy array
    data = np.vstack([tab['RA'], tab['Dec'], tab['Size']]).T

    # Initialize list of unique rows and mapping array
    unique_indices = []
    mapping = np.full(len(data), -1, dtype=int)
    used = np.zeros(len(data), dtype=bool)

    for i in range(len(data)):
        if used[i]:
            continue

        # Find all rows within tolerance of row i
        diff = np.abs(data - data[i])
        mask = np.all(diff < tol, axis=1)

        # Mark them as used and map them to the current unique group
        used[mask] = True
        group_idx = len(unique_indices)
        mapping[mask] = group_idx

        # Add representative row
        unique_indices.append(i)

    return tab[unique_indices], mapping


# Example usage:
# unique_tab, mapping = unique_rows_within_tol(my_table, tol=0.01)
#
# To get back to original table rows:
# for i, row in enumerate(my_table):
#     unique_idx = mapping[i]
#     print(f"Row {i} maps to unique row {unique_idx}")
#
# To find all original rows that map to a specific unique row:
# unique_row_idx = 5
# original_indices = np.where(mapping == unique_row_idx)[0]

def xdo_fig(xtab,outroot):

    outdir='./Figs_phot%s' %  XDIR

    os.makedirs(outdir,exist_ok=True)
    plt.figure(1,(12,6))
    plt.clf()
    plt.subplot(1,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.plot(xtab['G'],xtab['phot_mag'],'.',alpha=.05)
    plt.plot(xtab['G'],-xtab['phot_mag'],'.',alpha=.05)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')
    plt.text(13,20,outroot)

    plt.ylim(11,24)
    plt.xlim(11,24) 



    plt.tight_layout()
    plt.subplot(1,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.plot(xtab['R'],xtab['phot_mag'],'.',alpha=.05)
    plt.plot(xtab['R'],-xtab['phot_mag'],'.',alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')
    plt.ylim(11,24)
    plt.xlim(11,24)  
    plt.tight_layout()
    plt.savefig('%s/%s.png' % (outdir,outroot))



def do_fig(xtab,outroot=''):

    outdir='./Figs_phot%s' %  XDIR

    os.makedirs(outdir,exist_ok=True)
    plt.figure(1,(9,8))
    plt.clf()
    plt.subplot(2,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['G'],xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['G'],-xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label('G-R')
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0) 
    else:
        plt.scatter(xtab['G'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['G'],-xtab['phot_mag'],marker='.',alpha=.05)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')

    plt.ylim(14,22)
    plt.xlim(14,22) 



    plt.tight_layout()
    plt.subplot(2,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    if 'G' in xtab.colnames:
        sc=plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        sc=plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
        cbar=plt.colorbar(sc)
        cbar.set_label('G-R')
        # Make colorbar solid (ignore scatter alpha)
        if hasattr(cbar, "solids") and cbar.solids is not None:
            cbar.solids.set_alpha(1.0) 
    else:
        plt.scatter(xtab['R'],xtab['phot_mag'],marker='.',alpha=.05)
        plt.scatter(xtab['R'],-xtab['phot_mag'],marker='.',alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[11,24],'k-')
    plt.ylim(14,22)
    plt.xlim(14,22)  



    plt.subplot(2,2,3)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']-xtab['G'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['G'],xtab['phot_mag']+xtab['G'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label('G-R')
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0) 
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')

    plt.ylim(-2,2) 
    plt.xlim(14,22) 

    plt.subplot(2,2,4)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    under=xtab[xtab['phot_mag']>0]
    # plt.text(16,1.5,'Under %d Over %d' % (len(under),len(xtab)-len(under)))
    sc=plt.scatter(xtab['R'],xtab['phot_mag']-xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    sc=plt.scatter(xtab['R'],xtab['phot_mag']+xtab['R'],marker='.',alpha=.01,c=xtab['G']-xtab['R'],cmap='plasma',vmin=-1,vmax=1)
    cbar=plt.colorbar(sc)
    cbar.set_label('G-R')
    # Make colorbar solid (ignore scatter alpha)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_alpha(1.0) 
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')
    plt.ylim(-2,2) 
    plt.xlim(14,22)  

    plt.suptitle(outroot)
    # OK now we can save
    plt.tight_layout()

    plt.savefig('%s/%s.png' % (outdir,outroot))


def do_fig_diff(xtab,outroot):

    outdir='./Figs_phot%s' %  XDIR
    os.makedirs(outdir,exist_ok=True) 


    plt.figure(1,(12,6))
    plt.clf()
    plt.subplot(1,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.plot(xtab['G'],xtab['phot_mag']-xtab['G'],'.',alpha=.01)
    plt.plot(xtab['G'],xtab['phot_mag']+xtab['G'],'.',alpha=.01)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')
    plt.text(13,2,outroot)

    plt.ylim(-5,5) 
    plt.xlim(11,22) 


    under=xtab[xtab['phot_mag']>0]


    plt.tight_layout()
    plt.subplot(1,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.text(13,4,'Under %d Over %d' % (len(under),len(xtab)-len(under)))
    plt.plot(xtab['R'],xtab['phot_mag']-xtab['R'],'.',alpha=.01)
    plt.plot(xtab['R'],xtab['phot_mag']+xtab['R'],'.',alpha=.01)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,24],[0,0],'k-')
    plt.ylim(-5,5) 
    plt.xlim(11,22)  
    plt.tight_layout()
    plt.savefig('%s/%s.png' % (outdir,outroot))

    
def get_objects_from_image(filename='LMC_c48_T08.r.t060.fits',outroot=''):
    
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
        sources[col].info.format = '%.8g'  # for consistent table output

    outname='%s/%s_sources.txt' % (tab_dir,outroot)
    sources.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    return outname


def locate_first_image_extension(xx):
    '''
    where xx is an already open fits file

    This is to deal with the fact that for the CCD images, we often have the imge in extenstion 1

    '''

    i=0
    while i<len(xx):
        if isinstance(xx[i], (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)) and xx[i].data is not None:
            return i
        i+=1
    return -1


def do_forced_photometry(filename='LMC_c48_T08.r.t060.fits',object_file='objects.txt',nrows_max=-1,outroot='',rstar=6,b_in=8,b_out=12):
    '''
    Do forced photometry based on ra and decs, where the object file contains a set of source positions


    NOTE - this version needs to be replaced by the version in MefPhot, but that one does not write out
    the data within ther routine and this needs to be fixed both for forced and unforced photmetry
    '''
    
    try:
        x=fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'


    image_ext=locate_first_image_extension(x)
    if image_ext<0:
        raise IOError('No Image extension in %s' % filename)


    # ra,dec,size_deg=get_size(filename)
    # print('calculated ',ra,dec,size_deg)

        
    image_wcs=WCS(x[image_ext].header)
    image=x[image_ext].data
    image-=np.median(image)
    image_mask = (image == 0) |  ~np.isfinite(image)
    NAXIS1=x[image_ext].header['NAXIS1']
    NAXIS2=x[image_ext].header['NAXIS2']


    xexptime=x['PRIMARY'].header['EXPTIME']
    try:
        xfilter=x['PRIMARY'].header['FILTER']
    except:
        words=filename.split('.')
        xfilter=words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter,filename))


    sources=read_table(object_file)
    if 'G' in sources.colnames:
        good = ~sources['R'].mask      # True where FLUX is NOT masked
        sources=sources[good]
        # sources['G'] = sources['G'].filled()


    coords = SkyCoord(ra=sources['RA']*u.deg, dec=sources['Dec']*u.deg)

    # print('There are %d sources' % len(sources))
    # print('The range in ra and dec is  :', np.min(sources['RA']),np.max(sources['RA']),np.min(sources['Dec']),np.max(sources['Dec']))
    # print('The range in ra and dec is  :', np.median(sources['RA']),np.average(sources['RA']),np.median(sources['Dec']),np.average(sources['Dec']))



    # Initialize with NaNs
    sources['xcentroid'] = np.nan
    sources['ycentroid'] = np.nan

    try:
        x, y = image_wcs.world_to_pixel(coords)
        sources['xcentroid'] = x
        sources['ycentroid'] = y
    except (NoConvergence, InvalidCoordinateError) as e:
        if isinstance(e, NoConvergence):
            sources['xcentroid'] = e.best_solution[0]
            sources['ycentroid'] = e.best_solution[1]
            print(f"Warning: {len(e.divergent)} coordinates failed to converge")
        else:
            print(f"Warning: Severe coordinate transformation error - skipping bad coordinates")
            # Set all to NaN, mask will filter them out

    mask = (
        np.isfinite(sources['xcentroid']) &
        np.isfinite(sources['ycentroid']) &
        (sources['xcentroid'] >= 0) &
        (sources['xcentroid'] < NAXIS1) &
        (sources['ycentroid'] >= 0) &
        (sources['ycentroid'] < NAXIS2)
    )

    sources = sources[mask]
    print(f"Kept {np.sum(mask)} sources on detector out of {len(mask)}")

    center_coord = image_wcs.pixel_to_world(NAXIS1/2, NAXIS2/2)
    # print(f"Calculated center: RA={center_coord.ra.deg:.6f}, Dec={center_coord.dec.deg:.6f}")
    # print(f"CRVAL from header: RA={image_wcs.wcs.crval[0]:.6f}, Dec={image_wcs.wcs.crval[1]:.6f}")
    # print(f"CRPIX from header: {image_wcs.wcs.crpix}")


    npossible=len(sources)

    if nrows_max>0 and len(sources)>nrows_max:
        sources=random_rows(sources, nrows=nrows_max, seed=None)

    print('Forced photometry of %d of %d possible sources' % (len(sources),npossible))
    
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  

    apertures = CircularAperture(positions, r=rstar)  
    annulus_apertures=CircularAnnulus(positions,r_in=b_in, r_out=b_out)

    phot_table = aperture_photometry(image, apertures)  
    aper_stats=ApertureStats(image,apertures,sigma_clip=None,mask=image_mask)
    sigclip=SigmaClip(sigma=3,maxiters=10)
    bkg_stats=ApertureStats(image,annulus_apertures,sigma_clip=sigclip,mask=image_mask)
    total_background=bkg_stats.mean*aper_stats.sum_aper_area.value
    net=aper_stats.sum - total_background

    # Error estimation from background
    bkg_std_per_pixel = bkg_stats.std  # std dev per pixel in annulus
    n_aper_pixels = aper_stats.sum_aper_area.value
    n_bkg_pixels = bkg_stats.sum_aper_area.value

    # Total error includes:
    # 1. Poisson noise from source (approximated by the net flux)
    # 2. Background noise in aperture
    # 3. Uncertainty in background estimate
    error = np.sqrt(
        np.abs(net) +  # Poisson from source (assumes gain=1, ADU=electrons)
        n_aper_pixels * bkg_std_per_pixel**2 +  # Background noise in aperture
        n_aper_pixels**2 * bkg_std_per_pixel**2 / n_bkg_pixels  # Background estimation error
        )

    phot_table['Raw']=aper_stats.sum 
    phot_table['Bkg']=total_background
    phot_table['Net']=net
    phot_table['ErrNet']=error

    # 28th mag is correct

    phot_table['phot_mag']= 28-2.5*np.log10(np.fabs(phot_table['Net']))
    phot_table['phot_mag_simple']= 28-2.5*np.log10(np.fabs(phot_table['aperture_sum']))

    phot_table['phot_mag']=np.select([phot_table['Net']>0],[phot_table['phot_mag']],default=-phot_table['phot_mag'])
    phot_table['phot_mag_simple']=np.select([phot_table['Net']>0],[phot_table['phot_mag_simple']],default=-phot_table['phot_mag_simple'])


    for col in phot_table.colnames:  

        phot_table[col].info.format = '%.8g'  # for consistent table output
        
    pos=image_wcs.pixel_to_world(phot_table['xcenter'],phot_table['ycenter'])
    names=[]
    for one in phot_table:
        names.append('x%05d' % one['id'])
    phot_table['Source_name']=names
    phot_table['RA']=pos.ra.degree
    phot_table['Dec']=pos.dec.degree
    phot_table['File']=outroot
    phot_table['Filter']=xfilter
    phot_table['Exptime']=xexptime



    tab_dir='./TabPhot%s' % XDIR

    os.makedirs(tab_dir,exist_ok=True)

    
    if outroot=='':
        words=filename.split('/')
        outroot=words[-1].replace('.fits','')
        
    
    outfile='%s/%s_phot.txt' % (tab_dir,outroot)
    phot_table.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d objects' % (outfile,len(phot_table)))
    return outfile





def do_photometry(filename='LMC_c48_T08.r.t060.fits',outroot='',rstar=6,b_in=8,b_out=12):
    '''
    Locate and measure fluxes from source in an image
    '''


    
    try:
        x=fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

        

    image_ext=locate_first_image_extension(x)
    if x<0:
        raise IOError('No Image extension in %s' % filename)

    image_wcs=WCS(x[image_ext].header)
    image=x[image_ext].data

    # print(np.median(image))

    image-=np.median(image)

    xexptime=x['PRIMARY'].header['EXPTIME']
    try:
        xfilter=x['PRIMARY'].header['FILTER']
    except:
        words=filename.split('.')
        xfilter=words[-3]
        print('Filter keyword is missing. Setting to %s for %s' % (xfilter,filename))

    try:
        sources=ascii.read(object_file)
    except:
        print('Error: do_photometry: could not read object file %s' % object_file)
        return 'Error'
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  

    apertures = CircularAperture(positions, r=rstar)  
    annulus_apertures=CircularAnnulus(positions,r_in=b_in, r_out=b_out)

    phot_table = aperture_photometry(image, apertures)  
    aper_stats=ApertureStats(image,apertures,sigma_clip=None)
    sigclip=SigmaClip(sigma=3,maxiters=10)
    bkg_stats=ApertureStats(image,annulus_apertures,sigma_clip=sigclip)
    total_background=bkg_stats.median*aper_stats.sum_aper_area.value
    net=aper_stats.sum - total_background

    phot_table['Raw']=aper_stats.sum 
    phot_table['Bkg']=total_background
    phot_table['Net']=net

    phot_table['phot_mag']= 27-2.5*np.log10(phot_table['Net'])
    phot_table['phot_mag_simple']= 27-2.5*np.log10(phot_table['aperture_sum'])



    for col in phot_table.colnames:  

        phot_table[col].info.format = '%.8g'  # for consistent table output
        
    pos=image_wcs.pixel_to_world(phot_table['xcenter'],phot_table['ycenter'])
    names=[]
    for one in phot_table:
        names.append('x%05d' % one['id'])
    phot_table['Source_name']=names
    phot_table['RA']=pos.ra.degree
    phot_table['Dec']=pos.dec.degree
    phot_table['File']=outroot
    phot_table['Filter']=xfilter
    phot_table['Exptime']=xexptime


    

    tab_dir='./TabPhot%s' % XDIR

    os.makedirs(tab_dir,exist_ok=True)

    
    if outroot=='':
        words=filename.split('/')
        outroot=words[-1].replace('.fits','')
        
    
    # print(phot_table)  
    outfile='%s/%s_phot.txt' % (tab_dir,outroot)
    phot_table.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d objects' % (outfile,len(phot_table)))
    return outfile




def find_closest_objects(table1_path, table2_path, max_sep=0.5):
    '''
    Find the objects with a given distance given two astropy tables.  The
    routine returns only the closest object that satisfies this criterion.


    240527 - this is a new version which useds KDTree
    '''
    # Read the two Astropy tables

    table1=read_table(table1_path)
    table2=read_table(table2_path)

    print('get_closest_objects: Beginning x-match of %s and %s' % (table1_path,table2_path))
    
    # Convert RA and Dec columns to SkyCoord objects
    coords1 = SkyCoord(ra=table1['RA'] * u.degree, dec=table1['Dec'] * u.degree)
    coords2 = SkyCoord(ra=table2['RA'] * u.degree, dec=table2['Dec'] * u.degree)

    # Convert SkyCoord to Cartesian coordinates
    cartesian_coords1 = np.array([coords1.cartesian.x.value, coords1.cartesian.y.value, coords1.cartesian.z.value]).T
    cartesian_coords2 = np.array([coords2.cartesian.x.value, coords2.cartesian.y.value, coords2.cartesian.z.value]).T

    # Build a KDTree for the second set of coordinates
    tree = KDTree(cartesian_coords2)



    # Query the KDTree for the closest neighbor in table2 for each object in table1
    distances, indices = tree.query(cartesian_coords1)

    # Extract the matching rows from table2. this resorts table2 in the order of table 1
    closest_matches = table2[indices]

        # Compute the separations
    closest_coords = SkyCoord(ra=closest_matches['RA']*u.degree, dec=closest_matches['Dec']*u.degree)
    separations = coords1.separation(closest_coords).arcsecond


    # Create a new table to store the closest objects and their separations
    table1['Sep']=separations
    del closest_matches['Source_name']
    del closest_matches['RA']
    del closest_matches['Dec']
    table1['Sep'].format='.3f'
    table1['RA'].format='.6f'
    table1['Dec'].foramt='6f'
    xtab=hstack([table1,closest_matches])
    
    xtab=xtab[xtab['Sep']<max_sep]

    
    print('Of %d objects in %s and %d objects in %s, found %d matches' % (len(table1),table1_path,len(table2),table2_path,len(xtab)))

    tab_dir='TabPhot%s' % XDIR
    
    if len(xtab):
        words=table1_path.split('/')
        one=words[-1].replace('.txt','')
        one=one.replace('.fits','')
        words=table2_path.split('/')
        two=words[-1].replace('.txt','')
        two=two.replace('.fits','')
        outfile='%s/%s_x_%s.txt' % (tab_dir,two,one)
        xtab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('Error: There are no objects that are closer thn %f arcsec' % max_sep)
        return []

    return xtab


def get_size(filename='LMC_c48_T08.r.t060.fits'):

    try:
        x=fits.open(filename)
    except:
        print('get_size: Could not open %s' % filename)
        raise IOError('get_size: Could not open %s' % filename)

    try:
        wcs = WCS(x[0].header)
        # Get the shape of the image
        naxis1 = x[0].header['NAXIS1']
        naxis2 = x[0].header['NAXIS2']
    except:
        try:
            wcs = WCS(x[1].header)
            # Get the shape of the image
            naxis1 = x[1].header['NAXIS1']
            naxis2 = x[1].header['NAXIS2']
        except:
            raise IOError('get_size: Could not get info for %s' % filename)

    # Calculate the pixel coordinates of the center
    center_pixel = (naxis1 / 2, naxis2 / 2)

    # Convert pixel coordinates to RA and Dec
    center_ra_dec = wcs.pixel_to_world(center_pixel[0], center_pixel[1])

    # Calculate the size of the image in degrees
    # The size is determined by the diagonal distance from the center to the corner of the image
    corner_pixel = (0, 0)
    corner_ra_dec = wcs.pixel_to_world(corner_pixel[0], corner_pixel[1])
    size_deg = center_ra_dec.separation(corner_ra_dec).to(u.degree).value
    ra=center_ra_dec.ra.deg
    dec=center_ra_dec.dec.deg
    return ra,dec,size_deg


def do_xphot(filename,gaia_file,forced,nrows_max,outroot):

    print('XXX - do_xphot  %s gaia %s' % (filename,gaia_file))
    

    if forced:
        object_file=gaia_file
        phot_file=do_forced_photometry(filename,object_file,nrows_max,outroot)
    else:
        object_file=get_objects_from_image(filename,outroot)
        phot_file=do_photometry(filename,object_file,outroot)


    
    closest_objects_table = find_closest_objects(gaia_file, phot_file)
    if len(closest_objects_table)==0:
        print('Error: There are no objects that were xmatched')
        return
    
    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')

    do_fig(closest_objects_table,outroot)

    


def do_one(filename='LMC_c48_T08.r.t060.fits',gaia_cat_file='',forced=False,nrows_max=-1,outroot=''):
    '''
    Compare photometry in an image to photometry from Gaia

    '''

    



    try:
        x=fits.open(filename)
    except:
        print('Could not open %s' % filename)
        raise ValueError
        

    if gaia_cat_file!='' and os.path.isfile(gaia_cat_file)==True:
        gaia_file=gaia_cat_file
        print('Using existing GaiaCat file: %s' % gaia_cat_file)
    else:
        ra,dec,size_deg=get_size(filename)
        print('Making new GaiaCat file - %.2f %.2f %.2f' % (ra,dec,size_deg))
        print('do_one - RA, Dec, size: ',ra,dec,size_deg)
        gaia_file=GaiaCat.get_gaia(ra, dec, size_deg,outroot)


    do_xphot(filename,gaia_file,forced,nrows_max,outroot)

    return



def do_many(filenames=['LMC_c48_T08.r.t060.fits'],gaia_cat_file='',forced=True,nrows_max=10000,outroot=''):
    '''
    Compare photometry in an image to photometry from Gaia
    '''

    xra=[]
    xdec=[]
    xsize=[]


    for filename in filenames:
        try:
            x=fits.open(filename)
        except:
            print('do_many: Could not open %s' % filename)
            raise IOError
        ra,dec,size=get_size(filename)
        xra.append(ra)
        xdec.append(dec)
        xsize.append(size)

    xpos=Table([filenames,xra,xdec,xsize],names=['filename','RA','Dec','Size'])
    zpos, mapping =unique_rows_within_tol(xpos, tol=0.01)

    print("Finished getting positions ")

    gaia_files=[]
    for one in zpos:
        gaia_file=GaiaCat.get_gaia(one['RA'], one['Dec'], one['Size'],outroot='')
        gaia_files.append(gaia_file)
    zpos['gaia_file']=gaia_files

    print('Finished getting gaia tables for %d files' % len(zpos))

    xpos.write('xpos.txt',format='ascii.fixed_width_two_line',overwrite=True)
    zpos.write('zpos.txt',format='ascii.fixed_width_two_line',overwrite=True)


    xpos['gaia_file']=zpos['gaia_file'][mapping]
    xpos.write('xxpos.txt',format='ascii.fixed_width_two_line',overwrite=True)

    # At this point all of the gaia files that we need should exist

    for one in xpos:
        print('ZZZ',one)
        do_xphot(one['filename'],one['gaia_file'],forced,nrows_max,outroot)

    return


def do_dir(xdir='DECam_SWARP2/LMC_c37/T16',nrows_max=30000,forced=True):
    '''
    Process all of the images in a directory, and its 
    subdirecories
    '''

    xtab=ImageSum.table_create(xdir,outname=None)
    print('Starting %d files' %  len(xtab))

    do_many(xtab['filename'],gaia_cat_file='',forced=forced,nrows_max=nrows_max)

    return


def steer(argv):
    '''
    Run the script given choices from the command line

    Usage: PhotCompare.py -h -for -unf -dir -nmax -gcat file1
    '''

    global XDIR

    gaia_cat_file=''
    forced=True
    nrows_max=30000
    files=[]
    xdir=''
    outdir=''
    
    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(__doc__)
            return
        elif argv[i][:4]=='-for':
            forced=True
        elif argv[i][:4]=='-unf':
            forced=False
        elif argv[i]=='-dir':
            i+=1
            xdir=argv[i]
        elif argv[i]=='-out':
            i+=1
            out=argv[i]
        elif argv[i]=='-nmax':
            i+=1
            nrows_max=int(argv[i])
        elif argv[i]=='-gcat':
            i+=1
            gaia_cat_file=argv[i]
        elif argv[i][0]=='-':
            print('Unknown switch ',argv)
            return
        else:
            files.append(argv[i])
        i+=1


    if xdir!='':
        XDIR='_%s' % (xdir.replace('/','-'))
        do_dir(xdir=xdir,nrows_max=nrows_max,forced=forced)
        return



    i=1
    for one in files:
        print('\nProcessing %s (%d/%d' % (one,i,len(files)))
        do_one(one,gaia_cat_file,forced,nrows_max)
        i+=1


    return



           

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)   
    else:
        print (__doc__ )
