#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Compare the brightneess of objects found in kred-prociessed images to Gaia
assuming our standard rescaling


Command line usage (if any):

    usage: PhotCompare.py file1 file2 ...

Description:  

    The routines processes one or more files comparing Gaia photometry
    and images produced with kred

Primary routines:

    doit

Notes:
                                       
History:

240318 ksl Coding begun

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


from astroquery.gaia import Gaia

import pathlib
import os.path as path
import requests
from gaiaxpy import calibrate




def get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./gaia'):
    """
    Load or download and load from cache the spectrum of a gaia star, converted to erg/s/cm^2/A

    Note that I have 'appropiated' this from the lvm drp
    """
    # create cache dir if it does not exist
    pathlib.Path(GAIA_CACHE_DIR).mkdir(parents=True, exist_ok=True)

    if path.exists(GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + ".csv") is True:
        # read the tables from our cache
        gaiaflux = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + ".csv", format="csv"
        )
        gaiawave = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + "_sampling.csv", format="csv"
        )
    else:
        # need to download from Gaia archive
        CSV_URL = (
            "https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=XP_CONTINUOUS&ID=Gaia+DR3+"
            + str(gaiaID)
            + "&format=CSV&DATA_STRUCTURE=RAW"
        )
        FILE = GAIA_CACHE_DIR + "/XP_" + str(gaiaID) + "_RAW.csv"

        with requests.get(CSV_URL, stream=True) as r:
            r.raise_for_status()
            if len(r.content) < 2:
                return []
            with open(FILE, "w") as f:
                f.write(r.content.decode("utf-8"))

        # convert coefficients to sampled spectrum
        _, _ = calibrate(
            FILE,
            output_path=GAIA_CACHE_DIR,
            output_file="gaia_spec_" + str(gaiaID),
            output_format="csv",
        )
        # read the flux and wavelength tables
        gaiaflux = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + ".csv", format="csv"
        )
        gaiawave = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + "_sampling.csv", format="csv"
        )

    # make numpy arrays from whatever weird objects the Gaia stuff creates
    wave = np.fromstring(gaiawave["pos"][0][1:-1], sep=",") * 10  # in Angstrom
    flux = (
        1e7 * 1e-1 * 1e-4 * np.fromstring(gaiaflux["flux"][0][1:-1], sep=",")
    )  # W/s/micron -> in erg/s/cm^2/A

    results=Table([wave,flux],names=['WAVE','FLUX'])
    return results     


def get_gaia_mag27_flux(xid=4658615927801509760, gmag=15, wavelength=6563,dlambda=160):
    '''
     Get the giaa flux of a star at a particular wavelength and calculate thoe
     total flux assuming in the bandpass from this if it were the same star
     at mag27`.  

     Note:
     Added 240503. This is one possible way of calculating 
    '''
    xtab=get_gaia_spec(xid)
    if len(xtab)==0:
        print('Error: Could not get gaia spectrum for gaia ID %s' % (xid))
        return None
    # xtab.info()
    # print(xtab)
    i=0
    while xtab['WAVE'][i] < wavelength and i<len(xtab):
        i+=1
    #print(xtab['WAVE'][i])
    frac=(wavelength-xtab['WAVE'][i-1])/(xtab['WAVE'][i]-xtab['WAVE'][i-1])
    # print(frac)
    flux=(1-frac) * xtab['FLUX'][i-1]+frac*xtab['FLUX'][i]
    # print(flux)
    flux27=flux*10**(-0.4*(27-gmag))*dlambda
    return flux27

def get_gaia_mag27_ave(xid=4658615927801509760, gmag=15, wavelength=6563,dlambda=160):
    '''
     Get the averge gaia flux of a star in a partcular wavelength band

     Note:
     Added 240503 - This is a variant of the routine above
    '''
    xtab=get_gaia_spec(xid)
    if len(xtab)==0:
        print('Error: Could not get gaia spectrum for gaia ID %s' % (xid))
        return None
    # xtab.info()
    # print(xtab)
    wmax=wavelength+dlambda/2.
    wmin=wavelength-dlambda/2
    z=xtab[xtab['WAVE'] < wmax]
    z=z[z['WAVE']>wmin]
    flux=np.average(z['FLUX'])
    flux27=flux*10**(-0.4*(27-gmag))*dlambda
    return flux27
                         
                   
                   


def get_gaia_flux(xid=4658604348568208768):
    '''
    Get the flux for a Gaia star as observed through the various filters

    240503 - I am not convince this is useful for anything, as this is
    not normalized an way, and it is unrelated to the flux calibraiton
    memo.
    '''
    xtab=get_gaia_spec(xid)
    if len(xtab)==0:
        print('Error: Could not get gaia spectrum for gaia ID %s' % (xid))
        return

    test_dir=os.path.dirname(__file__)
    print('test2',test_dir)


    data_dir=os.path.dirname(__file__).replace('py_progs','data')

    print('test',data_dir)

    xfilt=ascii.read('%s/%s' % (data_dir,'n662.txt'))
    xtab['HA_TRANS']= np.interp(xtab['WAVE'], xfilt['WAVE'], xfilt['TRANS'],
                                       left=0, right=0)


    xfilt=ascii.read('%s/%s' % (data_dir,'n673.txt'))
    xtab['S2_TRANS']= np.interp(xtab['WAVE'], xfilt['WAVE'], xfilt['TRANS'],
                                       left=0, right=0)
    xfilt=ascii.read('%s/%s' % (data_dir,'r.txt'))
    xtab['R_TRANS']= np.interp(xtab['WAVE'], xfilt['WAVE'], xfilt['TRANS'],
                                       left=0, right=0)


    xfilt=ascii.read('%s/%s' % (data_dir,'n708.txt'))
    xtab['N708_TRANS']= np.interp(xtab['WAVE'], xfilt['WAVE'], xfilt['TRANS'],
                                       left=0, right=0)

    xtab.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)

    dw=20.

    r_flux=np.dot(xtab['FLUX'],xtab['R_TRANS'])*dw
    ha_flux=np.dot(xtab['FLUX'],xtab['HA_TRANS'])*dw
    s2_flux=np.dot(xtab['FLUX'],xtab['S2_TRANS'])*dw
    n708_flux=np.dot(xtab['FLUX'],xtab['N708_TRANS'])*dw

    return ha_flux,s2_flux,r_flux,n708_flux



def get_gaia(ra=84.92500000000001, dec= -66.27416666666667, rad_deg=0.3,outroot='',nmax=-1,redo=False):
    '''
    Get data from the Gaia photometric catalog
    '''
    if outroot=='':
        outroot='%.1f_%.1f' % (ra,dec)
    
    outfile='Gaia.%s.txt' % outroot

    if redo==False and os.path.isfile(outfile)==True:
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile
    
    print('get_gaia: Getting data from RA Dec of  %.5f %.5f and size of %.2f' % (ra,dec,rad_deg))

    



    Gaia.ROW_LIMIT = nmax  # Ensure the default row limit.

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

    j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))


    r = j.get_results()
    
    if len(r)==0:
        print('Error: get_gaia: No objects were retrieved')
        return []

    # print(r.info())

    r.rename_column('ra','RA')
    r.rename_column('dec','Dec')
    r.rename_column('SOURCE_ID','Source_name')
    r.rename_column('phot_g_mean_mag','G')
    r.rename_column('phot_bp_mean_mag','B')
    r.rename_column('phot_rp_mean_mag','R')
    r.rename_column('teff_gspphot','teff')
    r.rename_column('logg_gspphot','log_g')
    r.rename_column('distance_gspphot','D')
    
    r['Source_name','RA','Dec','B','G','R','teff','log_g','D'].write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d objects' %(outfile,len(r)))
    return outfile





def get_photometry(filename='LMC_c48_T08.r.t060.fits',outroot=''):
    
    try:
        x=fits.open(filename)
    except:
        print('Error: get_photometry: could not open %s' % filename)
        return 'Error'

    print('get_photometry: Beginning photometry of %s' % filename)

    xexptime=x['PRIMARY'].header['EXPTIME']
    xfilter=x['PRIMARY'].header['FILTER']

    
    if outroot=='':
        words=filename.split('/')
        outroot=words[-1].replace('.fits','')
        
    image_wcs=WCS(x[0].header)
    image=x[0].data
    # print(np.median(image))

    image-=np.median(image)
    
    bkg_sigma = mad_std(image)  

    daofind = DAOStarFinder(fwhm=4.0, threshold=3.0 * bkg_sigma)  

    sources = daofind(image)  

    for col in sources.colnames:  

        sources[col].info.format = '%.8g'  # for consistent table output

    # print(sources) 
    sources.write('%s_sources.txt' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  

    apertures = CircularAperture(positions, r=4.0)  
    annulus_apertures=CircularAnnulus(positions,r_in=4, r_out=8)

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


    
    # print(phot_table)  
    outfile='%s_phot.txt' % outroot
    phot_table.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d objects' % (outfile,len(phot_table)))
    return outfile


def find_closest_objects(table1_path, table2_path, max_sep=0.5):
    '''
    Find the objects with a given distance given two astropy tables.  The
    routine returns only the closest object that satisvies this criterion.
    '''
    # Read the two Astropy tables
    try:
        table1 = Table.read(table1_path,format='ascii.fixed_width_two_line')
    except:
        print('Error: find_closest_objects: could not read %s' % table1_path)
        return []
    
    try:
        table2 = Table.read(table2_path,format='ascii.fixed_width_two_line')
    except:
        print('Error: find_closest_objects: could not read %s' % table2_path)
        return []

    print('get_closest_objects: Beginning x-match of %s and %s' % (table1_path,table2_path))
    
    # Convert RA and Dec columns to SkyCoord objects
    coords_table1 = SkyCoord(ra=table1['RA'] * u.degree, dec=table1['Dec'] * u.degree)
    coords_table2 = SkyCoord(ra=table2['RA'] * u.degree, dec=table2['Dec'] * u.degree)
    
    # Find the closest objects in table2 to each object in table1
    closest_indices = []
    closest_separations = []

    for coord1 in coords_table1:
        separations = coord1.separation(coords_table2)
        min_index = np.argmin(separations)
        closest_indices.append(min_index)
        closest_separations.append(separations[min_index].to(u.arcsec).value)

    # Create a new table to store the closest objects and their separations
    table1['Sep']=closest_separations*u.arcsec
    del table2['Source_name']
    del table2['RA']
    del table2['Dec']
    table1['Sep'].format='.3f'
    table1['RA'].format='.6f'
    table1['Dec'].foramt='6f'
    xtab=hstack([table1,table2[closest_indices]])
    
    xtab=xtab[xtab['Sep']<max_sep]
    
    print('Of %d objects in %s and %d objects in %s, found %d matches' % (len(table1),table1_path,len(table2),table2_path,len(xtab)))
    
    if len(xtab):
        words=table1_path.split('/')
        one=words[-1].replace('.txt','')
        words=table2_path.split('/')
        two=words[-1].replace('.txt','')
        outfile='xmatch_%s_%s.txt' % (one,two)
        xtab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('Error: There are no objects that are closer thn %f arcsec' % max_sep)
        return []

    return xtab



def do_all(filename='LMC_c48_T08.r.t060.fits',outroot=''):
    '''
    Compare photometry in an image to photometry from Gaia
    '''
    try:
        x=fits.open(filename)
    except:
        print('Could not open %s' % filename)
        return
    
    
    wcs = WCS(x[0].header)
         
    # Get the shape of the image
    naxis1 = x[0].header['NAXIS1']
    naxis2 = x[0].header['NAXIS2']
        
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
    
    # print(ra,dec,size_deg)
    
    gaia_file=get_gaia(ra, dec, size_deg,outroot,nmax=-1)
    
    phot_file=get_photometry(filename,outroot)
    
    closest_objects_table = find_closest_objects(gaia_file, phot_file)
    if len(closest_objects_table)==0:
        print('Errror: There are no objects that were xmatched')
        return
    
    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')

    do_fig(closest_objects_table,outroot)

    return


def do_fig(xtab,outroot):
    plt.figure(1,(12,6))
    plt.clf()
    plt.subplot(1,2,1)
    # plt.plot(xtab['G'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.plot(xtab['G'],xtab['phot_mag'],'.',alpha=.05)
    plt.xlabel('Gaia G mag')
    plt.ylabel('DECam mag')
    plt.plot([11,22],[11,22],'k-')
    plt.text(13,20,outroot)

    plt.ylim(11,22)
    plt.xlim(11,22) 



    plt.tight_layout()
    plt.subplot(1,2,2)
    # plt.plot(xtab['R'],27-2.5*np.log10(xtab['aperture_sum']),'.',alpha=.05)
    plt.plot(xtab['R'],xtab['phot_mag'],'.',alpha=.05)
    plt.xlabel('Gaia R mag')
    plt.ylabel('DECam mag')
    plt.plot([11,22],[11,22],'k-')
    plt.ylim(11,22)
    plt.xlim(11,22)  
    plt.tight_layout()
    plt.savefig('%s.png' % outroot)

    

def steer(argv):
    '''
    Run the script given choices from the command line
    '''

    files=[]
    
    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Unknown switch ',argv)
            return
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        print('Processing %s' % one)
        do_all(one)



           

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)   
    else:
        print (__doc__ )
