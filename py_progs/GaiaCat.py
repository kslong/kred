#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Routines to handle retreiving information from the GAIA data base, and to interact
with datat that has been retrieved from the database.

This is NOT intended to be run from the command line at present but rather 
contains routines that should be called from other rotuines


Command line usage (if any):

    usage: None                                                                              ..

Description:  


                                       
History:

240318 ksl Coding begun
240527 ksl Speed up the catalog matching.
251105 ksl Split finding sources in an image from doing photometry
            on the sources
251130 ksl Cleaned up so this is just a routine for interacting with
    the GAIA catolog

'''

import os
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table,join,hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import timeit
import time
import pathlib
import os.path as path
import requests
from gaiaxpy import calibrate
from http.client import IncompleteRead


def load_Gaia(probe_service: bool = True):
    """
    Return the astroquery.gaia.Gaia class.
    Raises RuntimeError with a message distinguishing:
      - astroquery not installed
      - external Gaia assets/services unavailable

    Parameters
    ----------
    probe_service : bool, default True
        If True, perform a tiny network-dependent check to detect
        external service unavailability immediately. If False, only
        import Gaia and let service failures occur later at first use.

    Returns
    -------
    Gaia : type
        The Gaia class from astroquery.gaia.

    Raises
    ------
    RuntimeError
        With a distinguishing message for the two failure cases.
    """
    # 1) Import-time: distinguish "not installed"
    try:
        from astroquery.gaia import Gaia
    except ImportError as e:
        raise RuntimeError(
            "astroquery (or 'astroquery.gaia') is not installed in this environment."
        ) from e

    # Optional: reduce noise from astropy later (not required)
    try:
        from astropy.logger import log
        log.setLevel('ERROR')
        from astropy.utils import iers
        iers.conf.auto_download = False
        iers.conf.auto_max_age = None
    except Exception:
        pass  # if astropy isn't present or settings change, just continue

    # 2) Service probe: distinguish "assets unavailable"
    if probe_service:
        try:
            # Minimal, fast probe. Uncomment one of these depending on your preference:
            # A) Light sync query (hits TAP, but tiny):
            Gaia.launch_job("SELECT 1", dump_to_file=False)
            # B) Alternatively, accessing Gaia.tap can initialize the service object:
            # _ = Gaia.tap
        except Exception as e:
            raise RuntimeError(
                "Gaia external assets/services appear unavailable or unreachable."
            ) from e

    return Gaia

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



def get_gaia_spec(gaiaID, GAIA_CACHE_DIR='./GaiaSpec'):
    """
    Load or download and load from cache the spectrum of a gaia star, converted to erg/s/cm^2/A

    Note:

    This was  'appropiated'  from the lvmdrp

    """
    # create cache dir if it does not exist
    pathlib.Path(GAIA_CACHE_DIR).mkdir(parents=True, exist_ok=True)

    if path.exists(GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + ".csv") is True:
        print('Star is in cache')
        # read the tables from our cache
        gaiaflux = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + ".csv", format="csv"
        )
        gaiawave = Table.read(
            GAIA_CACHE_DIR + "/gaia_spec_" + str(gaiaID) + "_sampling.csv", format="csv"
        )
    else:
        print('Star is must be retrieved')
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
        1e4 * np.fromstring(gaiaflux["flux"][0][1:-1], sep=",")
    )  # W/s/nm -> in erg/s/cm^2/A

    results=Table([wave,flux],names=['WAVE','FLUX'])
    return results     


def get_gaia_mag28_flux(xid=4658615927801509760, gmag=15, wavelength=6563,dlambda=160):
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
    flux28=flux*10**(-0.4*(28-gmag))*dlambda
    return flux28

def get_gaia_mag28_ave(xid=4658615927801509760, gmag=15, wavelength=6563,dlambda=160):
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
    flux28=flux*10**(-0.4*(28-gmag))*dlambda
    return flux28
                         
                   

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


def get_gaia_from_archive(ra=84.92500000000001, dec=-66.27416666666667, rad_deg=0.3,
             outroot='', nmax=-1, redo=False, max_retries=3, retry_delay=5):
    '''
    Get data from the Gaia photometric catalog with retry logic for network errors.
    by conducting a cone search



    THIS IS UNTESTED, AND IS ONLY NEEDED IF WE NEED MORE DATA FROM ESA
    Even if it works, some changes are needed, including changing the
    output format to fits.

    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    rad_deg : float
        Search radius in degrees
    outroot : str
        Output root path (defaults to 'RA_Dec' format)
    nmax : int
        Maximum number of rows (-1 for no limit)
    redo : bool
        Whether to redo the query even if file exists
    max_retries : int
        Maximum number of retry attempts for network errors (default: 3)
    retry_delay : float
        Delay in seconds between retries (default: 5)

    Returns
    -------
    outfile : str or list
        Path to output file, or empty list if no objects retrieved
    '''

    try:
        Gaia = load_Gaia(probe_service=True)  # or False if you prefer no network on import
    except RuntimeError as err:
        # The message tells you which issue it is:
        # - "astroquery ... is not installed ..."
        # - "Gaia external assets/services appear unavailable ..."
        print(f"Dependency/service error: {err}")
        raise  RuntimeError

    if outroot == '':
        outroot = '%05.1f_%05.1f' % (ra, dec)

    os.makedirs('Gaia', exist_ok=True)
    outfile = 'Gaia/Gaia.%s.txt' % outroot

    if redo == False and os.path.isfile(outfile) == True:
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile

    print('get_gaia: Getting data for RA Dec of  %.5f %.5f and size of %.2f' % (ra, dec, rad_deg))

    Gaia.ROW_LIMIT = nmax  # Ensure the default row limit.
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

    # Retry loop for handling IncompleteRead errors
    r = None
    for attempt in range(max_retries):
        try:
            if attempt > 0:
                print(f'get_gaia: Retry attempt {attempt + 1}/{max_retries}...')

            j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))
            r = j.get_results()

            # If we get here, the query succeeded
            break

        except IncompleteRead as e:
            print(f'get_gaia: IncompleteRead error on attempt {attempt + 1}: {e}')
            if attempt < max_retries - 1:
                print(f'get_gaia: Retrying in {retry_delay} seconds...')
                time.sleep(retry_delay)
            else:
                print('get_gaia: Max retries reached. Query failed.')
                raise

        except Exception as e:
            print(f'get_gaia: Unexpected error on attempt {attempt + 1}: {type(e).__name__}: {e}')
            if attempt < max_retries - 1:
                print(f'get_gaia: Retrying in {retry_delay} seconds...')
                time.sleep(retry_delay)
            else:
                print('get_gaia: Max retries reached. Query failed.')
                raise

    if r is None or len(r) == 0:
        print('Error: get_gaia: No objects were retrieved')
        return []

    # Process and rename columns
    r.rename_column('ra', 'RA')
    r.rename_column('dec', 'Dec')
    try:
        r.rename_column('source_id', 'Source_name')
    except:
        r.rename_column('SOURCE_ID', 'Source_name')
    r.rename_column('phot_g_mean_mag', 'G')
    r.rename_column('phot_bp_mean_mag', 'B')
    r.rename_column('phot_rp_mean_mag', 'R')
    r.rename_column('teff_gspphot', 'teff')
    r.rename_column('logg_gspphot', 'log_g')
    r.rename_column('distance_gspphot', 'D')

    r['Source_name', 'RA', 'Dec', 'B', 'G', 'R', 'teff', 'log_g', 'D'].write(
        outfile, format='ascii.fixed_width_two_line', overwrite=True
    )

    print('Wrote %s with %d objects' % (outfile, len(r)))
    return outfile

def get_gaia_from_archive_old(ra=84.92500000000001, dec= -66.27416666666667, rad_deg=0.3,outroot='',nmax=-1,redo=False):
    '''
    Get data from the Gaia photometric catalog, by exucuting a cone search on the Gaia archive

    Notes

    The routine returns a Runtime Error if the archive is not availalbe (or if astroquery is not installed)

    '''

    try:
        Gaia = load_Gaia(probe_service=True)  # or False if you prefer no network on import
    except RuntimeError as err:
        # The message tells you which issue it is:
        # - "astroquery ... is not installed ..."
        # - "Gaia external assets/services appear unavailable ..."
        print(f"Dependency/service error: {err}")
        raise  RuntimeError

    if outroot=='':
        outroot='%06.2f_%06.2f' % (ra,dec)
    
    os.makedirs('Gaia',exist_ok=True)
    outfile='Gaia/Gaia.%s.txt' % outroot

    if redo==False and os.path.isfile(outfile)==True:
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile
    
    print('get_gaia: Getting data for RA Dec of  %.5f %.5f and size of %.2f' % (ra,dec,rad_deg))

    

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
    try:
        r.rename_column('source_id','Source_name')
    except:
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


def get_gaia(ra=84.92500000000001, dec= -66.27416666666667, size_deg=0.3,outroot='',filename='Gaia_MagClouds.fits'):
    '''
    Retrieve entries from a local table containg information about stars that are in the Gaia catolog.

    Description

    This routine retrieves information a local table that must be installed either in the directroy
    from which the program is being run, or from a specific directory, namely kred/xdata.  



    Notes:

    The file that covers the SMC and LMC can be found on box
    
    Note that this routine does not query the GAIA archive unlike many of the routines in this module

    Unlike some other routines the file that is covers a  'square' regionin RA and Dec.


    '''

    # first locate the file
    
    if os.path.isfile(filename):
        xfilename=filename
    elif os.path.isfile('%s/%s' % ('Gaia',filename)):
        xfilename='%s/%s' % ('Gaia',filename)
    else:
        KRED = os.environ.get("KRED")
        if KRED is not None:
            if os.path.isfile('%s/%s/%s' % (KRED,'xdata',filename)):
                xfilename='%s/%s/%s' % (KRED,'xdata',filename)
            else:
                raise IOError('Could not locate %s' % filename)
        else:
              raise IOError('Enviroment variable KRED is not set')

    if xfilename.count('fits'):
        xtab=Table.read(xfilename)
    else:
        xtab=ascii.read(xfilename)

    dec_min=dec-0.5*size_deg
    dec_max=dec+0.5*size_deg
    xscale=np.cos(dec/57.29578)
    factor=0.5*size_deg/xscale
    ra_min=ra-factor
    ra_max=ra+factor

    # print('XX get_gaia', ra_min,ra_max,dec_min,dec_max)

    mask=((ra_min < xtab['RA']) & (xtab['RA']< ra_max) &  (dec_min < xtab['Dec']) &  (xtab['Dec'] < dec_max))

    ftab=xtab[mask]

    if outroot=='':
        outroot='%06.2f_%06.2f' % (ra,dec)
    
    os.makedirs('Gaia',exist_ok=True)
    outfile='Gaia/Gaia.%s.fits' % outroot

    ftab.write(outfile,format='fits',overwrite=True)

    return outfile
    





def steer(argv):
    '''
    Run the script given choices from the command line

    Usage: PhotCompare.py -h -for -unf -dir -nmax -gcat file1
    '''

    print('This is not a runtime routine (currently)')
    print(__doc__)
    #  i=1
    # while i<len(argv):
    #     if argv[i].count('-h'):
    #         print(__doc__)
    #         return
    #     elif argv[i][0]=='-':
    #         print('Unknown switch ',argv)
    #         return
    #     else:
    #         files.append(argv[i])
    #     i+=1


    return


           

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)   
    else:
        print (__doc__ )
