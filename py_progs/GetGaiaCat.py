#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve the Gaia3 catalog for tiles in the LMC and SMC fields


Command line usage (if any):

    usage: GetGaiaCat.py [-h] [-redo] [-all] Field1 ... Tile1  ..

    where

    -h prints this brief header info and quits
    -redo forces a new retrieval of information from the archive if the data has been retrieve previously
    -all calls data for all tiles to be retrievd

    and


    Field is the name of one or more fields in the LMC or SMC
    Rile is the name of one or more tiles


Description:  

    The routine retrieves catalog data from the GAIA archive and stores the data in a local directory GaiaCat
    The resulting data is labelled by the field and tile_no.

    The results are stored in a local GaiaCat directory.  If the data has already been retrieved then there
    will be not attempt to retrive the data again unless the -redo flag is included from the command line

Primary routines:

    do_one_tile - processes data for one tile in one field

Notes:

    The results of this routine are not part of of the kred git repostiroy because the files produced are 
    fairly large. However, we suggest that if one is likely to have several directory stuctures for running
    kred, say a production and a test directory, then one puts a symbolic link to a master GaiaCat 
    diretory to prevent duplication.  The routine takes about an hour to run on a field.
                                       
History:

240504 ksl Coding begun

'''
# # Get the GAIA3 catalog for one or more tiles 



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




def get_gaia_all(ra=84.92500000000001, dec= -66.27416666666667, rad_deg=0.3,outroot='',nmax=-1,redo=False):
    '''
    Get data from the Gaia photometric catalog
    '''
    if outroot=='':
        outroot='%.1f_%.1f' % (ra,dec)

    outfile='./GaiaCat/Gaia.%s.txt' % outroot
    
    os.makedirs('./GaiaCat',exist_ok=True)
    
    if redo==False and os.path.isfile(outfile)==True:
        print('Outfile %s exists, to redo set redo to True' % (outfile))
        return
        

    if redo==False and os.path.isfile(outfile)==True:
        print('get_gaia: %s exists so returning, use redo==True to redo' % outfile)
        return outfile

    print('get_gaia_all: Getting data from RA Dec of  %.5f %.5f and size of %.2f to be stored in %s ' % (ra,dec,rad_deg,outfile))



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





def do_one_tile(field='LMC_c45',tile='T08',redo=False):
    
    tile_file='%s/config/MC_tiles.txt' % os.environ['KRED']
    print(tile_file)
    try:
        tile_tab=ascii.read(tile_file)
    except:
        print('Error: Could not read tile file in %s. Is $KRED defined?' % tile_file)
        return False
    
    # tile_tab.info()
    # print(tile_tab)
    
    xtile=tile_tab[tile_tab['Field']==field]
    # print(xtile)
    if tile!='':
        xtile=xtile[xtile['Tile']==tile]
    # print(xtile)
    if len(xtile)>0:
        for one in xtile:
            get_gaia_all(one['RA'],one['Dec'],1.414/2.*one['Size'],outroot='%s_%s' % (one['Field'],one['Tile']),redo=redo)
        return True
    else:
        print('Error: do_one_tile. Nothing in %s for %s %s' % (tile_file,field,tile))
        return False
                   
    

def steer(argv):
    '''
    This is just a steering routine so that the program can be run conveniently from the command line
    '''
    i=1
    redo=False
    xall=False
    tiles=[]
    fields=[]
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-redo':
            redo=True
        elif argv[i]=='-all':
            xall=True
        elif argv[i][0]=='-':
            print('Error: Problem with command line: ',argv)
            return
        elif argv[i][0]=='T':
            tiles.append(argv[i])
        else:
            fields.append(argv[i])
        i+=1

    for one_field in fields:
        if xall==True:
            do_one_tile(one_field,'',redo)
        else:
            for one_tile in tiles:
                do_one_tile(one_field,one_tile,redo)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__ )
