#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve Smash catolog stars


Command line usage (if any):

    usage:  Smash.py [-h] [-r size_deg] [-ex] [-out root]  ra dec

    where

    -h prints this documentation and exits
    =r size changes the radius of the cone search to a size in deg
    -ex runs and example without further inputs
    -out whatever sets the root name of the output table
    ra, and dec in degrees

Description:  

    Retrieve smash catalog stars from a circular region on the sky

Primary routines:

    doit

Notes:
                                       
History:

260108 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
from astropy.table import Table,join
import numpy as np
import matplotlib.pyplot as plt
import os

from dl import queryClient as qc

def smash_cone_search(ra, dec, radius_deg, limit=10000000):
    """
    Perform a cone search on SMASH DR2 using dl.queryClient.

    Parameters
    ----------
    ra : float
        Right Ascension of the cone center in degrees.
    dec : float
        Declination of the cone center in degrees.
    radius_deg : float
        Radius of the cone in degrees.
    limit : int
        Maximum number of rows to return.

    Returns
    -------
    Astropy Table
        Table of objects within the cone.
    """
    sql = f"""
    SELECT *
    FROM smash_dr2.object
    WHERE q3c_radial_query(
        ra, dec,
        {ra}, {dec}, {radius_deg}
    )
    LIMIT {limit}
    """
    # Query returns an Astropy Table directly
    return qc.query(sql=sql, fmt="table")

def steer(argv):
    '''
    usage:  Smash.py -h -r size -ex -out root  ra dec
    '''

    ra=-1.
    dec=-1.
    radius     = 0.5  # degrees
    outroot='smash_cat'

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:2]=='-r':
            i+=1
            radius = eval(argv[i])
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        elif argv[i][:3]=='-ex':
            ra  = 13.1867
            dec = -72.8286
            break
        elif argv[i][0]=='-' and  ra==-1:
            print('Error: Unknow option :',argv)
            return
        elif ra==-1.: 
            ra=eval(argv[i])
        elif dec==-1.:
            dec=eval(argv[i])
        else:
            print('Error: Cannot parse command line: ',argv)
        i+=1

    # Example: Small Magellanic Cloud region

    print("Querying SMASH DR2...")
    table = smash_cone_search(ra, dec, radius)
    print(f"Retrieved {len(table)} rows")
    table = table[(table['prob'] > 0.98) & (table['rmag'] < 22)]
    print(f"After filtering (prob > 0.98, rmag < 22): {len(table)} rows")
    table.pprint(max_lines=10)
    table.rename_column('ra','RA')
    table.rename_column('dec','Dec')
    table.rename_column('umag','U')
    table.rename_column('gmag','G')
    table.rename_column('rmag','R')
    table.rename_column('zmag','Z')
    if outroot=='':
        outroot='Smash_%06.2f_%0.6.2f_%5.2' % (ra,dec,radius)
    table.write(outroot+'.fits',format='fits',overwrite=True)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)

