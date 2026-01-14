#!/usr/bin/env python
# coding: utf-8

"""Smash - SMASH DR2 Catalog Retrieval

Space Telescope Science Institute

Synopsis
--------

Retrieve stars from the SMASH (Survey of the MAgellanic Stellar History) DR2
catalog for a specified region on the sky.

Command Line Usage
------------------

::

    Usage: Smash.py [-h] [-r size_deg] [-ex] [-out root] ra dec

where:

    -h          Print this help message and exit
    -r size     Set the cone search radius in degrees (default: 0.5)
    -ex         Run example query (SMC region) without further inputs
    -out root   Set the output filename root (default: smash_cat)
    ra          Right Ascension of cone center in degrees
    dec         Declination of cone center in degrees

Description
-----------

This module retrieves photometric data from the SMASH DR2 catalog using
the NOAO Data Lab query client. It performs a cone search around a specified
position and returns filtered sources with high probability stellar
classifications.

The output table includes:

- Position (RA, Dec)
- Photometry in U, G, R, Z bands
- Probability of being a star (prob)

Sources are filtered to include only high-confidence stars (prob > 0.98)
brighter than r = 22 mag.

Primary Routines
----------------

smash_cone_search
    Perform a cone search on SMASH DR2 catalog.

steer
    Command line interface for catalog retrieval.

Notes
-----

Requires the NOAO Data Lab client library (dl).

History:

260108 ksl Coding begun

"""

import sys
from astropy.io import ascii, fits
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import os

from dl import queryClient as qc


def smash_cone_search(ra, dec, radius_deg, limit=10000000):
    """
    Perform a cone search on SMASH DR2 using Data Lab query client.

    Queries the SMASH DR2 object catalog for all sources within a
    specified angular radius of the given sky position.

    Parameters
    ----------
    ra : float
        Right Ascension of the cone center in degrees.
    dec : float
        Declination of the cone center in degrees.
    radius_deg : float
        Radius of the cone search in degrees.
    limit : int, optional
        Maximum number of rows to return (default: 10000000).

    Returns
    -------
    astropy.table.Table
        Table of SMASH DR2 objects within the search cone.
        Columns include position, photometry, and classification info.

    Examples
    --------
    >>> table = smash_cone_search(13.1867, -72.8286, 0.5)
    >>> print(f"Found {len(table)} sources")
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
    """
    Command line interface for SMASH catalog retrieval.

    Parses command line arguments and performs a cone search on the
    SMASH DR2 catalog, filtering results and writing to a FITS file.

    Parameters
    ----------
    argv : list
        Command line arguments (typically sys.argv).

    Notes
    -----
    Output is filtered to include only sources with:
    - prob > 0.98 (high confidence stellar classification)
    - rmag < 22 (bright enough for reliable photometry)

    Column names are standardized: ra->RA, dec->Dec, umag->U, etc.
    """

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

