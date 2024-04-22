#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Read an astropy table containg RA's and Dec's
and return the field and tile located cose
to each object


Command line usage (if any):

    usage: MCLocate.py [-h] [-out whatever] filename

    where -h prints this documentation and quits
    and  -out whatever prints the output to a file named whatever.
    (if this is missing then the output file will be X+filename

Description:  

    The routine is hardwired to look at MC_tiles.txt in 
    the kred/config directory.  

    The input table needs to have fields with columns
    named 'RA', and 'Dec' with the positions of objects
    in degrees.

    The routine finds the field and tile closest to
    each position and adds that (and the separation to the output
    table, which will be printed to the scrren,
    and written to a file. 

Primary routines:

    doit

Notes:

    The routine does not check that that the object is actually
    in a tile; it just reports the distance in arcmin to the
    closest tile.
                                       
History:

240422 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack


def doit(filename='aida.txt',output=''):

    try:
        obj_tab=ascii.read(filename)
    except:
        print('Error: Could not read %s' % filename)
        return []
    
    tile_file='%s/config/MC_tiles.txt' % os.environ['KRED']
    print(tile_file)
    try:
        tile_tab=ascii.read(tile_file)
    except:
        print('Error: Could not read tile file in %s. Is $KRED defined?' % tile_file)
        return []

    # Convert RA and Dec columns to SkyCoord objects
    try:
        coords_table1 = SkyCoord(ra=obj_tab['RA'] * u.degree, dec=obj_tab['Dec'] * u.degree)
    except:
        print('Error: Something is wrong with the input table. Check for columns RA and Dec, with units of degrees')
        obj_tab.info()
        return []
    coords_table2 = SkyCoord(ra=tile_tab['RA'] * u.degree, dec=tile_tab['Dec'] * u.degree)

    # Find the closest objects in table2 to each object in table1
    closest_indices = []
    closest_separations = []

    for coord1 in coords_table1:
        separations = coord1.separation(coords_table2)
        min_index = np.argmin(separations)
        closest_indices.append(min_index)
        closest_separations.append(separations[min_index].to(u.arcmin).value)

    # Create a new table to store the closest objects and their separations
    obj_tab['Sep']=closest_separations*u.arcmin
    xtab=hstack([obj_tab,tile_tab[closest_indices]])

    return xtab


def steer(argv):
    '''
    Just a steering routine
    '''
    output=''
    filename=''

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        if argv[i][0:4]=='-out':
            i+=1
            output=argv[i]
        elif filename=='':
            filename=argv[i]
        else:
            print('Error: Badly formated inputs line: ',argv)
            return
        i+=1

    xtab=doit(filename,output)
    print(xtab)
    if output=='':
        output='X'+filename
    if len(xtab):
        xtab.write(output,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('Error: Nothing returned')




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        steer(sys.argv)
    else:
        print (__doc__)


