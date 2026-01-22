#!/usr/bin/env python
"""reg2master - Convert DS9 Region File to Master Table

Space Telescope Science Institute

Synopsis
--------

Convert a DS9 region file into a master table file with standardized
columns for source positions and region geometry.

Command Line Usage
------------------

::

    reg2master.py [-h] regionfile [masterfile]

**Required Arguments:**

regionfile
    Input DS9 region file containing source region definitions.

**Optional Arguments:**

masterfile
    Output master table filename. If not specified, defaults to
    ``<regionfile>.txt``.

-h
    Display this help message and exit.

Description
-----------

This program parses a DS9 region file and extracts source positions and
region geometry into a standardized master table format. The output can
be used as input for photometry routines or converted back to region
files using master2reg.py.

Supported region types:

* **circle** - Circular apertures
* **ellipse** - Elliptical apertures
* **box** - Rectangular regions
* **annulus** - Circular annuli

Output
------

An ASCII table in fixed_width_two_line format with columns:

* Source_name : str - Source identifier (from region label or auto-generated)
* RA : float - Right Ascension in degrees
* Dec : float - Declination in degrees
* RegType : str - Region type
* Major : float - Major axis/radius in arcseconds
* Minor : float - Minor axis in arcseconds
* Theta : float - Position angle in degrees
* Color : str - Region color from DS9 file

Notes
-----

* Source names are extracted from the ``text={}`` field in region definitions
* If no name is provided, names are auto-generated as 'zzz001', 'zzz002', etc.
* Colors are preserved from the region file for downstream use
* For boxes, Theta corresponds to the angle E of N for the Minor axis
* The same convention applies to ellipses

History
-------

090109 ksl
    Initial coding

111111 ksl
    Fixed for pydocs compatibility

190507 ksl
    Added Color column extraction

241230 ksl
    Added box region type support

.. moduleauthor:: KSL
"""

import sys
import os
from astropy.table import Table
import numpy


def radec2deg(ra, dec):
    """Convert RA/Dec strings to decimal degrees.

    Parameters
    ----------
    ra : str
        Right Ascension as 'HH:MM:SS' string or decimal degrees string.
    dec : str
        Declination as 'DD:MM:SS' string or decimal degrees string.

    Returns
    -------
    tuple of (float, float)
        RA and Dec in decimal degrees.

    Notes
    -----
    RA in sexagesimal format is assumed to be in hours and is
    converted to degrees (multiplied by 15).
    """

    r=ra.split(':')
    d=dec.split(':')

    rr=float(r[0])
    if len(r)>1:
        rr=rr+float(r[1])/60.
    if len(r)>2:
        rr=rr+float(r[2])/3600.
    if len(r)>1:
        rr=15.*rr  # Since we assume ra was in hms
    
    dd=float(d[0])
    x=0
    if len(d)>1:
        x=x+float(d[1])/60.
    if len(d)>2:
        x=x+float(d[2])/3600.
    
    if dd > 0:
        dd=dd+x
    else:
        dd=dd-x
    
    return rr,dd  


def size2arcsec(word):
    """Convert a size string to arcseconds.

    Parameters
    ----------
    word : str
        Size value with optional unit suffix ('' for arcsec, ' for arcmin).

    Returns
    -------
    float
        Size in arcseconds.

    Notes
    -----
    Recognizes:

    * ``"`` suffix - value is in arcseconds
    * ``'`` suffix - value is in arcminutes (converted to arcsec)
    * No suffix - value is returned as-is (assumed arcseconds)
    """
    if word.count('"')==1: # We have arcsec
        value=float(word.rstrip('"'))
    elif word.count("'")==1: # We have arcsec
        value=60.*float(word.rstrip("'"))
    else:
        print('Error : Could not parse ',word,' to arcsec')
        value=float(word)
    return value

    


def read_regions(filename):
    """Read and parse a DS9 region file.

    Parameters
    ----------
    filename : str
        Path to the DS9 region file.

    Returns
    -------
    tuple of (str, list)
        Coordinate type ('fk5', 'image', 'physical', or 'unknown') and
        list of region records. Each record is a list:
        [name, ra, dec, regtype, major, minor, theta, color].

    Notes
    -----
    Parses the following region types: circle, ellipse, box, annulus.
    Source names are extracted from ``text={}`` fields or auto-generated
    as 'zzz001', 'zzz002', etc. if not present.
    """
    f=open(filename,'r')
    type='unknown'
    xcolor='unknown'

    records=[]
    source_no=1

    lines=f.readlines()
    for line in lines:
        line=line.replace('(',' ')
        line=line.replace(')',' ')
        line=line.replace('{',' ')
        line=line.replace('}',' ')
        line=line.replace('=',' ')
        line=line.replace(',',' ')
        line=line.split()
        if len(line)==0:
            continue

        # process a generic line
        # print(line)
 

        if line[0]=='global':
            # print('Global:',line)
            j=0
            while j<len(line)-1:
                if line[j]=='color':
                    xcolor=line[j+1]
                j+=1

        # find the name
        j=0
        name='unknown'
        color='unknown'
        while j<len(line)-1:
            if line[j]=='text':
                name=line[j+1]
            if line[j]=='color':
                color=line[j+1]
            j=j+1
        if name=='unknown':
            name='zzz%03d' % source_no
        if color=='unknown':
            color=xcolor



        # Now process circles
        if len(line)==0:
            pass
        elif line[0]=='fk5':
            type='fk5'
        elif line[0]=='image':
            type='image'
        elif line[0]=='physical':
            type='physical'
        elif line[0]=='circle':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            record=[name,rr,dd,'circle',x1,0.0,0.0,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='ellipse':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=float(line[5])
            record=[name,rr,dd,'ellipse',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='box':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=float(line[5])
            record=[name,rr,dd,'box',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='annulus':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=0.0
            record=[name,rr,dd,'annulus',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        else:
            print('Did not use :',line)
        # print(line)
    f.close()
    return type,records




def write_masterfile(masterfile, records, source='unknown', type='unknown'):
    """Write a master table file from parsed region records.

    Parameters
    ----------
    masterfile : str
        Output filename for the master table.
    records : list
        List of region records from read_regions().
    source : str, optional
        Source file identifier. Default is 'unknown'.
    type : str, optional
        Coordinate type from region file. Default is 'unknown'.

    Returns
    -------
    None
        Writes table to disk in ascii.fixed_width_two_line format.

    Notes
    -----
    Output columns: Source_name, RA, Dec, RegType, Major, Minor, Theta, Color.
    Numeric columns are formatted with appropriate precision.
    """

    names=['Source_name','RA','Dec','RegType','Major','Minor','Theta','Color']

    x=Table()

    i=0
    while i<len(names):
        value=[]
        for record in records:
            # print(record)
            value.append(record[i])
        x[names[i]]=value
        i+=1

    x['RA'].format='10.6f'
    x['Dec'].format='10.6f'
    x['Major'].format='8.2f'
    x['Minor'].format='8.2f'
    x['Theta'].format='8.2f'
    x['Color'].format='10s'

    x.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)

    return

def doit(regionfile, masterfile=''):
    """Create a master table from a region file.

    Parameters
    ----------
    regionfile : str
        Input DS9 region file.
    masterfile : str, optional
        Output master table filename. If empty, defaults to
        ``<regionfile>.txt``.

    Returns
    -------
    None
        Writes master table to disk.
    """

    if masterfile=='':
        masterfile=regionfile+'.txt'
    elif masterfile[-4:].count('.txt')==0:
        masterfile=masterfile+'.txt'
    

    print('Regionfile          :', regionfile)
    print('Masterfile          :', masterfile)

    # Finished getting input information

    # read the regionfile

    type,records=read_regions(regionfile)

    print('No of records parsed ', len(records))

    write_masterfile(masterfile,records,regionfile,type)


def steer(argv):
    """Parse command line arguments and execute conversion.

    Parameters
    ----------
    argv : list of str
        Command line arguments (typically sys.argv).

    Returns
    -------
    None
        Calls doit() to perform the conversion.
    """
    regionfile = ''
    root = ''
    i = 1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='':
            print('Could not parse command line: ',argv)
            return
        elif regionfile=='':
            regionfile=argv[i]
        elif root=='':
            root=argv[i]
        else:
            print('Too many arguments in command line: ',argv[i])
            return

        i+=1

    doit(regionfile,masterfile=root)






# Subroutines and function calls are above.  The driving routine is below
if __name__ == "__main__":
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)


