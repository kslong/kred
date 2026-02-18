#!/usr/bin/env python
"""master2reg - Convert Master Table to DS9 Region File

Space Telescope Science Institute

Synopsis
--------

Convert a master table file containing source positions and region definitions
into a DS9-compatible region file for visualization and analysis.

Command Line Usage
------------------

::

    master2reg.py [-h] [-r radius] [-color color] masterfile [regionfile]

**Required Arguments:**

masterfile
    Input master table file (ASCII or FITS format) containing source
    definitions with columns for position and region geometry.

**Optional Arguments:**

regionfile
    Output DS9 region file name. If not specified, defaults to
    ``<masterfile>.reg``.

-h
    Display this help message and exit.

-r radius
    Default radius in arcseconds for circular regions when size is not
    specified in the input file. Default is 3.0 arcsec.

-color color
    Color for region outlines (e.g., red, green, blue, cyan, magenta).
    Default is red.

Description
-----------

This program reads a master table and converts it to a DS9 region file.
The master table can be in various formats:

**Full format (with region geometry):**

The table should contain the following columns:

* Source_name : str - Source identifier
* RA : float - Right Ascension in degrees
* Dec : float - Declination in degrees
* RegType : str - Region type ('circle', 'ellipse', 'box', or 'annulus')
* Major : float - Major axis or radius in arcseconds
* Minor : float - Minor axis in arcseconds (0 for circles)
* Theta : float - Position angle in degrees

**Minimal format (positions only):**

If only Source_name, RA, and Dec are provided, circular regions with
the default radius (-r option) will be created.

**RA/Dec only format:**

If only RA and Dec columns are present, source names will be auto-generated.

Output
------

A DS9 region file in FK5 coordinate format with:

* Region shapes (circle, ellipse, box, annulus)
* Source labels
* Color specifications

Notes
-----

Supported region types:

* **circle** - Circular aperture (Major = radius)
* **ellipse** - Elliptical aperture (Major, Minor = semi-axes, Theta = angle)
* **box** - Rectangular region (Major, Minor = dimensions, Theta = angle)
* **annulus** - Circular annulus (Major = outer radius, Minor = inner radius)

History
-------

090109 ksl
    Initial coding

111216 ksl
    Modified to produce region file from first 3 columns if geometry
    columns cannot be interpreted

121212 ksl
    Added -color option

130830 ksl
    Added -r option for default radius

150220 ksl
    Modified to use astropy for reading master files

251210 ksl
    Added support for FITS format input and Color column handling

.. moduleauthor:: KSL
"""

import sys
import os
from astropy.io import ascii
import numpy
from astropy.table import Table


def radec2deg(ra, dec):
    """Convert RA/Dec strings to decimal degrees.

    Handles both sexagesimal (HH:MM:SS, DD:MM:SS) and decimal degree formats.

    Parameters
    ----------
    ra : str or float
        Right Ascension as 'HH:MM:SS' string or decimal degrees.
    dec : str or float
        Declination as 'DD:MM:SS' string or decimal degrees.

    Returns
    -------
    tuple of (float, float)
        RA and Dec in decimal degrees.

    Notes
    -----
    If inputs are already floats, they are returned unchanged.
    RA in sexagesimal format is assumed to be in hours and is
    converted to degrees (multiplied by 15).
    """

    try:
        r = ra.split(":")
        d = dec.split(":")
    except AttributeError:
        return ra, dec

    rr = float(r[0])
    if len(r) > 1:
        rr = rr + float(r[1]) / 60.0
    if len(r) > 2:
        rr = rr + float(r[2]) / 3600.0
    if len(r) > 1:
        rr = 15.0 * rr  # Since we assume ra was in hms

    dd = float(d[0])
    x = 0
    if len(d) > 1:
        x = x + float(d[1]) / 60.0
    if len(d) > 2:
        x = x + float(d[2]) / 3600.0

    if dd > 0:
        dd = dd + x
    else:
        dd = dd - x

    return rr, dd


def read_masterfile(filename, xtype="circle", xmajor=3, xminor=3, xtheta=0.0):
    """Read a master table file and standardize column names.

    Reads ASCII or FITS format master tables and ensures all required
    columns are present, adding defaults where necessary.

    Parameters
    ----------
    filename : str
        Path to the master table file (ASCII or FITS format).
    xtype : str, optional
        Default region type if not specified in file. Default is 'circle'.
    xmajor : float, optional
        Default major axis/radius in arcseconds. Default is 3.
    xminor : float, optional
        Default minor axis in arcseconds. Default is 3.
    xtheta : float, optional
        Default position angle in degrees. Default is 0.0.

    Returns
    -------
    tuple of (astropy.table.Table, str)
        Standardized table with columns (Source_name, RA, Dec, RegType,
        Major, Minor, Theta, Color) and the coordinate type string.

    Notes
    -----
    The function handles various input formats:

    * Tables with full column headers
    * Generic tables with columns named col1, col2, etc.
    * Tables with only RA/Dec (source names auto-generated)
    * Tables with only positions and names (default geometry applied)
    """

    try:
        data = ascii.read(filename)
    except:
        try:
            data=Table.read(filename)
        except:
            print("Failed to read %s", filename)
            return 0

    colnames = data.colnames

    print(colnames)

    # Check whether the column names are defined or whether this a generic file

    if colnames[0] == "col1":
        if len(colnames) == 2:
            # Assume we have only a list of RA and Decs, in which
            # case we need to generate the source names
            data.rename_column("col1", "RA")
            data.rename_column("col2", "Dec")
            x = []
            i = 0
            while i < len(data["RA"]):
                x.append("x%03d" % i)
                i += 1
            data["Source_name"] = x
            # Add a check for nan's which cropped up in tweakreg files
            i = 0
            nan_row = []
            while i < len(data):
                if numpy.isnan(data["RA"][i]) or numpy.isnan(data["Dec"][i]):
                    print("gotcha %d" % i)
                    nan_row.append(i)
                i += 1
            if len(nan_row):
                data.remove_rows(nan_row)
        else:
            data.rename_column("col1", "Source_name")
            data.rename_column("col2", "RA")
            data.rename_column("col3", "Dec")
        if len(colnames) > 3:
            try:
                data.rename_column("col4", "RegType")
                data.rename_column("col5", "Major")
                data.rename_column("col6", "Minor")
                data.rename_column("col7", "Theta")
                complete = True
            except ValueError:
                complete = False
        else:
            complete = False
    else:  # This is an astropy table
        # Check that the essential files exist
        print("This is an astropy table with headers")
        ok = True
        if ("Source_name" in colnames) == False:
            print("There is no column named Source_name, so manufacturing one")
            names = []
            i = 0
            while i < len(data):
                names.append("x%03d" % (i + 1))
                i += 1
            data["Source_name"] = names
            # ok=False
        if ("RA" in colnames) == False:
            ok = False
        if ("Dec" in colnames) == False:
            ok = False
        if ok == False:
            print("Error: File read but column names are not correct:", colnames)
            print("       Minimally need Source_name,RA,Dec")
            return
        ok = True
        # print('test',colnames)
        # print('test','RegType' in colnames)

        if ("RegType" in colnames) == False:
            print("Did not get RegType")
            ok = False
        if ("Major" in colnames) == False:
            print("Did not get Major")
            ok = False
        if ("Minor" in colnames) == False:
            print("Did not get Minor")
            ok = False
        if ("Theta" in colnames) == False:
            print("Did not get Theta")
            ok = False

        if ("Color" in colnames) == False:
            print("Did not get  Color")
            data["Color"] = "Unknown"

        print("OK", ok)
        if ok == False:
            print("Colnames in table:\n", colnames)
            print(
                "Warning: Creating region file but with generic values for sizes:",
                colnames,
            )
            complete = False
        else:
            complete = True

    # So at this point we should have everything we need to continue

    print("Were we complete?", complete)

    if complete == False:
        # Create the region type
        data["RegType"] = xtype
        data["Major"] = xmajor
        data["Minor"] = 0.0
        data["Theta"] = 0.0
        # data["Color"] = "Unknown"

    # Make sure it is in the correct order
    xdata = data[
        "Source_name", "RA", "Dec", "RegType", "Major", "Minor", "Theta", "Color"
    ]

    print('XXX')
    print(xdata[:10])
    print('XXX')

    ttype = "unknown"
    return xdata, ttype


def write_regionfile(
    regionfile, records, source="unknown", type="unknown", color="red"
):
    """Write a DS9 region file from a master table.

    Creates a DS9-compatible region file with circles, ellipses, boxes,
    or annuli based on the input table records.

    Parameters
    ----------
    regionfile : str
        Output filename for the DS9 region file.
    records : astropy.table.Table
        Table containing region definitions with columns: Source_name,
        RA, Dec, RegType, Major, Minor, Theta, Color.
    source : str, optional
        Source identifier for header comment. Default is 'unknown'.
    type : str, optional
        Coordinate type ('fk5', 'physical', 'image'). Default is 'unknown'
        which outputs as 'fk5'.
    color : str, optional
        Default color for regions. Default is 'red'.

    Returns
    -------
    None
        Writes region file to disk.

    Notes
    -----
    If different rows have different colors in the Color column, each
    region will be written with its individual color. Otherwise, the
    default color is used for all regions.
    """

    # print(records[len(records)/2])

    # print('ok: got here')

    # First we need to understand how we are going to treat colors, namely whether
    # we need to write this out individually or collectively.  Our assumption is
    # that we will not change colors if the colors are different for different
    # lines of the master table

    records.info()

    xcolor = False
    i = 0
    while i < len(records):
        if records["Color"][i] != records["Color"][i - 1]:
            xcolor = True
            break
        i += 1

    if xcolor == True:
        print(
            "There are different colors for different region files, so we will preserve colors"
        )

    if color == "Unknown":
        color = "red"

    f = open(regionfile, "w")

    # write the header for the region file

    f.write("# Region file format: DS9 version 4.1\n")
    f.write("# Filename: %s\n" % regionfile)
    f.write(
        'global color=%s width=3 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        % color
    )
    if type != "unknown":
        f.write("%s\n" % type)
    else:
        f.write("fk5\n")

    i = 0
    while i < len(records):
        z = records[i]
        name = z["Source_name"]
        ra = z["RA"]
        dec = z["Dec"]
        regtype = z["RegType"]
        x1 = z["Major"]
        zcolor = z["Color"]
        if type == "physical" or type == "image":
            if z["RegType"] == "ellipse" and z["Minor"] > 0.0:
                x2 = z["Minor"]
                theta = z["Theta"]
                if xcolor:
                    f.write(
                        "ellipse(%f,%f,%.2f,%.2f,%.1f) # color=%s text={%s}\n"
                        % (ra, dec, x1, x2, theta, zcolor, name)
                    )
                else:
                    f.write(
                        "ellipse(%f,%f,%.2f,%.2f,%.1f) # text={%s}\n"
                        % (ra, dec, x1, x2, theta, name)
                    )
            else:
                if xcolor:
                    f.write(
                        "circle(%f,%f,%.2f)  # color=%s text={%s}\n"
                        % (ra, dec, x1, zcolor, name)
                    )
                else:
                    f.write("circle(%f,%f,%.2f)  # text={%s}\n" % (ra, dec, x1, name))
        else:  # We assume this is a normal masterfile with RA, DECs and sizes in arcsec
            if (z["RegType"] == "ellipse" or z["RegType"] == "box") and z["Minor"] > 0.0:
                zztype=z["RegType"]
                x2 = z["Minor"]
                theta = z["Theta"]
                try:
                    if xcolor:
                        f.write(
                            '%s(%f,%f,%.2f",%.2f",%.1f) # color=%s text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, zcolor, name)
                        )
                    else:
                        f.write(
                            '%s(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, name)
                        )
                except:
                    if xcolor:
                        f.write(
                            '%s(%s,%s,%.2f",%.2f",%.1f) # color=%s,text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, zcolor, name)
                        )
                    else:
                        f.write(
                            '%s(%s,%s,%.2f",%.2f",%.1f) # text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, name)
                        )
            elif z["RegType"] == "annulus":
                zztype=z["RegType"]
                x2 = z["Minor"]
                if xcolor:
                    f.write( '%s(%f,%f,%.2f",%.2f") # color=%s text={%s}\n'
                            % (zztype,ra, dec, x1, x2, zcolor, name))
                else:
                    f.write( '%s(%f,%f,%.2f",%.2f") #  text={%s}\n'
                            % (zztype,ra, dec, x1, x2, name))
            else:
                try:
                    if xcolor:
                        f.write(
                            'circle(%f,%f,%.2f")  # color=%s text={%s}\n'
                            % (ra, dec, x1, zcolor, name)
                        )
                    else:
                        f.write(
                            'circle(%f,%f,%.2f")  # text={%s}\n' % (ra, dec, x1, name)
                        )
                except TypeError:
                    if xcolor:
                        f.write(
                            'circle(%s,%s,%.2f")  # color=%s text={%s}\n'
                            % (ra, dec, x1, zcolor, name)
                        )
                    else:
                        f.write(
                            'circle(%s,%s,%.2f")  # text={%s}\n' % (ra, dec, x1, name)
                        )

        i = i + 1
    f.close()


def doit(argv):
    """Parse command line and execute master to region file conversion.

    Parameters
    ----------
    argv : list of str
        Command line arguments (typically sys.argv).

    Returns
    -------
    None
        Writes region file to disk.

    Notes
    -----
    Main driver function that:

    1. Parses command line options (-h, -r, -color)
    2. Reads the master table file
    3. Writes the DS9 region file
    """

    set_color = False
    color = "red"
    collated = "no"
    masterfile = ""
    regionfile = ""
    rad = 3.0

    i = 1
    while i < len(argv):
        if argv[i] == "-h":
            print(__doc__)
            return
        elif argv[i] == "-color":
            i = i + 1
            color = argv[i]
            set_color = True
        elif argv[i] == "-r":
            i = i + 1
            rad = eval(argv[i])
        elif argv[i] == "-collated":
            collated = "yes"
        else:
            if masterfile == "":
                masterfile = argv[i]
            elif regionfile == "":
                regionfile = argv[i]
        i = i + 1

    # Use the default region file name if the name is not yet assigned.
    if regionfile == "":
        regionfile = masterfile + ".reg"

    print("Masterfile          :", masterfile)
    print("Regionfile          :", regionfile)

    # Finished getting input information

    # read the masterfile

    data, ttype = read_masterfile(masterfile, xmajor=rad, xminor=rad)

    if data["Color"][0] == "Unknown":
        xcolor = color
    elif set_color:
        xcolor = color
    else:
        xcolor = data["Color"][0]

    # print('After return', data.colnames)
    # print('Got type',ttype)

    write_regionfile(regionfile, data, masterfile, ttype, xcolor)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":

    # Read the command line

    argc = len(sys.argv)

    if argc < 2:
        print(__doc__)
        sys.exit()

    doit(sys.argv)
