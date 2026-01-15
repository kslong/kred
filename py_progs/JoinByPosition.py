#!/usr/bin/env python
# coding: utf-8

"""JoinByPosition - Position-based table cross-matching

Space Telescope Science Institute

Synopsis
--------

Cross-match two tables based on sky position (RA, Dec). By default,
returns only matched rows. Use -all to preserve all rows from table1
(left join behavior).

Command Line Usage
------------------

::

    JoinByPosition.py [-h] [-sep ARCSEC] [-all] [-out OUTFILE]
                      [-ra1 COL] [-dec1 COL] [-ra2 COL] [-dec2 COL]
                      table1 table2

Arguments
---------

table1
    Path to the primary table

table2
    Path to the secondary table (matched rows joined to table1)

Options
-------

-h
    Print this help message and exit

-sep ARCSEC
    Maximum separation in arcseconds for a valid match. Default: 1.0

-all
    Return all rows from table1 (left join). Default: return only matches.

-out OUTFILE
    Output filename. Default: table1_x_table2.fits

-ra1 COL
    RA column name in table1. Default: RA

-dec1 COL
    Dec column name in table1. Default: Dec

-ra2 COL
    RA column name in table2. Default: RA

-dec2 COL
    Dec column name in table2. Default: Dec

Description
-----------

Uses KDTree algorithm for efficient cross-matching. For each object
in table1, finds the closest object in table2. If the separation is
less than the maximum allowed, the table2 columns are joined.

The output table contains:
- All columns from table1 (matched rows only, unless -all specified)
- All columns from table2 (with '_2' suffix if names conflict)
- Sep: separation in arcseconds

Notes
-----

History:

250115 ksl Coding begun, based on find_closest_objects from PhotCompare

"""

import sys
import os
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import KDTree


def read_table(filename):
    """
    Read a table from FITS or ASCII format.

    Parameters
    ----------
    filename : str
        Path to table file

    Returns
    -------
    Table
        Astropy Table object

    Raises
    ------
    IOError
        If file does not exist or cannot be read
    """
    if not os.path.isfile(filename):
        raise IOError('read_table: %s does not exist' % filename)

    try:
        xtable = Table.read(filename)
    except:
        try:
            xtable = ascii.read(filename)
        except:
            raise IOError('read_table: Could not read %s' % filename)

    return xtable


def do_join(table1, table2, max_sep=1.0, keep_all=False,
            ra1='RA', dec1='Dec', ra2='RA', dec2='Dec'):
    """
    Cross-match two tables based on sky position.

    Parameters
    ----------
    table1 : Table
        Primary table
    table2 : Table
        Secondary table (matched rows joined)
    max_sep : float, optional
        Maximum separation in arcseconds. Default: 1.0
    keep_all : bool, optional
        If True, keep all rows from table1 (left join).
        If False, keep only matched rows. Default: False
    ra1, dec1 : str, optional
        Column names for RA/Dec in table1. Default: 'RA', 'Dec'
    ra2, dec2 : str, optional
        Column names for RA/Dec in table2. Default: 'RA', 'Dec'

    Returns
    -------
    Table
        Joined table with Sep column added
    """
    # Validate columns exist
    if ra1 not in table1.colnames or dec1 not in table1.colnames:
        raise ValueError('table1 must have columns %s and %s' % (ra1, dec1))
    if ra2 not in table2.colnames or dec2 not in table2.colnames:
        raise ValueError('table2 must have columns %s and %s' % (ra2, dec2))

    # Convert to SkyCoord
    coords1 = SkyCoord(ra=table1[ra1] * u.degree, dec=table1[dec1] * u.degree)
    coords2 = SkyCoord(ra=table2[ra2] * u.degree, dec=table2[dec2] * u.degree)

    # Convert to Cartesian for KDTree
    cart1 = np.array([coords1.cartesian.x.value,
                      coords1.cartesian.y.value,
                      coords1.cartesian.z.value]).T
    cart2 = np.array([coords2.cartesian.x.value,
                      coords2.cartesian.y.value,
                      coords2.cartesian.z.value]).T

    # Build KDTree and query
    tree = KDTree(cart2)
    distances, indices = tree.query(cart1)

    # Get matched rows from table2
    matched_table2 = table2[indices]

    # Compute actual separations on sky
    matched_coords = SkyCoord(ra=matched_table2[ra2] * u.degree,
                              dec=matched_table2[dec2] * u.degree)
    separations = coords1.separation(matched_coords).arcsecond

    # Create mask for valid matches
    matched = separations < max_sep

    # Prepare table2 columns - rename duplicates
    table2_copy = matched_table2.copy()
    for col in table2_copy.colnames:
        if col in table1.colnames:
            table2_copy.rename_column(col, col + '_2')

    # Add separation column
    sep_col = separations.copy()

    # Build output table
    result = table1.copy()
    result['Sep'] = sep_col
    result['Sep'].format = '.3f'

    # Join with table2 columns
    result = hstack([result, table2_copy])

    if keep_all:
        # Set unmatched rows to NaN
        for col in table2_copy.colnames:
            if result[col].dtype.kind in ('f', 'i'):
                result[col] = result[col].astype(float)
                result[col][~matched] = np.nan
            elif result[col].dtype.kind in ('U', 'S', 'O'):
                result[col][~matched] = ''
        result['Sep'][~matched] = np.nan
    else:
        # Keep only matched rows
        result = result[matched]

    return result


def do_one(table1_path, table2_path, max_sep=1.0, keep_all=False,
           outfile='', ra1='RA', dec1='Dec', ra2='RA', dec2='Dec'):
    """
    Cross-match two table files and write result.

    Parameters
    ----------
    table1_path : str
        Path to primary table
    table2_path : str
        Path to secondary table
    max_sep : float, optional
        Maximum separation in arcseconds. Default: 1.0
    keep_all : bool, optional
        If True, keep all rows from table1. Default: False
    outfile : str, optional
        Output filename. If empty, auto-generated.
    ra1, dec1, ra2, dec2 : str, optional
        Column names for coordinates

    Returns
    -------
    Table
        The joined table
    """
    # Read tables
    print('Reading %s' % table1_path)
    table1 = read_table(table1_path)
    print('  %d rows' % len(table1))

    print('Reading %s' % table2_path)
    table2 = read_table(table2_path)
    print('  %d rows' % len(table2))

    # Perform join
    print('Cross-matching with max_sep = %.2f arcsec' % max_sep)
    result = do_join(table1, table2, max_sep=max_sep, keep_all=keep_all,
                     ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2)

    if keep_all:
        n_matched = np.sum(np.isfinite(result['Sep']))
        print('Matched %d of %d rows (%.1f%%), returning all rows' %
              (n_matched, len(table1), 100.0 * n_matched / len(table1)))
    else:
        print('Returning %d matched rows of %d in table1 (%.1f%%)' %
              (len(result), len(table1), 100.0 * len(result) / len(table1)))

    # Generate output filename if not provided
    if outfile == '':
        name1 = os.path.basename(table1_path).replace('.fits', '').replace('.txt', '')
        name2 = os.path.basename(table2_path).replace('.fits', '').replace('.txt', '')
        outfile = '%s_x_%s.fits' % (name1, name2)

    # Ensure .fits extension
    if not outfile.endswith('.fits'):
        outfile = outfile + '.fits'

    # Write output
    result.write(outfile, format='fits', overwrite=True)
    print('Wrote %s' % outfile)

    return result


def steer(argv):
    """
    Parse command line arguments and execute cross-match.

    Parameters
    ----------
    argv : list
        Command line arguments (sys.argv)
    """
    max_sep = 1.0
    keep_all = False
    outfile = ''
    ra1 = 'RA'
    dec1 = 'Dec'
    ra2 = 'RA'
    dec2 = 'Dec'
    table1_path = ''
    table2_path = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i][:4] == '-sep':
            i += 1
            max_sep = float(argv[i])
        elif argv[i][:4] == '-all':
            keep_all = True
        elif argv[i][:4] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i][:4] == '-ra1':
            i += 1
            ra1 = argv[i]
        elif argv[i][:5] == '-dec1':
            i += 1
            dec1 = argv[i]
        elif argv[i][:4] == '-ra2':
            i += 1
            ra2 = argv[i]
        elif argv[i][:5] == '-dec2':
            i += 1
            dec2 = argv[i]
        elif argv[i][0] == '-':
            print('Error: Unknown option: %s' % argv[i])
            return
        elif table1_path == '':
            table1_path = argv[i]
        elif table2_path == '':
            table2_path = argv[i]
        else:
            print('Error: Too many arguments: %s' % argv[i])
            return
        i += 1

    if table1_path == '' or table2_path == '':
        print('Error: Must provide two table files')
        print('Usage: JoinByPosition.py [-sep ARCSEC] [-all] table1 table2')
        return

    do_one(table1_path, table2_path, max_sep=max_sep, keep_all=keep_all,
           outfile=outfile, ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
