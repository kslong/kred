#!/usr/bin/env python 

"""
Space Telescope Science Institute - Tile Setup Module

This module sets up directories and identifies files required for initial
processing steps to run PrepFiles on DECam data.

Overview
--------

The module reads tile definition files and image definition files to organize
exposures for processing. It searches for CCD images that should be prepped
to construct images in tiles of a field, then writes the results to tables
in the Summary directory.

Command Line Usage
------------------
::

    SetupTile.py [-h] [-all] [-S7] [-seeing_max 1.3] [-xtab myfile.txt] 
                 [-ptab my_image.txt] [-use_all_data] field tile1 tile2 ...

**Required Arguments:**

field
    Name of a field (e.g., LMC_c42)
    
tile1, tile2, ...
    Names of tiles as specified in the MC_tiles.txt file

**Optional Arguments:**

-h
    Print help message
    
-all
    Create tables containing relevant information for all tiles in a field.
    Note: This is DIFFERENT from some earlier routines where -all refers
    to all files.
    
-xtab myfile.txt
    Use a different tile definition file
    
-ptab my_image.txt
    Use a different filter exposure table for creating images. This table
    allows combining images with different exposure times for a filter.
    If a filter is missing, that filter will not be processed.
    
-S7
    Ignore the S7 chip
    
-seeing_max VALUE
    Ignore images with seeing greater than VALUE
    
-use_all_data
    Include all data regardless of field (as long as that field has been
    MefPrepped)

Output
------
The program creates one table per tile in the Summary directory with filenames
like ``LMC_c30_T16.txt``. These tables contain only files which can be 
processed further, excluding:

* Images with seeing larger than the specified threshold
* Images from the problematic S7 chip (if requested)
* Images that have not been MefPrepped

Configuration Files
-------------------

The program searches for configuration files in:

1. Current directory
2. ``$KRED/config`` directory

Default files:

* Tile definition: ``MC_tiles.txt``
* Image definition: ``DeMCELS_images.txt``

Notes
-----
* If the directory exists where links to files will be placed, all existing
  file links are removed so only files from the current run remain.
* The routine reads tables in the Summary directory to decide which images
  to include.
* All relevant fields must be processed through MefPrep before running this
  script.

Examples
--------
Process a single tile::

    $ SetupTile.py LMC_c42 T07

Process multiple tiles with seeing constraint::

    $ SetupTile.py -seeing_max 1.3 LMC_c42 T07 T08 T09

Process all tiles in a field, excluding S7 chip::

    $ SetupTile.py -all -S7 LMC_c42

Use all available data from any field::

    $ SetupTile.py -use_all_data LMC_c42 T07

History
-------
230513 ksl
    Coding begun
"""

import sys
from astropy.io import ascii
import numpy as np
import os
from astropy.table import vstack, join
from glob import glob

from astropy.io import ascii
from log import *

CWD = os.getcwd()
DATADIR = '%s/DECam_CCD/' % (CWD)
PREPDIR = '%s/DECam_PREP/' % (CWD)


def find_overlapping_images(field='LMC_c45', tile='T07', ra_center=81.1, 
                           dec_center=-66.17, size_deg=.1):
    """
    Find all CCD images with pixels within the final image area.
    
    This routine identifies all possible CCD images that have pixels within
    the specified region. It does not perform additional checks and does not
    eliminate images that have not been processed by MefPrep.
    
    Parameters
    ----------
    field : str, optional
        Name of the field (e.g., 'LMC_c45'). Default: 'LMC_c45'.
    tile : str, optional
        Tile identifier (e.g., 'T07'). Default: 'T07'.
    ra_center : float, optional
        Right ascension of the tile center in degrees. Default: 81.1.
    dec_center : float, optional
        Declination of the tile center in degrees. Default: -66.17.
    size_deg : float, optional
        Size of the tile in degrees. Default: 0.1.
    
    Returns
    -------
    Table or empty list
        Astropy Table containing information about overlapping images with
        columns including Field, Filename, Root, EXTNAME, coordinates,
        FILTER, EXPTIME, MAGZERO, and SEEING. Returns empty list if no
        files are found or if required input files are missing.
    
    Notes
    -----
    The function performs a geometric search by:
    
    1. Reading detector and MEF summary tables for the field
    2. Computing RA/Dec boundaries accounting for cosine declination
    3. Filtering images whose corners overlap the tile region
    4. Computing offsets from tile center for each image
    
    The search accounts for spherical geometry when computing RA boundaries.
    
    Examples
    --------
    >>> images = find_overlapping_images('LMC_c42', 'T07', 81.1, -66.2, 0.67)
    >>> print(f"Found {len(images)} overlapping images")
    """
    
    tab_file = 'Summary/%s_det.tab' % (field)
    if os.path.isfile(tab_file):
        det = ascii.read(tab_file)
    else:
        print('SetupTile: Cannot find %s' % tab_file)
        return []

    mef_file = 'Summary/%s_mef.tab' % (field)

    if os.path.isfile(mef_file):
        mef = ascii.read(mef_file)
    else:
        print('SetupTile: Cannot find %s' % mef_file)
        return []

    x = join(det, mef, join_type='left')

    print('SetupTile: Looking in %s for CCD images with RA and DEC of %.5f %.5f and size of %.2f deg' % 
          (field, ra_center, dec_center, size_deg))

    # Compute search boundaries
    dec_min = dec_center - 0.5 * size_deg
    dec_max = dec_center + 0.5 * size_deg
    delta_ra = size_deg / np.cos(dec_center / 57.29578)
    ra_min = ra_center - 0.5 * delta_ra
    ra_max = ra_center + 0.5 * delta_ra
    
    # Filter by Dec boundaries
    x['Dec_max'] = np.max([x['COR1DEC1'], x['COR2DEC1'], x['COR3DEC1'], x['COR4DEC1']], axis=0)
    x['Dec_min'] = np.min([x['COR1DEC1'], x['COR2DEC1'], x['COR3DEC1'], x['COR4DEC1']], axis=0)
    
    x = x[x['Dec_min'] < dec_max]
    if len(x) == 0:
        print('SetupTile:   Failed to find any images with Dec less than %.5f ' % dec_max)
        return []

    x = x[x['Dec_max'] > dec_min]
    if len(x) == 0:
        print('SetupTile:   Failed to find any images with Dec greater than %.5f ' % dec_min)
        return []

    # Filter by RA boundaries
    x['RA_max'] = np.max([x['COR1RA1'], x['COR2RA1'], x['COR3RA1'], x['COR4RA1']], axis=0)
    x['RA_min'] = np.min([x['COR1RA1'], x['COR2RA1'], x['COR3RA1'], x['COR4RA1']], axis=0)
    
    x = x[x['RA_min'] < ra_max]
    if len(x) == 0:
        print('SetupTile:   Failed to find any images with RA less than %.5f ' % ra_max)
        return []

    x = x[x['RA_max'] > ra_min]
    if len(x) == 0:
        print('SetupTile:   Failed to find any images with RA greater than %.5f ' % ra_min)
        return []

    # Compute offsets from tile center
    x['Delta_Dec'] = x['CENDEC1'] - dec_center
    x['Delta_RA'] = (x['CENRA1'] - ra_center) / np.cos(dec_center / 57.29578)
    
    x['Delta_Dec'].format = '.2f'
    x['Delta_RA'].format = '.2f'

    # Create filename list
    xfiles = []
    for one in x:
        one_name = '%s_%s.fits' % (one['Root'], one['EXTNAME'])
        xfiles.append(one_name)

    x['Filename'] = xfiles
    x['Field'] = field

    print('SetupTile:   Found %d relevant files in %s' % (len(x), field))
    
    return x['Field', 'Filename', 'Root', 'EXTNAME', 'CENRA1', 'CENDEC1', 
             'COR1RA1', 'COR1DEC1', 'COR2RA1', 'COR2DEC1', 'COR3RA1', 'COR3DEC1', 
             'COR4RA1', 'COR4DEC1', 'Delta_Dec', 'Delta_RA',
             'FILTER', 'EXPTIME', 'MAGZERO', 'SEEING']


def locate_prepped_files(xtab, field, tile):
    """
    Check whether necessary files have been processed through MefPrep.
    
    This function verifies that files identified as overlapping the tile
    region actually exist in the expected data directory after MefPrep
    processing.
    
    Parameters
    ----------
    xtab : Table
        Astropy Table containing file information with at least 'Field' and
        'Filename' columns.
    field : str
        Name of the field.
    tile : str
        Tile identifier.
    
    Returns
    -------
    Table
        Modified input table with added columns:
        
        * Prepped : str - 'Yes' if file exists, 'No' if missing
        * PrepFile : str - Full path to the prepped file
    
    Raises
    ------
    IOError
        If the data directory for the field does not exist.
    
    Notes
    -----
    This function adds status information to help identify which files are
    ready for further processing. Files marked 'No' typically indicate that
    MefPrep has not been run on all relevant fields.
    
    Examples
    --------
    >>> xtab = find_overlapping_images('LMC_c42', 'T07', 81.1, -66.2, 0.67)
    >>> xtab = locate_prepped_files(xtab, 'LMC_c42', 'T07')
    >>> ready = xtab[xtab['Prepped'] == 'Yes']
    >>> print(f"{len(ready)} files are ready for processing")
    """
    
    data_dir = '%s/%s/data' % (DATADIR, field)

    if os.path.isdir(data_dir) == False:
        print('SetupTile: Error: %s does not appear to exist' % data_dir)
        print('Run PrepMef on this field first')
        raise IOError

    xtab['Prepped'] = 'Unknown'
    xlocation = []

    nerrors = 0
    nfound = 0
    for one in xtab:
        one_dir = '%s/%s/data' % (DATADIR, one['Field'])
        one_file = one['Filename']
        data_file = '%s/%s' % (one_dir, one_file)
        xlocation.append(data_file)
        if os.path.isfile(data_file) == False:
            one['Prepped'] = 'No'
            nerrors += 1
        else:
            one['Prepped'] = 'Yes'
            nfound += 1
            
    print('SetupTile: Successfully found %d files in %s' % (nfound, data_dir))
    if nerrors:
        print('SetupTile: Failed to find %d of %d files in %s' % (nerrors, len(xtab), data_dir))
        print('SetupTile: Normally this is because MefPrep was not run on all relevant Fields')

    xtab['PrepFile'] = xlocation
    return xtab


def populate_tile_dir(xtab, field, tile):
    """
    Create tile directory and populate with symbolic links to data files.
    
    This function creates the tile-specific directory structure and populates
    it with symbolic links to the prepped FITS files. If the directory already
    exists, existing FITS files and links are removed before creating new ones.
    
    Parameters
    ----------
    xtab : Table
        Astropy Table containing file information with 'Field' and 'Filename'
        columns.
    field : str
        Name of the field.
    tile : str
        Tile identifier.
    
    Raises
    ------
    IOError
        If the data directory does not exist or if required files cannot be
        found (typically because MefPrep was not run).
    
    Notes
    -----
    **Directory Structure:**
    
    Creates: ``DECam_PREP/{field}/{tile}/``
    
    **Symbolic Links:**
    
    For each file in xtab, creates a symlink from the tile directory to the
    actual data file in ``DECam_CCD/{field}/data/``
    
    **Cleanup:**
    
    If the tile directory already exists, all existing FITS files and links
    are removed to ensure a clean state.
    
    Examples
    --------
    >>> xtab = ascii.read('Summary/LMC_c42_T07.txt')
    >>> populate_tile_dir(xtab, 'LMC_c42', 'T07')
    SetupTile: Successfully linked 145 files from DECam_CCD/LMC_c42/data to DECam_PREP/LMC_c42/T07/
    """
    
    data_dir = '%s/%s/data' % (DATADIR, field)

    if os.path.isdir(data_dir) == False:
        print('SetupTile: Error: %s does not appear to exist' % data_dir)
        print('Run PrepMef on this field first')
        raise IOError

    tile_dir = '%s/%s/%s/' % (PREPDIR, field, tile)

    if os.path.isdir(tile_dir) == False:
        print('SetupTile: Creating the Prep directory as %s' % tile_dir)
        os.makedirs(tile_dir)
    else:
        print('SetupTile: The Prep directory %s already exists' % tile_dir)
        xfiles = glob('%s/*fits*' % tile_dir)
        if len(xfiles):
            print('SetupTile: Removing existing fits files or links and reinitializing')
            for one in xfiles:
                os.remove(one)

    nerrors = 0
    nfound = 0
    for one in xtab:
        one_dir = '%s/%s/data' % (DATADIR, one['Field'])
        one_file = one['Filename']
        data_file = '%s/%s' % (one_dir, one_file)
        tile_file = '%s/%s' % (tile_dir, one_file)
        if os.path.isfile(data_file) == False:
            print('SetupTile: Failed to find %s in %s (Was MefPrep run?)' % (one_file, data_dir))
            nerrors += 1
        else:
            nfound += 1
            os.symlink(data_file, tile_file)
            
    print('SetupTile: Successfully linked %d files from %s to %s' % (nfound, data_dir, tile_dir))
    if nerrors:
        print('SetupTile: Failed to find %d of %d files in %s' % (nerrors, len(xtab), data_dir))
        print('SetupTile: Normally this is because MefPrep was not run on all relevant Fields')
        raise IOError


def get_tile_files(field='LMC_c42', tile='T07', ra=81.108313, dec=-66.177280, 
                   size_deg=0.67, s7=True, seeing_max=1000., use_all_data=False):
    """
    Find and filter chip images for a specific tile.
    
    This is the main workhorse function that identifies all chip images that
    should be used for a specific tile, applies quality filters, and writes
    the results to a summary table.
    
    Parameters
    ----------
    field : str, optional
        Name of the field (e.g., 'LMC_c42'). Default: 'LMC_c42'.
    tile : str, optional
        Tile identifier (e.g., 'T07'). Default: 'T07'.
    ra : float, optional
        Right ascension of tile center in degrees. Default: 81.108313.
    dec : float, optional
        Declination of tile center in degrees. Default: -66.177280.
    size_deg : float, optional
        Size of tile in degrees. Default: 0.67.
    s7 : bool, optional
        If False, eliminate the S7 chip. Default: True.
    seeing_max : float, optional
        Maximum allowed seeing in arcseconds. Images with worse seeing are
        excluded. Default: 1000.0 (effectively no limit).
    use_all_data : bool, optional
        If True, search for relevant data in all fields of the same galaxy,
        not just the specified field. Default: False.
    
    Returns
    -------
    str
        Filename of the output table written to the Summary directory.
    
    Raises
    ------
    IOError
        If no suitable files are found for the field/tile combination.
    
    Notes
    -----
    **Processing Steps:**
    
    1. Find all overlapping images (optionally from multiple fields)
    2. Check which images have been MefPrepped
    3. Filter out S7 chip if requested
    4. Filter out images with poor seeing if threshold is set
    5. Write results to ``Summary/{field}_{tile}.txt``
    
    **Output Table Columns:**
    
    The output table includes all columns from find_overlapping_images plus
    Prepped status and PrepFile path.
    
    Examples
    --------
    >>> # Standard single-field processing
    >>> outfile = get_tile_files('LMC_c42', 'T07', 81.1, -66.2, 0.67)
    
    >>> # Exclude S7 and limit seeing
    >>> outfile = get_tile_files('LMC_c42', 'T07', 81.1, -66.2, 0.67, 
    ...                          s7=False, seeing_max=1.3)
    
    >>> # Use all available galaxy data
    >>> outfile = get_tile_files('LMC_c42', 'T07', 81.1, -66.2, 0.67,
    ...                          use_all_data=True)
    """
        
    if use_all_data == True:
        words = field.split('_')
        galaxy = words[0]
        xfiles = glob('Summary/%s*det.tab' % galaxy)
        print('\n\nSetupTile: There are %d fields to survey' % len(xfiles))
        fields = []
        for one_file in xfiles:
            xx = one_file.replace('Summary/', '')
            xx = xx.replace('_det.tab', '')
            fields.append(xx)
        z = []
        for one in fields:
            one_result = find_overlapping_images(one, tile, ra, dec, size_deg)
            if len(one_result):
                one_result = locate_prepped_files(one_result, field, tile)
                z.append(one_result)
        x = vstack(z)
    else:
        x = find_overlapping_images(field, tile, ra, dec, size_deg)
        if len(x):
            x = locate_prepped_files(x, field, tile)

    # Eliminate files that were not prepped
    nn = len(x)
    x = x[x['Prepped'] == 'Yes']
    print('\n\nSetupTile: Of %d images that overlapped the field, %d have been processed with MefPrep' % 
          (nn, len(x)))

    if len(x) == 0:
        print('SetupTile: Error: no files found for field/tile  %s %s' % (field, tile))
        raise IOError

    # Apply quality filters
    original_length = len(x)
    delta_s7 = 0
    delta_seeing = 0
    
    if s7 == False:
        foo = x[x['EXTNAME'] == 'S7']
        delta_s7 = len(foo)
        x = x[x['EXTNAME'] != 'S7']
        
    if seeing_max < 100.:
        foo = x[x['SEEING'] > seeing_max]
        delta_seeing = len(foo)
        x = x[x['SEEING'] <= seeing_max]
        
    if len(x) < original_length:
        print('SetupTile: Eliminated %d s7 images and %d bad seeing images of %d that were possible' % 
              (delta_s7, delta_seeing, original_length))
    else:
        print('SetupTile: Using all %d images possible for this tile' % original_length)

    x.meta['comments'] = ['RA %f' % ra, 'DEC %f' % dec]

    if len(x) > 0:
        outfile = 'Summary/%s_%s.txt' % (field, tile)
        x.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)
    else:
        print('SetupTile: Failed to set up %s %s' % (field, tile))
        raise IOError

    return outfile


def setup_tiles(ztab, s7, seeing_max, use_all_data, process_tab): 
    """
    Setup one or more tiles for processing with Swarp and other tools.
    
    This function orchestrates the complete setup process for multiple tiles,
    including file identification, filtering, and directory population.
    
    Parameters
    ----------
    ztab : Table
        Astropy Table containing tile definitions with columns: Field, Tile,
        RA, Dec, and Size.
    s7 : bool
        If True, include chip S7; if False, exclude it.
    seeing_max : float
        Maximum allowed seeing in arcseconds. Values >= 100 effectively
        disable this filter.
    use_all_data : bool
        If True, search for relevant data in all fields of the same galaxy.
    process_tab : Table
        Table defining image processing parameters, typically containing
        filter and exposure time information.
    
    Returns
    -------
    None
        Results are written to files and directories.
    
    Notes
    -----
    **Processing for Each Tile:**
    
    1. Call get_tile_files() to identify and filter images
    2. Join with processing table to add processing parameters
    3. Write updated summary table
    4. Call populate_tile_dir() to create symlinks
    
    **Error Handling:**
    
    If any tile fails, an error message is printed suggesting to run
    ``MefPrep.py -finish`` before retrying.
    
    Examples
    --------
    >>> # Assuming ztab and ptab are loaded
    >>> setup_tiles(ztab, s7=False, seeing_max=1.3, 
    ...            use_all_data=False, process_tab=ptab)
    """
    
    for one in ztab:
        field = one['Field']
        tile = one['Tile']
        try:
            outfile_name = get_tile_files(one['Field'], one['Tile'], one['RA'], 
                                         one['Dec'], one['Size'], s7, seeing_max, 
                                         use_all_data)
        except IOError:
            print('SetupTile: Error: Something is wrong, and must be sorted before continuing')
            print('SetupTile: Try rerunning MefPrep.py -finish, and then repeat this step')
            return

        xtab = ascii.read(outfile_name)
        xtab = join(process_tab, xtab, join_type='left')
        
        # Eliminate lines with no match
        xtab = xtab[~xtab['Filename'].mask]
        xtab.write(outfile_name, format='ascii.fixed_width_two_line', overwrite=True)
        
        try:
            populate_tile_dir(xtab, field, tile)
        except IOError:
            print('SetupTile: Problem with populating tile dir for %s %s' % (field, tile))

    return


def steer(argv):
    """
    Parse command-line arguments and execute tile setup.
    
    This is the main entry point that handles command-line argument parsing
    and orchestrates the tile setup process.
    
    Parameters
    ----------
    argv : list
        Command-line argument list (typically sys.argv).
    
    Returns
    -------
    None or int
        Returns without value on success, returns 0 on certain error conditions.
    
    Command-Line Arguments
    ----------------------
    See module docstring for complete argument documentation.
    
    Configuration Files
    -------------------
    Default files (searched in current directory, then $KRED/config):
    
    * MC_tiles.txt - Tile definitions
    * DeMCELS_images.txt - Image processing parameters
    
    Notes
    -----
    **Execution Flow:**
    
    1. Parse command-line arguments
    2. Load configuration files
    3. Filter tiles based on field and tile names
    4. Call setup_tiles() to perform setup
    5. Log operations to field-specific log files
    
    **Logging:**
    
    Creates log entries in ``{field}.log`` for each tile processed.
    
    Examples
    --------
    Command line usage::
    
        $ python SetupTile.py LMC_c42 T07 T08
        $ python SetupTile.py -all -S7 -seeing_max 1.3 LMC_c42
        $ python SetupTile.py -use_all_data LMC_c42 T07
    """
    
    field = ''
    tiles = []
    xall = False
    table = ''
    process_table = ''
    seeing_max = 1000.
    use_s7 = True
    use_all_data = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-all':
            xall = True
        elif argv[i] == '-S7' or argv[i] == '-s7':
            use_s7 = False
        elif argv[i] == '-use_all_data':
            use_all_data = True 
        elif argv[i] == '-xtab':
            i += 1
            table = argv[i]
        elif argv[i] == '-ptab':
            i += 1
            process_table = argv[i]
        elif argv[i] == '-seeing_max':
            i += 1
            seeing_max = eval(argv[i])
        elif argv[i][0] == '-':
            print('SetupTile: Error: Unknown switch %s ' % argv[i])
            return
        elif field == '':
            field = argv[i]
        else:
            if field == '':
                field = argv[i]
            else:
                tiles.append(argv[i])
        i += 1

    if field == '':
        print('SetupTile: Sorry: there seems to be nothing to do')
        print('SetupTile: A field must be provided in the command line')
        return

    if xall == False and len(tiles) == 0:
        print('SetupTile: Sorry: there seems to be nothing to do')
        print('-all not set and no tiles to set up provided')
        return

    # Set default configuration files
    if table == '':
        table = 'MC_tiles.txt'
    if process_table == '':
        process_table = 'DeMCELS_images.txt'

    # Load tile configuration
    kred = os.getenv('KRED')
    if os.path.isfile(table):
        xtab = ascii.read(table)
    elif os.path.isfile(kred + '/config/' + table):
        xtab = ascii.read(kred + '/config/' + table)
    else:
        print('SetupTile: Error: Could not find tile config %s in local directory or in kred/config' % table)
        return

    # Load image processing configuration
    if os.path.isfile(process_table):
        ptab = ascii.read(process_table)
        print('Read local process_table')
    elif os.path.isfile(kred + '/config/' + process_table):
        ptab = ascii.read(kred + '/config/' + process_table)
        print('Read process_table in kred/config directory')
    else:
        print('SetupTile: Error: Could not find image config %s in local directory or in kred/config' % process_table)
        return

    # Filter for requested field
    xtab = xtab[xtab['Field'] == field]

    if len(xtab) == 0:
        print('SetupTile: Could not find field %s in %s' % (field, table))
        return 0

    # Select requested tiles
    if xall == False:
        i = 0         
        for one_tile in tiles:
            q = xtab[xtab['Tile'] == one_tile]
            if i == 0:
                xtiles = q.copy()
                i += 1
            else:
                xtiles = vstack([xtiles, q])
    else:
        xtiles = xtab

    if len(xtiles) == 0:
        print('SetupTile: Sorry: there seems to be nothing to do')
        print('Looked for the following tiles: ', tiles)
        return

    # Perform tile setup
    setup_tiles(xtiles, use_s7, seeing_max, use_all_data, ptab)

    # Log operations
    fields = np.unique(xtiles['Field'])
    for one in fields:
        open_log('%s.log' % one)
        foo = xtiles[xtiles['Field'] == one]
        for one_tile in foo:
            log_message('SetupTile %s %s' % (one_tile['Field'], one_tile['Tile']))
        close_log()

    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
