#!/usr/bin/env python
"""BackPrep - Background Preparation using SWarp.

Space Telescope Science Institute

Synopsis
--------


Create background images from DECam data using SWarp for improved
source detection and photometry.

Command Line Usage
------------------

::

    Basic usage::

    BackPrep.py [-all] [-setup] [-run] [-indir directory] field tile

    Arguments
    ^^^^^^^^^

    field : str
    Field name (e.g., 'LMC_c42')
    tile : str
    Tile identifier (e.g., 'T07')

    Options
    ^^^^^^^

    -all
    Process all tiles found for the specified field
    -setup
    Only set up SWarp configuration files without running (default behavior)
    -run
    Set up and execute all SWarp commands
    -redo
    Regenerate background files even if they already exist
    -indir directory
    Use a different input directory (e.g., 'DECam_PREP2' for diagnostics).
    Outputs will be placed in a 'XXX' subdirectory.

Description
-----------


BackPrep.py performs two main tasks:

1. **Setup mode** (default): Creates SWarp configuration files and run scripts
   for background image generation without executing them.

2. **Run mode** (with -run): Sets up configuration files and executes all
   SWarp commands to generate the background images.

The module processes astronomical images organized by field and tile, creating
resampled background images that can be used for background subtraction and
source matching.

Version History
---------------

230504 ksl
    Coding begun

Examples
--------

Set up SWarp for a single tile::

    $ python BackPrep.py LMC_c42 T07

Set up and run SWarp for a single tile::

    $ python BackPrep.py -run LMC_c42 T07

Process all tiles in a field::

    $ python BackPrep.py -all -run LMC_c42

Use a custom input directory::

    $ python BackPrep.py -indir DECam_PREP2 -run LMC_c42 T07

Version History
---------------


230504 ksl

    Coding begun

Examples
--------

Set up SWarp for a single tile::

    $ python BackPrep.py LMC_c42 T07

Set up and run SWarp for a single tile::

    $ python BackPrep.py -run LMC_c42 T07

Process all tiles in a field::

    $ python BackPrep.py -all -run LMC_c42

Use a custom input directory::

    $ python BackPrep.py -indir DECam_PREP2 -run LMC_c42 T07
"""


import os
import stat

from glob import glob
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats
from log import *

import timeit
import subprocess

SWARPDIR = 'DECam_BACK'
PREPDIR = os.path.abspath('DECam_PREP')
exec_dir = os.getcwd()


def create_swarp_dir(field='LMC_c42', tile='T07'):
    """Create a directory for SWarp outputs if it does not exist.

    Parameters
    ----------
    field : str, optional
        Field name. Default is 'LMC_c42'.
    tile : str, optional
        Tile identifier. Default is 'T07'.

    Returns
    -------
    outdir : str
        Path to the created or existing output directory.

    """
    outdir = '%s/%s/%s/' % (SWARPDIR, field, tile)
    if os.path.isdir(outdir) == False:
        os.makedirs(outdir)
    return outdir


# Default SWarp configuration template
default = '''
# Default configuration file for SWarp 2.38.0
# EB 2019-07-10
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME         foo.fits      # Output filename
WEIGHTOUT_NAME      foo.weight.fits # Output weight-map filename

HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
CENTER                %.4f, %.4f       # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            2            # Pixel scale
IMAGE_SIZE          2000,2000          # Image size (0 = AUTOMATIC)

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          N               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE             2048             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)

#------------------------------ Miscellaneous ---------------------------------
DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
'''


def prep_one_tile(field='LMC_c42', tile='T07', redo=False):
    """Prepare SWarp configuration and run script for one tile.

    This function reads the summary file for a field/tile combination,
    creates the necessary output directories, generates a SWarp configuration
    file, and creates a run script with SWarp commands for each image.

    Parameters
    ----------
    field : str, optional
        Field name. Default is 'LMC_c42'.
    tile : str, optional
        Tile identifier. Default is 'T07'.
    redo : bool, optional
        If True, regenerate all background files even if they exist.
        If False (default), skip files that already exist.

    Returns
    -------
    n : int
        Number of background files that need to be created.

    Raises
    ------
    IOError
        If the summary file or input directory cannot be found.

    Notes
    -----
    The function creates:
    - Output directory: ``SWARPDIR/field/tile/``
    - SWarp config: ``SWARPDIR/field/tile/swarp.default``
    - Run script: ``SWARPDIR/field/tile/RunSwarp``

    The run script is made executable and contains one SWarp command per
    input image.

    Examples
    --------
    Prepare configuration for a single tile::

        >>> n = prep_one_tile('LMC_c42', 'T07')
        >>> print(f"{n} files need processing")

    Force regeneration of all files::

        >>> n = prep_one_tile('LMC_c42', 'T07', redo=True)

    """
    try:
        sumfile = 'Summary/%s_%s.txt' % (field, tile)
        x = ascii.read(sumfile)
    except:
        print('BackPrep: Could not find %s' % sumfile)
        raise IOError

    indir = '%s/%s/%s' % (PREPDIR, field, tile)
    if os.path.exists(indir) == False:
        print('BackPrep: Could not find the directory %s' % indir)
        raise IOError

    # Now create a directory if it does not exist
    outdir = '%s/%s/%s' % (SWARPDIR, field, tile)
    if os.path.exists(outdir) == False:
        os.makedirs(outdir)
        print('BackPrep: The output directory did not exist, so it is being created')
    else:
        print('BackPrep: The output directory exists')

    ra = np.average(x['CENRA1'])
    dec = np.average(x['CENDEC1'])

    # This needs to be relative to the outdir
    indir = '%s/%s/%s/' % (PREPDIR, field, tile)

    f = open('%s/swarp.default' % outdir, 'w')
    f.write(default % (ra, dec))
    f.close()

    runfile = '%s/RunSwarp' % outdir
    f = open(runfile, 'w')
    n = 0
    for one in x:
        root = one['Filename'].replace('.fits', '')
        if redo == False:
            xfile = '%s/%s.fits' % (outdir, root)
            if os.path.isfile(xfile):
                # print('Found %s already, so skipping' % xfile)
                continue
        f.write('swarp %s/%s.fits -IMAGEOUT_NAME %s.fits  -WEIGHTOUT_NAME  /dev/null -c swarp.default\n' % (indir, root, root))
        n += 1
        # gen_one_swarp_command(one['Filename'],ra,dec,indir,outdir)
    f.close()
    os.chmod(runfile, stat.S_IRWXU)

    print('BackPrep: There are %d background files of %d for field %s and tile %s that need creating' % (n, len(x), field, tile))

    return n


def run_back(field='LMC_c42', tile='T07'):
    """Execute SWarp for background field processing.

    This function changes to the appropriate output directory and executes
    the RunSwarp script to generate background images using SWarp.

    Parameters
    ----------
    field : str, optional
        Field name. Default is 'LMC_c42'.
    tile : str, optional
        Tile identifier. Default is 'T07'.

    Returns
    -------
    None

    Notes
    -----
    The function:
    - Changes to the SWarp output directory
    - Executes the RunSwarp script
    - Captures output to ``Swarp.out.txt``
    - Prints the last 10 lines of output
    - Reports total execution time
    - Returns to the original execution directory

    The function filters out progress indicators from SWarp output to
    keep the log file clean.

    Examples
    --------
    Run SWarp for a single tile::

        >>> run_back('LMC_c42', 'T07')
        ***BackPrep: Completely done for tile T07 in 123.4 s

    """
    xstart = qstart = start_time = timeit.default_timer()
    run_dir = '%s/%s/%s/' % (SWARPDIR, field, tile)

    try:
        os.chdir(run_dir)
    except:
        print('BackPrep: Error: Could not cd to %s' % run_dir)
        os.chdir(exec_dir)

    command = 'RunSwarp'
    outfile = 'Swarp.out.txt'
    xout = subprocess.run(command, shell=True, capture_output=True)

    z = open(outfile, 'w')
    zz = xout.stderr.decode()
    lines = zz.split('\n')
    xlines = []
    for line in lines:
        if line.count('\x1b[1M>') == 0:
            xlines.append(line)
            z.write('%s\n' % line)
    z.close()

    for xline in xlines[-10:]:
        print(xline)

    current_time = timeit.default_timer()
    print('\n***BackPrep: Completely done for tile %s in %.1f s\n' % (tile, current_time - xstart))

    os.chdir(exec_dir)
    return


def steer(argv):
    """Parse command-line arguments and coordinate BackPrep operations.

    This steering function processes command-line arguments and orchestrates
    the setup and execution of SWarp background processing for specified
    fields and tiles.

    Parameters
    ----------
    argv : list of str
        Command-line arguments including script name.

    Returns
    -------
    None

    Notes
    -----
    The function processes the following workflow:

    1. Parse command-line arguments
    2. Determine which tiles to process
    3. Set up SWarp configuration for each tile
    4. Optionally execute SWarp commands (if -run specified)
    5. Log all operations to a field-specific log file

    Global variables SWARPDIR and PREPDIR may be modified if -indir is used.

    Examples
    --------
    From command line::

        $ python BackPrep.py -run LMC_c42 T07

    From within Python::

        >>> import sys
        >>> steer(['BackPrep.py', '-run', 'LMC_c42', 'T07'])

    See Also
    --------
    prep_one_tile : Set up SWarp for a single tile
    run_back : Execute SWarp processing

    """
    global SWARPDIR
    global PREPDIR

    field = ''
    tiles = []
    xall = False
    xrun = False
    xredo = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-all':
            xall = True
        elif argv[i] == '-run':
            xrun = True
        elif argv[i] == '-setup':
            xrun = False
        elif argv[i] == '-redo':
            xredo = True
        elif argv[i] == '-indir':
            i += 1
            PREPDIR = os.path.abspath(argv[i])
            SWARPDIR = 'XXX/'
        elif argv[i][0] == '-':
            print('BackPrep: Error: Unknown switch %s' % argv[i])
            return
        elif field == '':
            field = argv[i]
        else:
            tiles.append(argv[i])
        i += 1

    if xrun:
        print('BackPrep: Preparing runs and executing')
    else:
        print('BackPrep: Preparing run, but NOT executing')

    if xall:
        # Assume all directories with fits files should be searched
        # Assume we could have both background subtracted and non
        # background subtracted data to deal with
        xfiles = glob('%s/%s/*/*.fits' % (PREPDIR, field))
        xdirs = []
        for one in xfiles:
            words = one.split('/')
            xdirs.append(words[-2].replace('_b', ''))

        tiles = np.unique(xdirs)

    if len(tiles) == 0:
        print('BackPrep: The tiles to be processed must be listed after the field, unless -all is invoked')

    open_log('%s.log' % field, reinitialize=False)
    number = []
    for one in tiles:
        n = prep_one_tile(field, one, xredo)
        log_message('BackPrep: Setup %s %s ' % (field, one))
        number.append(n)

    if xrun:
        for one in tiles:
            log_message('BackPrep: Starting to run %s %s' % (field, one))
            run_back(field, one)
            log_message('BackPrep: Finished %s %s' % (field, one))

    close_log()

    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
