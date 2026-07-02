SetupTile
=========

.. py:module:: SetupTile

.. autoapi-nested-parse::

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



Attributes
----------

.. autoapisummary::

   SetupTile.CWD
   SetupTile.DATADIR
   SetupTile.PREPDIR


Functions
---------

.. autoapisummary::

   SetupTile.find_overlapping_images
   SetupTile.get_tile_files
   SetupTile.locate_prepped_files
   SetupTile.populate_tile_dir
   SetupTile.setup_tiles
   SetupTile.steer


Module Contents
---------------

.. py:data:: CWD

.. py:data:: DATADIR

.. py:data:: PREPDIR

.. py:function:: find_overlapping_images(field='LMC_c45', tile='T07', ra_center=81.1, dec_center=-66.17, size_deg=0.1)

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


.. py:function:: get_tile_files(field='LMC_c42', tile='T07', ra=81.108313, dec=-66.17728, size_deg=0.67, s7=True, seeing_max=1000.0, use_all_data=False)

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


.. py:function:: locate_prepped_files(xtab, field, tile)

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


.. py:function:: populate_tile_dir(xtab, field, tile)

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


.. py:function:: setup_tiles(ztab, s7, seeing_max, use_all_data, process_tab)

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


.. py:function:: steer(argv)

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


