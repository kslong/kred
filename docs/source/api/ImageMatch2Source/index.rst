ImageMatch2Source
=================

.. py:module:: ImageMatch2Source

.. autoapi-nested-parse::

   ImageMatch2Source - Match sources to images by position

   Space Telescope Science Institute

   Synopsis
   --------

   Given a source catalog and an image catalog (from ImageSum.py), find which
   images cover each source position. Returns the closest N images of a specified
   type for each source.

   Command Line Usage
   ------------------

   ::

       ImageMatch2Source.py [-h] [-out outname] [-sep arcmin] [-n_closest N]
                            file_tab source_tab image_type

   Arguments
   ---------

   file_tab
       Image catalog table produced by ImageSum.py. Must contain columns:
       filename, RA, Dec, Image_type.

   source_tab
       Source catalog table with columns: Source_name, RA, Dec (at minimum).

   image_type
       Type of image to match (e.g., 'ha', 's2_sub_r'). Must match values in
       the Image_type column of file_tab.

   **Optional Arguments:**

   -h
       Print this documentation and exit.

   -out outname
       Set output filename. Default: ``XX_{image_type}.{source_tab}``

   -sep arcmin
       Maximum separation between source and image center in arcmin.
       Default: 33 arcmin.

   -n_closest N
       Return up to N closest matches per source. Default: 1.

   Description
   -----------

   This script performs spatial matching between a list of astronomical sources
   and a catalog of available images. For each source, it finds images where the
   source falls within a specified distance of the image center.

   This is useful for:

   * Finding which tile(s) contain a given source
   * Selecting the best image for each source (closest to center)
   * Identifying sources covered by multiple overlapping images

   Workflow
   --------

   1. Run ImageSum.py to produce an image catalog for a directory
   2. Run this script to match sources to available images
   3. Use the output with XSnap.py to create snapshots centered on each source

   Output Columns
   --------------

   Source_name        Source identifier from source_tab
   RA                 Right ascension of the source (degrees)
   Dec                Declination of the source (degrees)
   filename           Path to the matching image file
   separation_arcmin  Distance from source to image center (arcmin)
   rank               Match rank (1 = closest, 2 = second closest, etc.)

   Examples
   --------

   Find the closest image for each source::

       ImageMatch2Source.py Image_Sum_DECam_SUB2.txt sources.txt s2_sub_r

   Find up to 3 closest images within 20 arcmin::

       ImageMatch2Source.py -sep 20 -n_closest 3 Image_Sum_DECam_SUB2.txt sources.txt s2_sub_r

   Specify custom output filename::

       ImageMatch2Source.py -out my_matches.txt Image_Sum_DECam_SUB2.txt sources.txt s2_sub_r

   History
   -------

   240813 ksl  Coding begun



Functions
---------

.. autoapisummary::

   ImageMatch2Source.doit
   ImageMatch2Source.find_closest
   ImageMatch2Source.steer


Module Contents
---------------

.. py:function:: doit(source_tab='smc_snr_cotton24.txt', file_tab='Image_Sum_DECam_SUB2.txt', image_type='s2_sub_r', nmatch=1, max_sep=20.0)

   Match sources to images of a specific type.

   This is the main processing function that reads source and file catalogs,
   filters for the desired image type, and performs the spatial matching.

   Parameters
   ----------
   source_tab : str, optional
       Filename of source catalog table. Must contain columns: Source_name,
       RA, Dec. Default: 'smc_snr_cotton24.txt'.
   file_tab : str, optional
       Filename of file catalog table (created by ImageSum.py). Must contain
       columns: filename, RA, Dec, Image_type. Default: 'Image_Sum_DECam_SUB2.txt'.
   image_type : str, optional
       Type of image to match (e.g., 's2_sub_r', 's2_r', etc.). Must match
       values in the Image_type column of file_tab. Default: 's2_sub_r'.
   nmatch : int, optional
       Number of closest matches to return for each source. Default: 1.
   max_sep : float, optional
       Maximum separation in arcminutes between source and image center.
       Default: 20.0.

   Returns
   -------
   Table or empty list
       Astropy Table containing matched sources and files with columns:
       
       * Source_name : str - Source identifier
       * RA : float - Source RA (formatted to 5 decimals)
       * Dec : float - Source Dec (formatted to 5 decimals)
       * filename : str - Matched image filename
       * separation_arcmin : float - Separation in arcminutes
       * rank : int - Match rank
       
       Returns empty list if errors occur or no matches are found.

   Notes
   -----
   **Error Handling:**

   The function will print error messages and return an empty list if:

   * Either input file cannot be read
   * No images of the specified type exist in file_tab
   * No sources fall within max_sep of any image

   If no images of the specified type are found, the function prints a list
   of available image types.

   **Processing Steps:**

   1. Read file catalog and filter by image_type
   2. Read source catalog
   3. Call find_closest() to perform matching
   4. Join results with original source coordinates
   5. Format output columns

   Examples
   --------
   >>> # Match sources to r-band subtracted images
   >>> matches = doit('sources.txt', 'Image_Sum.txt', 's2_sub_r')
   Of 125 sources in sources.txt, matched 98 in images

   >>> # Find 3 closest matches within 30 arcmin
   >>> matches = doit('sources.txt', 'Image_Sum.txt', 's2_sub_r',
   ...                nmatch=3, max_sep=30.0)
   Of 125 sources in sources.txt, matched 98 in images


.. py:function:: find_closest(source_table, file_table, max_sep_arcmin, n)

   Find the closest N files to each source within a maximum separation.

   This function performs a spatial match between sources and image files,
   returning the N closest matches for each source that fall within the
   specified separation limit.

   Parameters
   ----------
   source_table : Table
       Astropy Table containing sources with columns:
       
       * Source_name : str - Source identifier
       * RA : float - Right ascension in degrees
       * Dec : float - Declination in degrees
       
   file_table : Table
       Astropy Table containing files with columns:
       
       * filename : str - Image filename
       * RA : float - Right ascension of image center in degrees
       * Dec : float - Declination of image center in degrees
       
   max_sep_arcmin : float
       Maximum separation allowed between source and image center in
       arcminutes. Sources beyond this distance will not be matched.
   n : int
       Maximum number of files to return for each source. If fewer than
       N files are within max_sep, all qualifying files are returned.

   Returns
   -------
   Table
       Astropy Table with columns:
       
       * Source_name : str - Source identifier
       * filename : str - Matched image filename
       * separation_arcmin : float - Separation in arcminutes (formatted to 3 decimals)
       * rank : int - Match rank (1=closest, 2=second closest, etc.)
       
       If no matches are found, returns an empty table.

   Notes
   -----
   **Algorithm:**

   1. Convert RA/Dec coordinates to SkyCoord objects
   2. For each source, compute angular separations to all files
   3. Filter files within max_sep
   4. Sort by separation (closest first)
   5. Select up to N closest matches
   6. Assign ranks based on separation

   **Performance:**

   The function uses vectorized astropy operations for efficient coordinate
   matching. For large catalogs (>10000 sources or files), consider chunking
   the input tables.

   Examples
   --------
   >>> from astropy.table import Table
   >>> sources = Table({'Source_name': ['SN1', 'SN2'], 
   ...                  'RA': [10.0, 20.0], 
   ...                  'Dec': [-30.0, -40.0]})
   >>> files = Table({'filename': ['img1.fits', 'img2.fits'],
   ...                'RA': [10.1, 20.0],
   ...                'Dec': [-30.1, -40.0]})
   >>> matches = find_closest(sources, files, max_sep_arcmin=30.0, n=1)
   >>> print(matches)
   Source_name  filename    separation_arcmin  rank
   ----------- ----------  -----------------  ----
           SN1  img1.fits             8.485     1
           SN2  img2.fits             0.000     1


.. py:function:: steer(argv)

   Parse command-line arguments and execute source-to-image matching.

   This is the main entry point for command-line execution. It parses
   arguments and calls doit() to perform the matching.

   Parameters
   ----------
   argv : list
       Command-line argument list (typically sys.argv).

   Returns
   -------
   None
       Results are written to a file. The function prints status messages
       and errors to stdout.

   Command-Line Format
   -------------------
   ::

       ImageMatch2Source.py [-h] [-out name] [-sep arcmin] [-n_closest N]
                            file_tab source_tab image_type

   See module docstring for detailed argument descriptions.

   Output Files
   ------------
   **Default naming:** ``XX_{image_type}.{source_tab}``

   **Custom naming:** Specified with -out option

   **Format:** ASCII fixed-width two-line format (Astropy standard)

   Examples
   --------
   Basic usage from command line::

       $ python ImageMatch2Source.py files.txt sources.txt s2_sub_r

   Advanced usage::

       $ python ImageMatch2Source.py -sep 25 -n_closest 2 -out matches.txt 
               files.txt sources.txt s2_sub_r

   Notes
   -----
   **Error Handling:**

   * Prints help if -h is specified
   * Validates command-line parsing
   * Checks for required arguments
   * Reports unknown switches

   If no matches are found or errors occur during processing, appropriate
   error messages are printed and the function returns without creating
   an output file.


