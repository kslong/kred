master2reg
==========

.. py:module:: master2reg

.. autoapi-nested-parse::

   master2reg - Convert Master Table to DS9 Region File

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



Attributes
----------

.. autoapisummary::

   master2reg.argc


Functions
---------

.. autoapisummary::

   master2reg.doit
   master2reg.radec2deg
   master2reg.read_masterfile
   master2reg.write_regionfile


Module Contents
---------------

.. py:data:: argc

.. py:function:: doit(argv)

   Parse command line and execute master to region file conversion.

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


.. py:function:: radec2deg(ra, dec)

   Convert RA/Dec strings to decimal degrees.

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


.. py:function:: read_masterfile(filename, xtype='circle', xmajor=3, xminor=3, xtheta=0.0)

   Read a master table file and standardize column names.

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


.. py:function:: write_regionfile(regionfile, records, source='unknown', type='unknown', color='red')

   Write a DS9 region file from a master table.

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


