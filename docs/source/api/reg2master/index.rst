reg2master
==========

.. py:module:: reg2master

.. autoapi-nested-parse::

   reg2master - Convert DS9 Region File to Master Table

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



Functions
---------

.. autoapisummary::

   reg2master.doit
   reg2master.radec2deg
   reg2master.read_regions
   reg2master.size2arcsec
   reg2master.steer
   reg2master.write_masterfile


Module Contents
---------------

.. py:function:: doit(regionfile, masterfile='')

   Create a master table from a region file.

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


.. py:function:: radec2deg(ra, dec)

   Convert RA/Dec strings to decimal degrees.

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


.. py:function:: read_regions(filename)

   Read and parse a DS9 region file.

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


.. py:function:: size2arcsec(word)

   Convert a size string to arcseconds.

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


.. py:function:: steer(argv)

   Parse command line arguments and execute conversion.

   Parameters
   ----------
   argv : list of str
       Command line arguments (typically sys.argv).

   Returns
   -------
   None
       Calls doit() to perform the conversion.


.. py:function:: write_masterfile(masterfile, records, source='unknown', type='unknown')

   Write a master table file from parsed region records.

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


