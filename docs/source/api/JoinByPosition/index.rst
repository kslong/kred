JoinByPosition
==============

.. py:module:: JoinByPosition

.. autoapi-nested-parse::

   JoinByPosition - Position-based table cross-matching

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



Functions
---------

.. autoapisummary::

   JoinByPosition.do_join
   JoinByPosition.do_one
   JoinByPosition.read_table
   JoinByPosition.steer


Module Contents
---------------

.. py:function:: do_join(table1, table2, max_sep=1.0, keep_all=False, ra1='RA', dec1='Dec', ra2='RA', dec2='Dec')

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


.. py:function:: do_one(table1_path, table2_path, max_sep=1.0, keep_all=False, outfile='', ra1='RA', dec1='Dec', ra2='RA', dec2='Dec')

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


.. py:function:: read_table(filename)

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


.. py:function:: steer(argv)

   Parse command line arguments and execute cross-match.

   Parameters
   ----------
   argv : list
       Command line arguments (sys.argv)


