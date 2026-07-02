SetupSpecial
============

.. py:module:: SetupSpecial

.. autoapi-nested-parse::

   Setup the directories and identify the files that are required to carry

   Space Telescope Science Institute

   Synopsis
   --------

   Setup the directories and identify the files that are required to carry
   out the initial processing steps to run PrepFiles on a set of data

   Command Line Usage
   ------------------

   ::

       usage: SetupSpecial  -h [-all]  [-field LMC_c40] [-field LMC_c42] source_table.txt  [source_name1] [Source_name]

       where:

       source_table.txt is an astropy table with the source names and positions

       and the optional prameters

       -field is that name of a field, e.g LMC_c42.  For multiple files one can add
       tile1, tiles are the names of tiles as specified in the MC_tiles.txt file
       If no fields are provided all of the data will be searched.

       -h prints this help
       -all  causes the program to create 'tiles'  containing the relevant
       infomation for all of the sources in the table

   Description
   -----------

   This program reads an astropy table containing at least the following columns

       Source_name  - a string (one_word no space) naming the object
       RA, Dec       - the right ascension and declination of the object in degrees

       and constructs 'tiles' of one or more of the objects from some of all of the
       available data.

       The astropy table can be local, or it can be stored in the config directory

       If the table also contains a column

       Size   - a size in degrees

       then the that size will be used in selecting the images that can contribute
       to the total.  If it does not then a default is assumed.

       At present the program does not use the size for anything except what images
       to include in the processing of that tile.

       The program works by reading the det.tab files created by MefSum.py

   Primary Routines
   ----------------

   doit

   Notes
   -----

   History:

   230513 ksl Coding begun

   Version History
   ---------------

   230513 ksl
       Coding begun



Attributes
----------

.. autoapisummary::

   SetupSpecial.CWD
   SetupSpecial.DATADIR
   SetupSpecial.PREPDIR


Functions
---------

.. autoapisummary::

   SetupSpecial.find_files_to_prep
   SetupSpecial.setup_objects
   SetupSpecial.setup_one_tile
   SetupSpecial.setup_tiles
   SetupSpecial.steer


Module Contents
---------------

.. py:data:: CWD

.. py:data:: DATADIR

.. py:data:: PREPDIR

.. py:function:: find_files_to_prep(field='LMC_c45', ra_center=81.1, dec_center=-66.17, size_deg=0.1)

.. py:function:: setup_objects(source_names, table_name, fields, xdir, process_table)

       


.. py:function:: setup_one_tile(field=['LMC_c42'], xmaster='snr', tile='foo', ra=81.108313, dec=-66.17728, size_deg=0.67, process_table='')

.. py:function:: setup_tiles(xtab)

   Run setup_one_tile multiple times, with the possibility of parallel
   processing


.. py:function:: steer(argv)

   This is just a steering routine

   The command line it interprets is

   SetupSpecial  [-all]  [-field LMC_c40] [-field LMC_c42] source_table source_name1 [Source_name]

   where -all says to do all of the sources in the xtab file


