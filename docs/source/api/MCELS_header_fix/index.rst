MCELS_header_fix
================

.. py:module:: MCELS_header_fix

.. autoapi-nested-parse::

   The MCELS images (both for the LMC and SMC) have

   Space Telescope Science Institute

   Synopsis
   --------

   The MCELS images (both for the LMC and SMC) have
   some header info that is depracated, notably:

   Command Line Usage
   ------------------

   ::

       usage: head_fix.py filename

   Description
   -----------

   Primary routines:

       doit

   Primary Routines
   ----------------

   doit

   Notes
   -----

   WARNING - while this does fix the actual problem one will likely still see warnings
       apparently there are problems with the astroy wcs library.  To suppress these warnings
       add to your code the following:

       import warnings
       from astropy.wcs import WCS, FITSFixedWarning
       warnings.filterwarnings('ignore', category=FITSFixedWarning, message=".*datfix.*")

   Version History
   ---------------

   250731 ksl
       Coding begun



Functions
---------

.. autoapisummary::

   MCELS_header_fix.fix_fits_header
   MCELS_header_fix.parse_non_iso_date
   MCELS_header_fix.steer


Module Contents
---------------

.. py:function:: fix_fits_header(filename)

.. py:function:: parse_non_iso_date(date_str)

   Try to convert 'DD/MM/YY' → 'YYYY-MM-DD'


.. py:function:: steer(argv)

   Just steer


