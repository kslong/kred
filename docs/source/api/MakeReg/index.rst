MakeReg
=======

.. py:module:: MakeReg

.. autoapi-nested-parse::

   Create region files for fits files in a directory,

   Space Telescope Science Institute

   Synopsis
   --------

   Create region files for fits files in a directory,
   assuming that the images have the corners of the images
   store in the headers.

   Command Line Usage
   ------------------

   ::

       usage: MakeReg.py xdir

       where xdir is a directory that contains fits files with
       the corners of the images given

   Description
   -----------

   The routine produces a set of region files, one
       for each unique filter and exposure time for
       the files that are in the directory

   Primary Routines
   ----------------

   gen_regions

   Notes
   -----

   History:

   231017 ksl Coding begun

   Version History
   ---------------

   231017 ksl
       Coding begun



Attributes
----------

.. autoapisummary::

   MakeReg.reg_header


Functions
---------

.. autoapisummary::

   MakeReg.gen_regions
   MakeReg.steer


Module Contents
---------------

.. py:function:: gen_regions(xdir='.', out_root='test')

   Locate the fits files in a directory and
   create region files for each unique
   filter and exposure time.


.. py:data:: reg_header
   :value: Multiline-String

   .. raw:: html

      <details><summary>Show Value</summary>

   .. code-block:: python

      """# Region file format: DS9 version 4.1
      global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
      fk5
      """

   .. raw:: html

      </details>



.. py:function:: steer(argv)

   Setup the inputs for this routine


