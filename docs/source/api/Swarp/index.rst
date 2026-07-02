Swarp
=====

.. py:module:: Swarp

.. autoapi-nested-parse::

   Create Combined images of the different filters in the field using swarp

   Space Telescope Science Institute

   Command Line Usage
   ------------------

   ::

       Usage:   Swarp.py [-all] [-bsub] field [tiles]

       where -all will cause swarp to be run on all 16 tiles.

       and   -bsub will cause swarp to be run in the DECAM_SWARP2 directories, and should
       be operating on files that are in background subtracted DECAM_PREP2 directories

       If one wants to run only 1 or a few tiles then the command will be something like

       Swarp.py LMC_c42  T01 T03 T05

   Notes
   -----

   It would be sensible to combine SetupSwarp and Swarp



Attributes
----------

.. autoapisummary::

   Swarp.PREPDIR
   Swarp.SWARPDIR
   Swarp.exec_dir


Functions
---------

.. autoapisummary::

   Swarp.add_header_keywords
   Swarp.run_swarp
   Swarp.steer


Module Contents
---------------

.. py:data:: PREPDIR
   :value: b'.'


.. py:data:: SWARPDIR
   :value: 'DECam_SWARP'


.. py:function:: add_header_keywords(root, run_dir='./', field='', tile='')

   Add header keywords to an the results of Swarp


.. py:data:: exec_dir

.. py:function:: run_swarp(field='LMC_c42', tile='T07', bsub=False)

   run_all of the swarp commands in a particular field and tile

   The normal outputs from swarp are writeen to a .txt file,
   a few are written to the screen here so that one can see 
   that the routines have run correctly.  


.. py:function:: steer(argv)

   This is just a steering routine for running swarp on one or more
   tiles from the command line



