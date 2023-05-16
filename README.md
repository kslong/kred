This repository contains tools developed to reduce a set
of DECam images acquired to map SII and Ha emission in the 
Magellanic Clouds


To use these tools this directory should be placed both
in PATH and PYTHONPATh

If one has placed kred in one's home directory, the following
should be added to ones .bashrc or .bash\_profile

    export KRED=$HOME/kred/
    PATH=$PATH:$KRED/py_progs/
    PYTHONPATH=$PYTHONPATH:$KRED/py_progs/


There is an Jupyter script that performs very basic tests
of the installation in the exmaples directory.

The tools require a set of relatively commmon modules, such
as astropy, and normally are run from within an appropriate 
conda enviromnment.

The routines assume a specfic directory structure THAT IS
CURRENTLY UNDER DEVELOPMENT as we switch from a PHOTPIPE
model to one of our own.  At present the new structure
presumes:

The programs are all run from a toplevel directory:

* a directory that is (or is symbolically linked) to the
(some or all of) DECam\_MEF directory stucture.  The top level
direcory must named DECam\_MEF


At present two programs exist that use the new directory
structure

* MEFSum.py - This summarizes the mef files in a directory 
that is created called Summary.  

* SetupTile.py - This is intended to identify CCD images
that need to be processed.  It produces 
tables in the Summary directory that identify what 
CCDs need to be processed to produce images of a single 
field.

* PrepFiles.py - This creates a location to put the files
that contribute to a single tile and then  rescales the 
individual CCD images that contritute to a tile and 
stores them.   The top level directory for storing
the results of PrepFiles.py is called DECam\_prep.

-----

Portions of the README file beyound this point have not
yet been updated.

-----




The top level directory can be named whatever one wishes,
within this directory one should create two directories


* a directory from which to run all of the programs, e.g. progs.
At present the easiest thing to do, is to copy the .py in py\_progs
and any of the ipynb files in the examples directory there


Running the programs will create several other directories
that parallel that would be creted by Armin Rest's Photopipe routines

* workspace - which has subdirecotires with all of the outputs, and
* workspace/Summary - which will have summaries of the characteristics
of the files that go to make up a tile

At present there are 3 main scripts

* PrepFiles.py  -  This reads from the rawdata and creates background
and rescaled raw data files for the various tiles in a field.  This can
be used to process one or more tiles (or all of them) in a field with 
a single command.

* SumFiles.py -   This reads the results of PrepFiles runs and creates a
set of summary tables, one for each tile in a field.  Again with 
the appropriate inputs one can process one or more tiles, or all of the 
tiles in a field.

The inputs require for either of the above commands can be found by entering

Prepfiles.py -h 

or

SumFiles.py -h


* Swarp.py - This routine is used to create inputs to allow swarp to run
on various images.  It should proably be run from a Jupyter script,
and example of which can be found in the examples directory




