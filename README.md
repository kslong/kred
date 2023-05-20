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

* PrepSum.py - This summarizes the files that have been
Prepped prior to combining all of the images in 
the various filters with swarp

* Swarp.py - This contains code to setup the swarp
directiories in a standard way and to create .run
.default files and input file lists with a routine
create\_swarp\_command, that is currently best run
from a Jupyter notebook.  One can run the Swarp
commands from inside the notebook or from the command 
line.  As is the case for the rest of the scripts
in this sequence, the notebook and the command
line execution should be run from whatevery one
is using as the toplevel directory for the data
reduction.


The CheckInstall notebook in examples has been updated
to reflect the changes to the various routines, and
to show one how to run through the various scripts
to produce a swapred Ha image for LMC\_c42 tile T08

