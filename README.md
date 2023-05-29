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
that contribute to a single tile and then  writes the 
individual CCD images that contritute to a tile and 
stores them.   The top level directory for storing
the results of PrepFiles.py is called DECam\_prep.
The images are scaled to the same magnitude scale. 
There is also a switch which allows one to attempt
to subtract a background from the images, but in the absence
of this there is no background subtraction

* PrepSum.py - This summarizes the files that have been
Prepped prior to combining all of the images in 
the various filters with swarp

* SetupSwarp.py - This creates inputs (the.run and .default files)
for running swarp
to create the tile images.  It assumes we
want to make tile images for each filter and for each 
exposure time.  It largely replaces the portion of Swarp.py
below that requires one to create Jupyter notebooks

(Note that at present this routine does not use the field
centers defined in the .config file, but takes the center
of all of the images in a tile to be the center of the
swarp field.  This is something that still needs updating)

* Swarp.py -  The primary purpose of Swarp.py when
run from the command line is to run all of the commands that
have been crated with SetupSwarp.  

The script also contains
code to create special 
This contains code to setup the swarp
directories in a standard way and to create .run
.default files and input file lists with a routine
create\_swarp\_command, that is currently best run
from a Jupyter notebook.  One can also run the Swarp
commands from inside the notebook or from the command 
line.   

As is the case for the rest of the scripts
in this sequence, the notebook and the command
line execution should be run from whatevery one
is using as the toplevel directory for the data
reduction.

The various scripts abave allows on to produce swarped imaged
for each tile, but they do not include any attempt to
match backgrounds between images.

Work on this is still underway, but some preperatory scripts
exist (as of 230529).

* FindOverlaps.py uses information, namely the RAs and DECs of the corners of
the indiviual CCDs, contained in the headers of the images to
identify what images overlap with others.  The routine assumes that the
only overlaps we are concerned with are those between images with the
same filter and exposure time.  

* PrepBack.py creates a directory structure DECam\_BACK parallel to the other
directory structures for storing images created for the purpose estimaging
background in the overlap regions of images.  PrepBack also creates Swarp
commands for projecting all of the images (which have gone through the PrepFiles
stage) onto the same WCS.  This WCS has larger pixels than the original WCS
in order to keep the files sizes plausible.

At present, one must actually run the swarp commands from inside the directories,
somehtin that needs to be modificed.

* GetStats.py determines the background in the overlap regions of images as
identified in FindOverlaps, using the impages created by PrepBack.  The
routine produces a single file for each tile that contains estimates of the
background.

At present there is no completed routine to estimate the best offsets
to produce composite images that are smooth across image boundaries.  That is, in NASA
jargon, 'in work'




The CheckInstall notebook in examples has been updated
to reflect the changes to the various routines, and
to show one how to run through the various scripts
to produce a swapred Ha image for LMC\_c42 tile T08

