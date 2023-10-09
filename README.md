This repository contains tools developed to reduce a set
of DECam images acquired to map SII and Ha emission in the 
Magellanic Clouds


To use these tools this directory should be placed both
in PATH and PYTHONPATH

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


**The version of astropy that is used is relatively recent**

If you encounter problems, one should create a specific conda
enviromment with something like the following commands:

````
conda create --name kred  scipy
conda install astropy
conda install matplotlib
````

and run in this enviroment which should have versions of these
three packages greater than or equal to:

````
scipy 1.10.1
astropy 5.1
numpy 1.25.0
````

The routines assume a specific directory structure 
to one of our own.  At present the new structure
presumes:

The programs are all run from a toplevel directory.  Underneath
this directory there must be:

* a directory that is (or is symbolically linked) to the
(some or all of) DECam\_MEF directory stucture.  The top level
direcory must named DECam\_MEF

Additional directories will be created as various
programs are run.  

To get help on the options of the various routines, one should
type the program name, with a -h option.  

**Several of the programs have a flag -bsub that change 
the data on which the operate.  The first discussion
of each routine assumes the -bsub flag is NOT invoked
so that users can get a clearer idea of the flow. 
The -bsub option will be explained later**

## Standard Processing

To process the data without attempting to improve over
the background options in MefPrep, the programs should
be run in the following order:

* (MefCheck.py - This should not normally be necessary, but it
will check if there are problems with the headers in the individual
Mef files.  It creates a table called \_mef\_qual.tab which indicates
which Mef files are problematic, but is not used by any of the
downstream programs.)

* MefSum.py - This summarizes the mef files in a field.
If it does not already exist this program creates a directory 
that is created called Summary and writes two astropy
tables that contain the summary information


Note that MefSum.py -all wil attempt to create a summary of all of the files
that one has in the DECam\_MEF directory.   This routine is parallelized;
each thread handles an individual mef file.  (In a previous version,
the routine was parallelized so that each thread handled a different field;
the new version is faster for processing a single field, but not if
you are processing a large number of fields.)
The program does not
check to see if certain Fields have already been analyzed, 
and so if one adds a single Field one should not use the -all option.

Note - In the current version of the program the most important
statistic is produced is the Mode, which is nominally what is 
normall taken for the background.  The mean, median, and 
standard deviation, which are also calculated are not the
sigma clipped versions.  There is an option to create sigma
clipped versions of the median, and mean, but this takes about
a factor of 3 times longer to run.)

* MefPrep.py - This routine rescales all of the images in 
a field or fields, subtracts background if desired, and
stores the results in "data" directories located with the
DECam\_PREP directories.   For example,
if one processes field LMC\_c45, the data will be in
DECam\_CCD/LMC\_c45/data.

Note - this routine has switches that allow one 
to subtract background or not from the original
images.  The default currently is to subtract a single
background from all of the CCDS associated with a given
exposure.  The default is to subtract the median value 
of the mode for the individual detectors.


Also, like MefSum.py one can run this on all of the fields
with the -all switch.


* SetupTile.py - This is intended to identify CCD images
that need to be processed.  It produces 
tables in the Summary directory that identify what 
CCDs need to be processed to produce images of a single 
field. It also cretes directories DECam\_PREP/LMC\_c45/Tile01 etc
that contain links to the appropriate files in the data directories.

Note - This routine uses the config file.  If one wanted to create 
data that was for a different region, then one would need to point 
this to a different configuration file.  Note though, that at
present each field is separate.


* SwarpSetup.py - This creates directories inputs (the.run and .default files)
for running swarp to create the tile images on the tile images.  It assumes we
want to make tile images for each filter and for each 
exposure time.    


Note - At present this routine does not use the field
centers defined in the .config file, but takes the center
of all of the images in a tile to be the center of the
swarp field.  

 When run without the -bsub flag, the routine creates directories
 in DECam\_SWARP and sets up the inputs so the files which will
 be 'swarp'ed come from the subdirctories of DECam\_CCD.


* Swarp.py -  The purpose of Swarp.py 
is to run all of the commands that
have been created with SetupSwarp.  


* SwarpEval.py - This routine creates a set of evaluation
figures for the Swarped images.  The scaling is intended
to highlight flaws.  The images for a field 
are stored in the directory DECam\_SWARP/field/eval
where field is the field name, e.g. LMC\_c42

* CleanStars.py - This routines attempts to create continuum
subtracted images by first creating emission line
fre images of the N708 and r-band images and then 
subtracting these from the emission line images.  The
results are stored in the directory DECam\_SUB/field/tile


**All of the scripts in this sequence should
be run from the toplevel directory for the data
reduction.**


## Using overlaps regions to improve the background matching

The various scripts above allow one to produce swarped imaged
for each tile, but they do not include any attempt to
match backgrounds between images.


In order to improve backrouund matching one
should run the scripts below. They depend on MefPrep having
been run on a field, BUT the other steps in the standard
processing are unnecessary.


The routines should be run in the following order.


* FindOverlaps.py uses information, namely the RAs and DECs of the corners of
the indiviual CCDs, contained in the headers of the images to
identify what images overlap with others.  The routine assumes that the
only overlaps we are concerned with are those between images with the
same filter and exposure time.  

* BackPrep.py creates a directory structure DECam\_BACK parallel to the other
directory structures for storing images created for the purpose estimating
background in the overlap regions of images.  BackPrep also creates Swarp
commands for projecting all of the images (which have gone through the PrepFiles
stage) onto the same WCS.  This WCS has larger pixels than the original WCS
in order to keep the files sizes plausible. When run with the -run option 
BackPrep then Backprep will run Swarp to create these images.

Note -  The files created by BackPrep are large.  If one runs BackPrep on
all of LMC\_c42, for example, about 420 GB of images will be created. 

* BackStats.py determines the background in the overlap regions of images as
identified in FindOverlaps, using the impages created by BackPrep.  The
routine produces a single file for each tile that contains estimates of the
background. 

Note _ Optionally, with the -rm option, BackStats removes the files created by BackPrep,
to retain diskspace, since BackPrep can generate up to a half a TB of
data.

* BackCalc.py uses the differences calculated with BackStats between 
images to determine an optimal set of offsets to add or 
subtract from the images to produce better image matching.  

* BackSub.py uses the results of BackCalc to create new images with 
the backgounds values produced by BalCalc.py subtracted.  The outputs
appear in a directories, with names of  DECam\_Prep2/field/tile.  

* SwarpSetup.py is run (again) but with the -bsub flag in order to 
create directories of the form  DECam\_Swarp/field/tile\_b and commands
in these directories that point to the files cerated by BackSub.py

* Swarp.py  is run again, but with the -bsub flag to run the comamnds
from SetupSwarp.py

* SwarpEval.py  is run to make images (pngs) of for both types of 
processing. 

* CleanStars.py is run with the bsub option.  It now uses data in 
the DECam\_SWARP2 directory, and it writes the data to the 
DECam\_SUB2 directory

The progams generally write information to a log file, one log file
for each field. The file is initialised by running MefSum, and then
each program has the ability to write to it.  

When trying to process a significant amount of data, it is useful 
to carry out the first two steps together, e.g.

````
MEFSum.py  -np 8 LMC_c42
MefPrep.py -np 8 LMC_c42
````

on the field or fields one wants to analyze.

At that point one can create a file to process
either a single tile, or a field completely.  As an example
one could then process all of LMC\_c42, with the following 
set of commands:

````

SetupTile.py -all LMC_c42
SwarpSetup.py -all  LMC_c42
Swarp.py -all LMC_c42
SwarpEval.py -all DECam_SWARP LMC_c42

FindOverlaps.py -all LMC_c42
BackPrep.py -all -run LMC_c42
BackStats.py -all -np 8 -rm LMC_c42 
BackCalc.py -all LMC_c42
BackSub.py -all -np 8 LMC_c42
SwarpSetup.py -all -bsub LMC_c42 
Swarp.py -all -bsub LMC_c42 T02
SwarpEval.py -all DECam_SWARP2 LMC_c42
CleanStars.py -all LMC_c42
SwarpEval.py -all DECam_SUB2 LMC_c42
````

Creating a command file and commenting out the sections 
that one has completed is a sensible approch to precesing.

Note that the first 4 commands here are not normall necessary
and are only usdful if one wants to see the immprovement
from a more sophistcated image mataching of backgrounds


## Additional tools

Two other programs are currently available.  Both make use of
astropy tables with (minimally) the columns: Source\_name,RA,Dec
where RA and Dec must be in degrees.  (Examples can be found
in the config folder, for LMC and SMC SNRs.)

* Cutout.py can be used to created cutout images at the positions of a set of objects.  Depending
on the command line options, one can create cutouts of all the fits images in a particular directory
or in all of the images in that directory or any of its subdirectories.

* SetupSpecial.py can be used to create 'tiles' at an arbitrary set of positions, which can then
be processed in the standard manner.  Depending on the command line switches, the program will
attempt to use only data from specific fields, or all of the data that exists  on one's machine.

