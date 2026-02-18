===============================
Snapshots - How to create them
===============================

Kred contains several tools for creating snapshots/subimages of objects in situations
where one knows in which image the object lies, or as in the case of the DECam images
objects could be in one of many images, and one needs to identify which FITS file
contains the object of interest.

For a single image such as MCELS one typically needs only to run a single program,
namely :doc:`XSnap.py <api/XSnap/index>`.

XSnap is designed to work with a "master table" that contains, minimally, the RA and Dec
for objects of interest. It produces a FITS file of the snapshot in the ``xdata/``
directory and a figure in the ``ximage/`` directory.

Like nearly all routines in the repository, the options can be displayed with::

    XSnap.py -h

The documentation in XSnap is fairly complete.

Single-image workflow
---------------------

A typical example of a command line for a single file would be::

    XSnap.py -size 10 -type ha -min -1 -max 20 DECam_SWARP/LMC_c42_T01.ha.fits lmc_snr.txt

This will produce 10 arcmin sized snapshots of the LMC SNRs in this particular image,
with plots scaled between -1 and 20. The output files will have the type "ha" in the
filenames so that one can write files from different filters to the same output directory.


Multi-image workflow
--------------------

The second situation for which XSnap is intended is when there are many images and one
wants to select a specific image from which to create each snapshot.

In this case there are 3 steps:

#. Summarize the images with :doc:`ImageSum.py <api/ImageSum/index>`
#. Select the best image for the snapshot of each object with
   :doc:`ImageMatch2Source.py <api/ImageMatch2Source/index>`
#. Make the snapshots with :doc:`XSnap.py <api/XSnap/index>`

Example::

    ImageSum.py DECam_SUB2
    ImageMatch2Source.py Image_Sum_DECam_SUB2.txt smc_snr_cotton24.txt s2_sub_r
    XSnap.py -size 10 -type s2 XX_s2_sub_r.smc_snr_cotton24.txt smc_snr_cotton24.txt

These tools are largely designed for our kred reductions of the DECam data:

**ImageSum.py**
    Searches all of the FITS files in the ``DECam_SUB2/`` directory and any of its
    subfolders, and examines the files to find enough information to locate sources
    within them. It produces an output file, in this case ``Image_Sum_DECam_SUB2.txt``,
    that contains this information.

**ImageMatch2Source.py**
    Reads this file along with a master file containing all of the SNRs in the SMC,
    and locates the best image of type ``s2_sub_r`` from which to create a snapshot.
    This routine assumes that all files with ``s2_sub_r`` are effectively of the same
    group. It produces the file ``XX_s2_sub_r.smc_snr_cotton24.txt`` containing the
    results.

**XSnap.py**
    Reads ``XX_s2_sub_r.smc_snr_cotton24.txt``, which contains information about the
    best image for each source, and joins this with the table containing the object
    positions to produce the snapshots.
