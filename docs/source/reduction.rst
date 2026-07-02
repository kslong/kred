=========
Reduction
=========

This page describes the standard data reduction process for DECam images using
kred. All scripts should be run from the top-level working directory containing
the ``DECam_MEF/`` folder.

Overview
--------

The reduction pipeline processes multi-extension FITS (MEF) files through a
sequence of steps to produce calibrated, background-matched tile mosaics. There
are two main workflows:

#. **Standard processing**: Creates tile mosaics without background matching
#. **Background matching**: Improves background uniformity across overlapping images

Both workflows use the same initial preparation steps and share many scripts.

.. note::

   Several programs have a ``-bsub`` flag that changes which data they operate
   on. The descriptions below assume ``-bsub`` is NOT used unless explicitly
   stated. Background-subtracted processing is covered in a separate section.


Standard Processing
-------------------

To process data without background matching, run the following programs in order:

MefCheck.py
^^^^^^^^^^^

Run whenever files have been added to the data repository. Performs basic checks
on files to ensure they can be read and have required header keywords. Also
checks for duplicate files that could cause issues downstream.

::

    MefCheck.py LMC_c42

MefSum.py
^^^^^^^^^

Summarizes the MEF files in a field. Creates a ``Summary/`` directory and writes
two astropy tables containing summary information. Files missing required keywords
(e.g., MAGZERO) are placed in a separate summary file and excluded from further
processing.

::

    MefSum.py -np 8 LMC_c42    # Process single field with 8 processes
    MefSum.py -all             # Process all fields

The routine also checks that filter/exposure time combinations appear in the
instrument configuration file (``config/DeMCELS_images.txt``).

MefPrep.py
^^^^^^^^^^

Rescales all images in a field, optionally subtracts background, and stores
results in ``DECam_PREP/``. For example, processing ``LMC_c45`` creates data in
``DECam_PREP/LMC_c45/data/``.

::

    MefPrep.py -np 8 LMC_c42

By default, subtracts the median value of the mode across all CCDs from each
exposure.  The flux scaling uses the ``MAGZERO`` keyword from each MEF header
as supplied by the NOIRLAB community pipeline, placing all images on a scale
where 1 DN corresponds to magnitude 28.  This is the standard approach and is
sufficient in most cases.

If you need to verify that the ``MAGZERO`` values are mutually consistent
across exposures, or to replace them with empirically-derived zero points, see
:doc:`relative_photometry`.

.. important::

   If you plan to use the ``-use_all_data`` option in SetupTile (to include all
   possible images in a tile), you must run MefPrep on all data first. Otherwise,
   only images processed by MefPrep will be included.

SetupTile.py
^^^^^^^^^^^^

Identifies CCD images that need to be processed for each tile. Creates tables in
``Summary/`` identifying which CCDs contribute to each tile, and creates
directories like ``DECam_PREP/LMC_c45/Tile01/`` containing links to the
appropriate data files.

::

    SetupTile.py -all LMC_c42

This routine uses two configuration files from the ``config/`` directory:

**MC_tiles.txt**
    Defines tile centers and geometries. Modify this to create tiles for
    different regions.

**DeMCELS_images.txt**
    Defines which filter/exposure combinations to include in final images.
    Multiple exposure times can be combined for the same filter.

The ``-use_all_data`` option examines all CCD images in ``DECam_PREP/`` to find
any that overlap the tile area, not just those from the field's exposure sequence.

SwarpSetup.py
^^^^^^^^^^^^^

Creates directories and input files (``.run`` and ``.default`` files) for running
SWarp to create tile mosaics. Creates one tile image per filter and exposure time.

::

    SwarpSetup.py -all LMC_c42

Output is written to ``DECam_SWARP/``.

Swarp.py
^^^^^^^^

Runs all SWarp commands created by SwarpSetup.

::

    Swarp.py -all LMC_c42

SwarpEval.py
^^^^^^^^^^^^

Creates evaluation figures for the swarped images. Scaling is designed to
highlight flaws. Images are stored in ``DECam_SWARP/field/eval/``.

::

    SwarpEval.py -all DECam_SWARP LMC_c42

CleanStars.py
^^^^^^^^^^^^^

Creates continuum-subtracted images by generating emission-line-free images from
N708 and r-band data, then subtracting these from the emission line images.
Results are stored in ``DECam_SUB/field/tile/``.

::

    CleanStars.py -all LMC_c42


Background Matching
-------------------

The standard processing does not match backgrounds between overlapping images.
For improved background uniformity, run the following additional scripts after
MefPrep (the other standard processing steps are not required first).

FindOverlaps.py
^^^^^^^^^^^^^^^

Uses header information (corner coordinates) to identify which images overlap.
Only considers overlaps between images with the same filter and exposure time.

::

    FindOverlaps.py -all LMC_c42

BackPrep.py
^^^^^^^^^^^

Creates a ``DECam_BACK/`` directory structure for storing background estimation
images. Projects all prepared images onto a common WCS with larger pixels to
keep file sizes manageable. With the ``-run`` option, executes SWarp to create
these images.

::

    BackPrep.py -all -run LMC_c42

.. warning::

   The files created by BackPrep are large. Processing all of ``LMC_c42`` creates
   approximately 420 GB of images.

BackStats.py
^^^^^^^^^^^^

Determines the background in overlap regions using images from BackPrep. Produces
a single file per tile containing background estimates.

::

    BackStats.py -all -np 8 -rm LMC_c42

The ``-rm`` option removes BackPrep files after processing to recover disk space.

BackCalc.py
^^^^^^^^^^^

Uses the differences calculated by BackStats to determine optimal offsets for
each image to produce better background matching.

::

    BackCalc.py -all LMC_c42

BackSub.py
^^^^^^^^^^

Applies the offsets from BackCalc to create new images with corrected backgrounds.
Output appears in ``DECam_PREP2/field/tile/``.

::

    BackSub.py -all -np 8 LMC_c42

Re-running SwarpSetup and Swarp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After background subtraction, re-run SwarpSetup and Swarp with the ``-bsub`` flag
to create mosaics from the background-corrected images::

    SwarpSetup.py -all -bsub LMC_c42
    Swarp.py -all -bsub LMC_c42

This creates output in ``DECam_SWARP2/``.

Finally, run CleanStars with the ``-bsub`` flag to create continuum-subtracted
images from the background-matched data::

    CleanStars.py -all -bsub LMC_c42

Output is written to ``DECam_SUB2/``.


Complete Example
----------------

The following command sequence processes a field with full background matching:

.. code-block:: bash

    # Initial preparation
    MefSum.py -np 8 LMC_c42
    MefPrep.py -np 8 LMC_c42

    # Setup tiles
    SetupTile.py -all LMC_c42
    SwarpSetup.py -all LMC_c42

    # Optional: create mosaics without background matching
    # Swarp.py -all LMC_c42
    # SwarpEval.py -all DECam_SWARP LMC_c42

    # Background matching
    FindOverlaps.py -all LMC_c42
    BackPrep.py -all -run LMC_c42
    BackStats.py -all -np 8 -rm LMC_c42
    BackCalc.py -all LMC_c42
    BackSub.py -all -np 8 LMC_c42

    # Create background-matched mosaics
    SwarpSetup.py -all -bsub LMC_c42
    Swarp.py -all -bsub LMC_c42
    SwarpEval.py -all DECam_SWARP2 LMC_c42

    # Continuum subtraction
    CleanStars.py -all -bsub LMC_c42
    SwarpEval.py -all DECam_SUB2 LMC_c42

.. tip::

   Save these commands in a script file and comment out completed sections as
   you progress through the reduction.
