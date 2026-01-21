============
Installation
============

This page describes how to install and configure kred for reducing DECam data.

Requirements
------------

Kred requires Python 3.8 or later with the following package versions:

* scipy >= 1.10.1
* astropy >= 5.1
* numpy >= 1.25.0
* matplotlib
* photutils
* astroquery
* gaiaxpy

Additionally, the external program **SWarp** must be installed separately for
image mosaicing.

Installation
------------

Clone the repository
^^^^^^^^^^^^^^^^^^^^

Clone the kred repository to your home directory::

    cd ~
    git clone https://github.com/kslong/kred.git

Create a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

Create and activate a new conda environment with the required packages::

    conda create --name kred python=3.10 scipy astropy numpy matplotlib
    conda activate kred
    conda install -c conda-forge photutils astroquery gaiaxpy

Alternatively, install with pip after creating the environment::

    conda create --name kred python=3.10
    conda activate kred
    pip install scipy astropy numpy matplotlib photutils astroquery gaiaxpy

Set environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Add the following to your shell configuration file (e.g., ``~/.bashrc`` or
``~/.zshrc``)::

    export KRED=$HOME/kred/
    export PATH=$PATH:$KRED/py_progs/
    export PYTHONPATH=$PYTHONPATH:$KRED/py_progs/

Then reload your shell configuration::

    source ~/.bashrc   # or ~/.zshrc

Install SWarp
^^^^^^^^^^^^^

SWarp is required for creating tile mosaics. Install it using your system
package manager:

**macOS (Homebrew)**::

    brew install swarp

**Ubuntu/Debian**::

    sudo apt-get install swarp

**Conda**::

    conda install -c conda-forge astromatic-swarp

Verify installation
^^^^^^^^^^^^^^^^^^^

Verify that kred is correctly installed::

    # Check that scripts are accessible
    MefSum.py -h

    # Check Python imports
    python -c "import MefSum; print('OK')"

    # Check SWarp
    swarp -v

Directory structure
-------------------

Kred expects a specific directory structure for processing. Create a working
directory and ensure it contains::

    Working Directory/
    ├── DECam_MEF/          # Raw multi-extension FITS files
    └── config/             # Configuration files (optional)

Output directories are created automatically during processing:

==================  ==========================================
Directory           Contents
==================  ==========================================
``Summary/``        Field and tile inventory (MefSum)
``DECam_PREP/``     Rescaled images (MefPrep)
``DECam_SWARP/``    Tile mosaics without background matching
``DECam_BACK/``     Temporary background images
``DECam_PREP2/``    Background-corrected images (BackSub)
``DECam_SWARP2/``   Tile mosaics with background matching
``DECam_SUB/``      Continuum-subtracted images (CleanStars)
``TabPhot/``        Photometry tables
``Figs_phot/``      Photometry figures
``GAIA/``           Cached GAIA catalogs
==================  ==========================================

Getting started
---------------

Once installed, see the main documentation for the standard processing pipeline.
A quick test to verify everything works::

    cd /path/to/working/directory
    MefCheck.py LMC_c42
    MefSum.py -np 4 LMC_c42
