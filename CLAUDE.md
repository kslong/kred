# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KRED is a Python-based DECam data reduction pipeline for imaging the Magellanic Clouds in emission lines (Ha and [SII]). It processes multi-extension FITS (MEF) files through a sequence of scripts to produce calibrated, background-matched tile mosaics.

## Environment Setup

```bash
export KRED=$HOME/kred/
export PATH=$PATH:$KRED/py_progs/
export PYTHONPATH=$PYTHONPATH:$KRED/py_progs/

# Conda environment
conda create --name kred scipy
conda install astropy matplotlib photutils astroquery gaiaxpy
```

Required versions: scipy >= 1.10.1, astropy >= 5.1, numpy >= 1.25.0

## Running Scripts

All scripts are in `py_progs/` and run from the command line with `-h` for help:

```bash
MefSum.py -h                    # Show usage
MefSum.py -np 8 LMC_c42         # Process field with 8 parallel processes
MefSum.py -all                  # Process all fields
```

Common flags across scripts:
- `-h` - Show help (prints module docstring)
- `-np N` - Number of parallel processes (default varies)
- `-all` - Process all fields/tiles
- `-bsub` - Use background-subtracted data (SwarpSetup, Swarp, CleanStars)

## Standard Processing Pipeline

Run from the working directory containing `DECam_MEF/`:

```bash
# 1. Validate and summarize data
MefCheck.py LMC_c42
MefSum.py -np 8 LMC_c42

# 2. Prepare images
MefPrep.py -np 8 LMC_c42

# 3. Setup and create tile mosaics
SetupTile.py -all LMC_c42
SwarpSetup.py -all LMC_c42
Swarp.py -all LMC_c42

# 4. Background matching (optional but recommended)
FindOverlaps.py -all LMC_c42
BackPrep.py -all -run LMC_c42
BackStats.py -all -np 8 -rm LMC_c42
BackCalc.py -all LMC_c42
BackSub.py -all -np 8 LMC_c42
SwarpSetup.py -all -bsub LMC_c42
Swarp.py -all -bsub LMC_c42

# 5. Continuum subtraction
CleanStars.py -all LMC_c42
```

## Photometric Consistency Workflow (committed 2026-07-01)

Run `CheckPhot.py LMC_c42` to execute the full intra- and inter-filter
consistency check.  Key output in `Summary/{field}_phot_check_intra.fits`
and `Summary/{field}_phot_check_inter.fits`.  See
`docs/source/relative_photometry.rst` for full interpretation guide.

## Architecture

### Three-Layer Design

All major scripts follow this pattern:

1. **CLI Layer** (`steer(argv)`) - Manual `sys.argv` parsing, dispatches to processing
2. **Orchestration Layer** (`do_one()`, `do_many()`) - File I/O, error handling, parallelization
3. **Core Processing** (e.g., `do_forced_photometry()`) - Pure computation, returns Table or error string

Example pattern:
```python
def steer(argv):
    # Parse args manually (not argparse)
    if len(filenames) == 1 or np_proc < 2:
        for file in filenames:
            do_one(file)
    else:
        do_many(filenames, n_processes=np_proc)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
```

### Key Modules

| Module | Purpose |
|--------|---------|
| `log.py` | Timestamped logging with git commit tracking |
| `ImageSum.py` | FITS header utilities, WCS extraction |
| `GaiaCat.py` | GAIA DR3 catalog retrieval (`get_gaia()`) |
| `Smash.py` | SMASH DR2 catalog via NOAO Data Lab |

### Output Directory Structure

```
Working Directory/
├── Summary/           # MefSum output - field/tile inventory
├── DECam_PREP/        # MefPrep - rescaled images
├── DECam_SWARP/       # Swarp - tile mosaics (no background matching)
├── DECam_BACK/        # BackPrep - temporary background images (can be deleted)
├── DECam_PREP2/       # BackSub - background-corrected images
├── DECam_SWARP2/      # Swarp -bsub - tile mosaics (with background matching)
├── DECam_SUB/         # CleanStars - continuum-subtracted images
├── TabPhot/           # Photometry tables (MefPhot, PhotCompare)
├── Figs_phot/         # Photometry figures
└── GAIA/              # Cached GAIA catalogs
```

## Configuration Files (config/)

- `MC_tiles.txt` - Tile definitions (Field, Tile, RA, Dec, Size)
- `DeMCELS_images.txt` - Filter/exposure combinations for stacking
- `lmc_snr.txt`, `smc_snr.txt` - Object catalogs (Source_name, RA, Dec)

## Photometry Tools

```bash
# Forced aperture photometry
MefPhot.py image.fits -r 6 -b_in 8 -b_out 12

# Compare to GAIA magnitudes
PhotCompare.py image.fits

# PSF construction
StarFind.py image.fits
PsfBuild.py image.fits TabPhot/all_stars.fits
```

## External Dependencies

- **SWarp** - External command-line tool for image mosaicing (must be installed separately)
- **GAIA archive** - Online catalog access via astroquery
- **NOAO Data Lab** - For SMASH catalog queries (Smash.py)
