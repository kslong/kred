#!/usr/bin/env python 

"""
Space Telescope Science Institute PSF Photometry Module
========================================================

This module performs PSF photometry on astronomical images using model PSFs,
optimized for bright-star flux accuracy and crowded field performance.

Command Line Usage
------------------
::

    Usage: PsfPhot.py [-out whatever] [-np N] image.fits psf.fits stars.fits
    
    where:
        -out : Output filename root (default: 'psf_phot')
        -np  : Number of parallel processes (default: 1, 0=auto-detect)

Notes
-----
This parallelized version handles large images efficiently and includes:
    
* Saturation masking to prevent flux overestimation
* Post-fit peak consistency checks
* Robust source grouping with size limits
* Parallel processing with joblib
* Progress monitoring with heartbeat status messages
* Timing information for each processing chunk

The algorithm performs post-processing validation by comparing predicted
peak values with actual image values in each region.
"""

import sys
import time
import threading
from datetime import datetime
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import NDData, StdDevUncertainty 
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import SigmaClip
from photutils.psf import PSFPhotometry, SourceGrouper, ImagePSF
from photutils.background import LocalBackground, MMMBackground, Background2D, MedianBackground

try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False

# --- CONFIGURATION ---
SATURATION_LIMIT = 3000  #: Pixel saturation threshold for masking
READNOISE = 5.0          #: Detector read noise in electrons/ADU for uncertainty calculation

#: PSF photometry configuration parameters
PSF_CONFIG = {
    'oversampling': 1,              #: PSF oversampling factor
    'fit_shape': (31, 31),          #: Fitting box size (larger captures bright star wings)
    'fitter_maxiters': 300,         #: Maximum fitting iterations
    'aperture_radius': 4,           #: Initial flux aperture radius (smaller avoids neighbor bias)
    'pixel_scale': 0.27,            #: Pixel scale in arcsec/pixel
    'seeing_arcsec': 1.5,           #: Typical seeing FWHM in arcseconds
    'min_separation_factor': 1.5,   #: Minimum separation multiplier for grouping
    'max_group_size': 50,           #: Maximum sources per group (prevents slowdown)
}

#: Background estimation configuration
BACKGROUND_CONFIG = {
    'box_size': (50, 50),           #: Background mesh box size
    'filter_size': (3, 3),          #: Median filter size for background
    'sigma_clip_sigma': 3.0,        #: Sigma clipping threshold
    'local_inner_radius': 12,       #: Local background annulus inner radius
    'local_outer_radius': 22,       #: Local background annulus outer radius
}

# --- UTILITIES ---

class RobustSourceGrouper(SourceGrouper):
    """
    Enhanced SourceGrouper with maximum group size enforcement.
    
    Extends photutils SourceGrouper to prevent exponential slowdown
    in dense stellar clusters by splitting oversized groups.
    
    Parameters
    ----------
    min_separation : float
        Minimum separation distance for grouping sources (pixels).
    max_group_size : int, optional
        Maximum number of sources allowed in a single group (default: 50).
        Groups exceeding this size are split into smaller subgroups.
    
    Attributes
    ----------
    max_group_size : int
        Maximum allowed group size.
    
    Examples
    --------
    >>> grouper = RobustSourceGrouper(min_separation=5.0, max_group_size=30)
    >>> groups = grouper(x_positions, y_positions)
    """
    
    def __init__(self, min_separation, max_group_size=50):
        super().__init__(min_separation=min_separation)
        self.max_group_size = max_group_size

    def __call__(self, x, y):
        """
        Assign group IDs to sources, splitting oversized groups.
        
        Parameters
        ----------
        x : array_like
            X coordinates of sources.
        y : array_like
            Y coordinates of sources.
        
        Returns
        -------
        ndarray
            Group ID for each source.
        """
        groups = super().__call__(x, y)
        unique_groups, counts = np.unique(groups, return_counts=True)
        new_groups = groups.copy()
        next_id = np.max(groups) + 1
        
        for g_id, count in zip(unique_groups, counts):
            if count > self.max_group_size:
                idx = np.where(groups == g_id)[0]
                sort_idx = idx[np.argsort(x[idx])]
                for i in range(self.max_group_size, count, self.max_group_size):
                    sub_idx = sort_idx[i:i+self.max_group_size]
                    new_groups[sub_idx] = next_id
                    next_id += 1
        return new_groups


class HeartbeatPrinter(threading.Thread):
    """
    Background thread that prints periodic status messages.
    
    Provides visual confirmation that long-running processes haven't stalled
    by printing timestamped heartbeat messages at regular intervals.
    
    Parameters
    ----------
    chunk_id : int
        Identifier for the processing chunk being monitored.
    interval : int, optional
        Time between heartbeat messages in seconds (default: 900 = 15 minutes).
    
    Attributes
    ----------
    chunk_id : int
        The chunk identifier.
    interval : int
        Heartbeat interval in seconds.
    stop_event : threading.Event
        Event used to signal thread termination.
    
    Examples
    --------
    >>> heartbeat = HeartbeatPrinter(chunk_id=1, interval=600)
    >>> heartbeat.start()
    >>> # ... do work ...
    >>> heartbeat.stop()
    """
    
    def __init__(self, chunk_id, interval=900): 
        super().__init__()
        self.chunk_id = chunk_id
        self.interval = interval
        self.stop_event = threading.Event()
        self.daemon = True

    def run(self):
        """Execute the heartbeat loop until stopped."""
        while not self.stop_event.wait(self.interval):
            print(f"[{datetime.now().strftime('%H:%M:%S')}] HEARTBEAT: Chunk {self.chunk_id} is processing...")
            sys.stdout.flush()

    def stop(self):
        """Stop the heartbeat thread."""
        self.stop_event.set()


# --- CORE ROUTINES ---

def create_psf_photometry_objects(image_data, psf_data, global_bkg_median, readnoise=None):
    """
    Create and configure PSF photometry objects.
    
    Sets up the PSF model, NDData structure, background estimator, and
    source grouper needed for PSF photometry.
    
    Parameters
    ----------
    image_data : ndarray
        2D image data array.
    psf_data : ndarray
        2D PSF model array (should be normalized).
    global_bkg_median : float
        Median background level for Poisson uncertainty estimation.
    readnoise : float, optional
        Detector read noise in electrons/ADU (default: uses READNOISE config).
    
    Returns
    -------
    psf_phot : PSFPhotometry
        Configured PSF photometry object.
    nddata : NDData
        Image data wrapped with uncertainty information.
    
    Notes
    -----
    The function configures:
    
    * ImagePSF model with specified oversampling
    * Local background estimation using MMM algorithm
    * Robust source grouping with separation constraints
    * Levenberg-Marquardt least-squares fitter
    
    **Uncertainty Calculation**: Uses Poisson statistics based on the data
    values plus read noise: uncertainty = sqrt(data + readnoise²). This is
    appropriate when local background subtraction is performed, as it avoids
    double-counting background variations that would lead to unrealistically
    low reduced chi-squared values.
    """
    if readnoise is None:
        readnoise = READNOISE
        
    psf_model = ImagePSF(psf_data, oversampling=PSF_CONFIG['oversampling'])
    
    # Use Poisson + read noise for uncertainties (not background RMS)
    # This avoids double-counting background variations that are removed by local background subtraction
    error = np.sqrt(np.maximum(image_data - global_bkg_median, 0) + readnoise**2)
    nddata = NDData(data=image_data, uncertainty=StdDevUncertainty(error))
    
    bkg_estimator = LocalBackground(BACKGROUND_CONFIG['local_inner_radius'], 
                                    BACKGROUND_CONFIG['local_outer_radius'], 
                                    MMMBackground())
    
    min_sep = PSF_CONFIG['min_separation_factor'] * (PSF_CONFIG['seeing_arcsec'] / PSF_CONFIG['pixel_scale'])
    grouper = RobustSourceGrouper(min_separation=min_sep, max_group_size=PSF_CONFIG['max_group_size'])
    
    psf_phot = PSFPhotometry(
        psf_model=psf_model,
        fit_shape=PSF_CONFIG['fit_shape'],
        grouper=grouper,
        fitter=LevMarLSQFitter(),
        fitter_maxiters=PSF_CONFIG['fitter_maxiters'],
        localbkg_estimator=bkg_estimator,
        aperture_radius=PSF_CONFIG['aperture_radius']
    )
    return psf_phot, nddata


def process_single_chunk(chunk_info, image_data, psf_data, sources_x, sources_y, global_bkg_median, readnoise=None):
    """
    Process a single chunk of sources with PSF photometry.
    
    Performs PSF fitting on a subset of sources, including saturation masking
    and post-fit validation of peak consistency. Reports timing information.
    
    Parameters
    ----------
    chunk_info : tuple
        Tuple of (chunk_index, start_idx, end_idx, n_chunks) defining the chunk.
    image_data : ndarray
        2D image data array.
    psf_data : ndarray
        Normalized 2D PSF model array.
    sources_x : ndarray
        X coordinates of all sources.
    sources_y : ndarray
        Y coordinates of all sources.
    global_bkg_median : float
        Median background level for uncertainty calculation.
    readnoise : float, optional
        Detector read noise in electrons/ADU (default: uses READNOISE config).
    
    Returns
    -------
    Table or None
        Photometry results table with additional columns:
        
        * model_peak : Predicted peak pixel value from fitted flux
        * data_peak : Actual peak pixel value in image
        * peak_ratio : Ratio of model_peak to data_peak
        
        Returns None if processing fails.
    
    Notes
    -----
    **Timing Information**: The function prints:
    
    * Start time and source index range
    * End time and elapsed duration
    * Timing breakdown: initialization, fitting, and validation
    
    **Saturation Handling**: Pixels above SATURATION_LIMIT are masked
    during fitting to prevent flux overestimation from saturated cores.
    
    **Peak Validation**: Post-fit check compares the model-predicted peak
    with the actual data peak to identify problematic fits.
    """
    i, start_idx, end_idx, n_chunks = chunk_info
    chunk_id = i + 1
    
    # Start timing
    t_start = time.time()
    start_time_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    heartbeat = HeartbeatPrinter(chunk_id=chunk_id)
    heartbeat.start()
    
    print(f"\n>>> [CHUNK {chunk_id}/{n_chunks}] START: {start_idx} to {end_idx}")
    print(f"    Start time: {start_time_str}")
    
    try:
        # Phase 1: Initialization
        t_init_start = time.time()
        psf_phot, nddata = create_psf_photometry_objects(image_data, psf_data, global_bkg_median, readnoise)
        
        # MASKING: Ignore saturated pixels in the fit to improve wing-based flux estimation
        mask = image_data > SATURATION_LIMIT
        
        init_params = Table()
        init_params['x_init'] = sources_x[start_idx:end_idx]
        init_params['y_init'] = sources_y[start_idx:end_idx]
        t_init_elapsed = time.time() - t_init_start
        
        # Phase 2: PSF Fitting
        t_fit_start = time.time()
        phot_chunk = psf_phot(nddata, init_params=init_params, mask=mask)
        t_fit_elapsed = time.time() - t_fit_start

        # Phase 3: POST-FIT Peak Consistency Check
        t_validate_start = time.time()
        psf_peak_val = np.max(psf_data)
        peak_ratios = []
        data_peaks  = []
        model_peaks = []
        
        for row in phot_chunk:
            try:
                ix, iy = int(round(row['x_fit'])), int(round(row['y_fit']))
                data_peak = np.max(image_data[iy-1:iy+2, ix-1:ix+2])-row['local_bkg']
                data_peaks.append(data_peak)
                model_peak = row['flux_fit'] * psf_peak_val
                model_peaks.append(model_peak)
                peak_ratios.append(model_peak / data_peak)
            except:
                peak_ratios.append(1.0)
                data_peaks.append(-99.)
                model_peaks.append(-99.)
        
        phot_chunk['model_peak'] = model_peaks
        phot_chunk['data_peak'] = data_peaks
        phot_chunk['peak_ratio'] = peak_ratios
        t_validate_elapsed = time.time() - t_validate_start
        
        # Total timing
        t_end = time.time()
        t_total = t_end - t_start
        end_time_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        print(f"<<< [CHUNK {chunk_id}/{n_chunks}] COMPLETE")
        print(f"    End time: {end_time_str}")
        print(f"    Total elapsed: {t_total:.2f} seconds")
        print(f"    Timing breakdown: Init={t_init_elapsed:.2f}s, Fit={t_fit_elapsed:.2f}s, Validate={t_validate_elapsed:.2f}s")
        
        return phot_chunk

    except Exception as e:
        t_end = time.time()
        t_total = t_end - t_start
        end_time_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"!!! [CHUNK {chunk_id}] ERROR: {str(e)}")
        print(f"    Error occurred at: {end_time_str}")
        print(f"    Time before error: {t_total:.2f} seconds")
        return None
    finally:
        heartbeat.stop()


def do_crowded(image, psf, source_table, outroot='psf_phot', n_processes=1, readnoise=None):
    """
    Perform PSF photometry on a crowded field image.
    
    Main processing function that orchestrates parallel PSF photometry
    on large source catalogs, with chunked processing and optional
    multi-core parallelization.
    
    Parameters
    ----------
    image : str
        Path to FITS image file.
    psf : str
        Path to FITS PSF model file.
    source_table : str
        Path to FITS table containing source positions.
        Must have 'xcentroid' and 'ycentroid' columns.
    outroot : str, optional
        Output filename root (default: 'psf_phot').
        Results saved as ``{outroot}.fits``.
    n_processes : int, optional
        Number of parallel processes (default: 1).
        Set to 0 for auto-detection (uses CPU count - 1).
        Requires joblib for n_processes > 1.
    readnoise : float, optional
        Detector read noise in electrons/ADU (default: uses READNOISE config).
        Used for uncertainty calculation.
    
    Returns
    -------
    None
        Results are written to disk.
    
    Notes
    -----
    **Processing Steps**:
    
    1. Load image, PSF, and source catalog
    2. Estimate global background median (for uncertainty calculation)
    3. Split sources into chunks of 5000
    4. Process chunks in parallel (if enabled)
    5. Validate results with peak consistency checks
    6. Write combined results to FITS table
    
    **Uncertainty Calculation**: Uses Poisson statistics (sqrt(counts)) plus
    read noise, rather than background RMS. This is appropriate because 
    PSFPhotometry performs local background subtraction, which removes the
    spatial background variations that would otherwise be included in the RMS.
    Using background RMS would lead to unrealistically low reduced chi-squared
    values due to overestimated uncertainties.
    
    **Output Columns**:
    
    Standard PSFPhotometry columns plus:
    
    * model_peak : Model-predicted peak value
    * data_peak : Measured peak value
    * peak_ratio : Quality metric (model/data)
    
    **Performance**: Chunk size of 5000 balances memory usage
    and parallelization efficiency. Smaller chunks increase overhead;
    larger chunks increase memory footprint.
    
    Examples
    --------
    >>> # Single-threaded processing
    >>> do_crowded('image.fits', 'psf.fits', 'sources.fits')
    
    >>> # Parallel processing with 4 cores
    >>> do_crowded('image.fits', 'psf.fits', 'sources.fits', 
    ...            outroot='my_phot', n_processes=4)
    
    >>> # Auto-detect CPU count with custom read noise
    >>> do_crowded('image.fits', 'psf.fits', 'sources.fits', 
    ...            n_processes=0, readnoise=7.0)
    """
    image_data = fits.getdata(image)
    psf_data = fits.getdata(psf)
    psf_data /= np.sum(psf_data)
    sources = Table.read(source_table)
    
    if readnoise is None:
        readnoise = READNOISE
    
    # Global background estimation - use median for uncertainty calculation
    print("Estimating global background...")
    bkg = Background2D(image_data, BACKGROUND_CONFIG['box_size'], 
                      filter_size=BACKGROUND_CONFIG['filter_size'],
                      sigma_clip=SigmaClip(sigma=3.0), 
                      bkg_estimator=MedianBackground())
    global_bkg_median = np.median(bkg.background)
    print(f"Global background median: {global_bkg_median:.2f} ADU")
    print(f"Using read noise: {readnoise:.2f} e-/ADU")
    
    chunk_size = 5000 
    n_sources = len(sources)
    n_chunks = int(np.ceil(n_sources / chunk_size))
    chunk_infos = [(i, i*chunk_size, min((i+1)*chunk_size, n_sources), n_chunks) for i in range(n_chunks)]
    
    sources_x, sources_y = np.array(sources['xcentroid']), np.array(sources['ycentroid'])

    if n_processes == 0:
        import multiprocessing
        n_processes = max(1, multiprocessing.cpu_count() - 1)
        print(f"Auto-detected {n_processes} processes")

    print(f"\nProcessing {n_sources} sources in {n_chunks} chunks with {n_processes} process(es)")
    
    overall_start = time.time()
    
    if n_processes > 1 and HAS_JOBLIB:
        print("Using parallel processing with joblib")
        all_results = Parallel(n_jobs=n_processes)(
            delayed(process_single_chunk)(info, image_data, psf_data, sources_x, sources_y, global_bkg_median, readnoise)
            for info in chunk_infos
        )
    else:
        if n_processes > 1 and not HAS_JOBLIB:
            print("WARNING: joblib not available, falling back to serial processing")
        print("Using serial processing")
        all_results = [process_single_chunk(info, image_data, psf_data, sources_x, sources_y, global_bkg_median, readnoise) 
                      for info in chunk_infos]

    overall_elapsed = time.time() - overall_start

    valid_results = [r for r in all_results if r is not None]
    if valid_results:
        final_table = vstack(valid_results)
        out_name = f"{outroot if outroot else 'PhotResult'}.fits"
        final_table.write(out_name, overwrite=True)
        print(f"\n{'='*70}")
        print(f"DONE. Results saved to {out_name}")
        print(f"Total processing time: {overall_elapsed:.2f} seconds ({overall_elapsed/60:.2f} minutes)")
        print(f"Processed {len(final_table)} sources successfully")
        print(f"{'='*70}")


def steer(argv):
    """
    Parse command-line arguments and execute photometry.
    
    Parameters
    ----------
    argv : list
        Command-line argument list (typically sys.argv).
    
    Command-Line Options
    --------------------
    -out <root>
        Output filename root (default: 'psf_phot')
    -np <N>
        Number of parallel processes (default: 1, 0=auto)
    
    Positional Arguments
    --------------------
    image.fits
        Input image FITS file
    psf.fits
        PSF model FITS file
    stars.fits
        Source catalog FITS table
    
    Examples
    --------
    Command line usage::
    
        $ python PsfPhot.py image.fits psf.fits stars.fits
        $ python PsfPhot.py -out myresults -np 4 image.fits psf.fits stars.fits
        $ python PsfPhot.py -np 0 image.fits psf.fits stars.fits  # auto-detect CPUs
    """
    image, psf, stars, root = '', '', '', ''
    n_processes = 1
    i = 1
    while i < len(argv):
        if argv[i] == '-out': 
            i += 1
            root = argv[i]
        elif argv[i] == '-np': 
            i += 1
            n_processes = int(argv[i])
        elif image == '': 
            image = argv[i]
        elif psf == '': 
            psf = argv[i]
        elif stars == '': 
            stars = argv[i]
        i += 1
    
    if all([image, psf, stars]):
        do_crowded(image, psf, stars, outroot=root, n_processes=n_processes)
    else:
        print("ERROR: Missing required arguments")
        print("Usage: PsfPhot.py [-out root] [-np N] image.fits psf.fits stars.fits")


if __name__ == "__main__":
    if len(sys.argv) > 1: 
        steer(sys.argv)
    else: 
        print(__doc__)
        print("\nUsage: PsfPhot.py [-out root] [-np N] image.fits psf.fits stars.fits")
        print("\nOptions:")
        print("  -out <root>  Output filename root (default: 'psf_phot')")
        print("  -np <N>      Number of processes (default: 1, 0=auto-detect)")
