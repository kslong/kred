#!/usr/bin/env python 
# coding: utf-8

"""BackStats - Background Statistics Calculator

Space Telescope Science Institute

Command Line Usage
------------------

::

    BackStats.py [-np 8] [-all] [-rm] [-sigma_clipped] field T01 T02 ...

**Required Arguments:**

field
    Field name (e.g., 'LMC_c42')

T01 T02 ...
    One or more tile identifiers (unless -all is specified)

**Optional Arguments:**

-all
    Process all tiles (T01 through T16)

-rm
    Remove files created by BackPrep after successful completion. Only
    occurs at the end of the program so that if something fails, no
    background files will have been deleted.

-np N
    Perform calculations in parallel using N threads (default: 1)

-sigma_clipped
    Diagnostic mode that calculates sigma-clipped versions of mean,
    median, and std. Greatly increases run time but provides robust
    statistics.

Output
------
Results written to ``Summary/{field}_{tile}_xxx.txt`` containing:

* Overlapping file pairs
* Number of overlapping pixels
* Mean, median, std of differences
* Sigma-clipped statistics (if requested)
* Skewness and kurtosis
* Mode (half-sample algorithm)
* Filter and exposure time metadata

Performance Notes
This routine can take considerable time because:

* Each tile involves many files (often 100-1000)
* Must check all pairwise combinations (N × (N-1) matches)
* For 100 files: ~10,000 comparisons
* For 1000 files: ~1,000,000 comparisons

**Optimization:**

Running with maximum available threads is strongly recommended. Use ``-np 0``
for auto-detection or explicitly set thread count based on your system.

Background Level Determination
The mode (calculated using half-sample algorithm) is used for setting
background levels between images. This is more robust than mean or median
for skewed distributions common in astronomical images.

Notes
-----

This routine can take considerable time because:

* Each tile involves many files (often 100-1000)
* Must check all pairwise combinations (N × (N-1) matches)
* For 100 files: ~10,000 comparisons
* For 1000 files: ~1,000,000 comparisons

**Optimization:**

Running with maximum available threads is strongly recommended. Use ``-np 0``
for auto-detection or explicitly set thread count based on your system.

Background Level Determination
The mode (calculated using half-sample algorithm) is used for setting
background levels between images. This is more robust than mean or median
for skewed distributions common in astronomical images.

Version History
---------------

230505 ksl

    Coding begun and routine parallelized

230617 ksl

    Adapted to new version of routines

230704 ksl

    Mode and additional statistics added. Mode now used for background levels.

Author
Space Telescope Science Institute
"""


import sys
import os
import psutil
from glob import glob
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
from astropy.stats import sigma_clipped_stats
from multiprocessing import Pool
from log import *
import shutil
from medianrange import *
import timeit
import gc
import time

#: Global error flag for tracking errors during processing
ierror = False

#: Background directory name
BACKDIR = 'DECam_BACK'


def halfsamplemode(inputData, axis=None):
    """
    Compute mode of an array using the half-sample algorithm.
    
    This is a robust estimator of the mode that works well for
    approximately unimodal distributions.

    Parameters
    ----------
    inputData : ndarray
        Input array of values.
    axis : int, optional
        Axis along which to compute the mode. If None, computes mode
        of flattened array. Default: None.

    Returns
    -------
    float
        Mode estimate.
    
    Notes
    -----
    **Algorithm:**
    
    1. Sort the data
    2. Find the narrowest half-sample (shortest interval containing
       half the data)
    3. Recursively apply to the narrowest half
    4. Continue until 1-2 values remain
    
    This algorithm is more robust than simple binning approaches and
    works well for skewed distributions common in background differences.
    
    **Computational Complexity:**
    
    O(n log n) due to initial sort, then O(log n) iterations.
    
    Examples
    --------
    >>> data = np.random.normal(5, 1, 1000)
    >>> mode = halfsamplemode(data)
    >>> print(f"Mode: {mode:.2f}")
    """
    if axis is not None:
        return np.apply_along_axis(halfsamplemode, axis, inputData)

    vals = np.sort(inputData.ravel())
    n = len(vals)
    while n > 2:
        nhalf = (n+1)//2
        isub = (vals[n-nhalf:]-vals[:nhalf]).argmin()
        vals = vals[isub:isub+nhalf]
        n = nhalf
    if n == 2:
        return 0.5*(vals[0]+vals[1])
    else:
        return vals[0]


def quartile_stats(xdata):
    """
    Compute skewness and kurtosis with outlier rejection.
    
    Uses quartile-based outlier rejection before computing higher-order
    moments, providing robust statistics for distributions with outliers.

    Parameters
    ----------
    xdata : ndarray
        Input data array.

    Returns
    -------
    skew : float
        Skewness (third standardized moment).
    kurt : float
        Kurtosis (fourth standardized moment).
    
    Notes
    -----
    **Outlier Rejection:**
    
    Uses the IQR (interquartile range) method:
    
    * Q1 = 25th percentile
    * Q3 = 75th percentile
    * IQR = Q3 - Q1
    * Outliers: values < Q1 - 1.5×IQR or > Q3 + 1.5×IQR
    
    After outlier removal, computes standard moments.
    
    **Reference:**
    
    See: https://towardsdatascience.com/skewness-and-kurtosis-with-outliers-f43167532c69
    
    **Interpretation:**
    
    * Skewness = 0: symmetric distribution
    * Skewness > 0: right-skewed (tail extends right)
    * Skewness < 0: left-skewed (tail extends left)
    * Kurtosis = 3: normal distribution
    * Kurtosis > 3: heavy tails
    * Kurtosis < 3: light tails
    
    Examples
    --------
    >>> data = np.random.normal(0, 1, 10000)
    >>> skew, kurt = quartile_stats(data)
    >>> print(f"Skewness: {skew:.3f}, Kurtosis: {kurt:.3f}")
    """
    vals = np.sort(xdata.ravel())

    # Eliminate outliers using IQR method
    n = len(vals)
    q3 = vals[int(0.75*n)]
    q1 = vals[int(0.25*n)]
    iqr = q3 - q1
    xmin = q1 - 1.5*iqr
    xmax = q3 + 1.5*iqr

    xvals = vals[vals > xmin]
    xxvals = xvals[xvals < xmax]

    std = np.std(xxvals)
    mean = np.average(xxvals)
    med = np.median(xxvals)

    skew = (xxvals - mean)**3 / std**3
    skew = np.average(skew)

    kurt = (xxvals - mean)**4 / std**4
    kurt = np.average(kurt)

    return skew, kurt


def calculate_one(file, xmatch_files, indir, calc_sigma_clipped=False, npix_min=100):
    """
    Calculate background statistics for one file vs multiple matches.
    
    Computes statistics on the difference in overlap regions between one
    base file and a list of comparison files.

    Parameters
    ----------
    file : str
        Base FITS filename.
    xmatch_files : list of str
        List of FITS filenames to compare against base file.
    indir : str
        Directory containing the FITS files.
    calc_sigma_clipped : bool, optional
        If True, calculate sigma-clipped statistics. Increases runtime
        significantly. Default: False.
    npix_min : int, optional
        Minimum number of overlapping pixels required. Pairs with fewer
        pixels are skipped. Default: 100.

    Returns
    -------
    list of list
        List of records, where each record is a list containing:
        
        * file1, file2 : str - Filenames
        * npix : int - Number of overlapping pixels
        * mean, med, std : float - Unbiased statistics
        * mean_clipped, med_clipped, std_clipped : float - Sigma-clipped stats
        * skew, kurt : float - Higher moments
        * Delta : float - Mode of difference
    
    Notes
    -----
    **Processing Steps:**
    
    1. Open base file
    2. For each comparison file:
       
       - Find overlapping pixels (both non-zero)
       - Calculate difference in overlap region
       - Compute statistics (mean, median, std)
       - Optionally compute sigma-clipped stats
       - Calculate skewness and kurtosis
       - Calculate mode using half-sample algorithm
    
    3. Close files and clean up memory
    
    **Memory Management:**
    
    Uses explicit garbage collection and file closure to handle large
    datasets without memory leaks.
    
    **Error Handling:**
    
    * Sets global ierror flag on failures
    * Logs errors for problematic file pairs
    * Skips pairs with insufficient overlap or NaN results
    
    Examples
    --------
    >>> records = calculate_one('base.fits', ['comp1.fits', 'comp2.fits'],
    ...                         'DECam_BACK/LMC_c42/T07')
    >>> print(f"Processed {len(records)} comparisons")
    """
    global ierror

    one = fits.open('%s/%s' % (indir, file))
    records = []
        
    j = 0
    while j < len(xmatch_files):
        one_record = []
        two = fits.open('%s/%s' % (indir, xmatch_files[j]))
        try:
            overlap = np.select([one[0].data*two[0].data != 0], [1], default=0)
            nonzero = np.count_nonzero(overlap)
        except:
            log_message('Error: Problem comparing %s and %s' % (file, xmatch_files[j]))
            nonzero = 0
            ierror = True

        if nonzero > npix_min:
            ok = True
            one_record = [file, xmatch_files[j], nonzero]

            # Create mask (good values = 0, masked = 1)
            xmask = overlap - 1
            one_masked = np.ma.array(one[0].data, mask=xmask)
            two_masked = np.ma.array(two[0].data, mask=xmask)
            xxdelta = two_masked - one_masked
            xdelta = xxdelta.compressed()
            
            mean = np.ma.average(xdelta)
            med = np.ma.median(xdelta)
            std = np.ma.std(xdelta)
            if np.isnan(med):
                ok = False

            if ok:
                # Add unbiased values
                one_record = one_record + [mean, med, std]
                
                # Add clipped values
                if calc_sigma_clipped: 
                    mean1, median1, std1 = sigma_clipped_stats(xxdelta, sigma_lower=2, 
                                                                sigma_upper=2, grow=3)
                    if np.isnan(median1):
                        ok = False
                else:
                    mean1 = -99.
                    median1 = -99.
                    std1 = -99.

            # Add skew and kurtosis
            if ok:
                one_record = one_record + [mean1, median1, std1]
                skew, kurt = quartile_stats(xdelta)
                if np.isnan(skew):
                    ok = False

            # Add mode
            if ok:
                one_record = one_record + [skew, kurt]
                xmode = halfsamplemode(xdelta)
                if np.isnan(xmode):
                    ok = False

            if ok:
                one_record = one_record + [xmode]
                records.append(one_record)
            else:
                print('BackStats: Failed with (NaNs for) files %s and %s' % 
                      (file, xmatch_files[j]), flush=True)

        two.close()
        del two
        gc.collect()
        j += 1
        
    one.close()
    print('BackStats: Finished %3d x-matches for %s/%s' % (len(records), indir, file), flush=True)
    return records


def do_one_tile(field='LMC_c42', tile='T07', nproc=1, calc_sigma_clipped=False):
    """
    Process all overlapping file pairs for one tile.
    
    Orchestrates background statistics calculation for all file pairs in
    a tile, with optional parallel processing.

    Parameters
    ----------
    field : str, optional
        Field name. Default: 'LMC_c42'.
    tile : str, optional
        Tile identifier. Default: 'T07'.
    nproc : int, optional
        Number of parallel processes. If 1, uses serial processing.
        Default: 1.
    calc_sigma_clipped : bool, optional
        If True, calculate sigma-clipped statistics. Default: False.

    Returns
    -------
    Table
        Astropy Table containing background statistics for all file pairs.
    
    Raises
    ------
    IOError
        If required input files or directories cannot be found.
    
    Notes
    -----
    **Input Files:**
    
    * ``Summary/{field}_{tile}_overlap.txt`` - Overlap file pairs
    * ``Summary/{field}_{tile}.txt`` - Image metadata
    * ``DECam_BACK/{field}/{tile}/`` - Background FITS files
    
    **Processing:**
    
    1. Read overlap file listing all file pairs
    2. Create work queue for all base files
    3. Process in parallel (if nproc > 1) or serially
    4. Combine results into single table
    5. Rescale by pixel area (0.2631² → 2² arcsec²)
    6. Join with metadata (filter, exposure time)
    7. Write to Summary directory
    
    **Memory Management:**
    
    Uses explicit garbage collection and careful table handling to avoid
    memory leaks in astropy join operations.
    
    **Output Scaling:**
    
    Results are scaled by (0.2631/2)² to account for pixel size difference
    between native and binned images.
    
    Examples
    --------
    >>> # Serial processing
    >>> tab = do_one_tile('LMC_c42', 'T07')
    
    >>> # Parallel processing with 8 threads
    >>> tab = do_one_tile('LMC_c42', 'T07', nproc=8)
    
    >>> # With sigma-clipped statistics
    >>> tab = do_one_tile('LMC_c42', 'T07', nproc=8, calc_sigma_clipped=True)
    """
    sumfile = 'Summary/%s_%s_overlap.txt' % (field, tile)
    try:
        x = ascii.read(sumfile)
    except:
        print('BackStats: Could not find overlap file %s' % sumfile)
        raise IOError

    # Read summary file for exposure time
    imsumfile = 'Summary/%s_%s.txt' % (field, tile)
    try:
        imsum = ascii.read(imsumfile)
        imsum = join(imsum, x['Filename', 'XEXPTIME'], join_type='left')
    except:
        print('BackStats: Could not read imsum file: %s' % imsumfile)
        raise IOError

    indir = '%s/%s/%s' % (BACKDIR, field, tile)
    if os.path.exists(indir) == False:
        print('BackStats: Could not find the directory %s with the data' % indir)
        raise IOError

    qfiles = np.unique(x['Filename'])

    records = []
    all_inputs = []
    
    i = 0
    while i < len(qfiles):
        files = x[x['Filename'] == qfiles[i]]
        all_inputs.append([qfiles[i], files['XFilename'], indir, calc_sigma_clipped])
        i += 1
    
    log_message('BackStats: %s %s has %d files to process' % (field, tile, len(all_inputs)))

    if nproc < 2:
        i = 0
        while i < len(all_inputs):
            xrecords = calculate_one(all_inputs[i][0], all_inputs[i][1], 
                                    all_inputs[i][2], all_inputs[i][3])
            records = records + xrecords
            print('Debug: calculate_one complete (%d) for %s len %d' % 
                  (i, all_inputs[i][0], len(xrecords)))
            i += 1
    else:
        with Pool(nproc) as p:
            zrecords = p.starmap(calculate_one, all_inputs)
        records = []
        print('Debug: finished starmap', flush=True)
        for one in zrecords:
            records = records + one

    log_message('Backstats: calculate_one is complete for all files in %s %s' % (field, tile))
    xrecords = np.array(records)
    
    xtab = Table(xrecords, names=['file1', 'file2', 'npix', 'mean', 'med', 'std', 
                                  'mean_clipped', 'med_clipped', 'std_clipped', 
                                  'skew', 'kurt', 'Delta'])

    os.sync()

    # Write and re-read to ensure proper formats
    print('Debug: xtab with %d lines is complete' % len(xtab), flush=True)
    xtab.write('foo.txt', format='ascii.fixed_width_two_line', overwrite=True)
    time.sleep(1.0)
    xtab = ascii.read('foo.txt')
    print('Debug: read back foo.txt with %d lines' % (len(xtab)), flush=True)

    # Rescale to account for pixel size difference
    scale_factor = (0.2631*0.2631) / (2.*2.)

    xtab['mean'] *= scale_factor
    xtab['med'] *= scale_factor
    xtab['std'] *= scale_factor
    if calc_sigma_clipped:
        xtab['mean_clipped'] *= scale_factor
        xtab['med_clipped'] *= scale_factor
        xtab['std_clipped'] *= scale_factor
    xtab['Delta'] *= scale_factor

    print('Debug: rescaled', flush=True)

    xtab.info()

    colnames = xtab.colnames.copy()

    # Convert string columns to numeric where possible
    for colname in colnames:
        col = xtab[colname]
    
        # Skip if already numeric
        if col.dtype.kind in 'biufc':
            continue
    
        print(f'Debug: Processing column {colname}...', flush=True)
    
        try:
            # Try to convert to float
            numeric_col = col.astype(float)
            xtab.remove_column(colname)
            xtab[colname] = numeric_col
            print(f'Debug: Converted {colname} to float', flush=True)
        except (ValueError, TypeError):
            print(f'Debug: Kept {colname} as string', flush=True)
            continue

    print('Debug: column conversion complete', flush=True)        
    print('Debug: removed columns', flush=True)

    ztab = xtab[colnames]
    print('Debug: created ztab', flush=True)

    # Set column formats
    ztab['mean'].format = '.3f'
    ztab['med'].format = '.3f'
    ztab['std'].format = '.3f'
    ztab['med_clipped'].format = '.3f'
    ztab['mean_clipped'].format = '.3f'
    ztab['std_clipped'].format = '.3f'
    ztab['skew'].format = '.3f'
    ztab['kurt'].format = '.3f'
    ztab['Delta'].format = '.3f'

    print('Debug: assembled ztab', flush=True)

    # Join with metadata
    imsum.rename_column('Filename', 'file1')
    
    gc.collect()
    print('Debug: analyzing tables before join', flush=True)

    # Check table sizes and memory
    print(f'Debug: ztab: {len(ztab)} rows, {len(ztab.colnames)} columns')
    print(f'Debug: imsum subset: {len(imsum)} rows, 5 columns')

    # Check for indices
    print(f'Debug: ztab indices: {len(ztab.indices) if hasattr(ztab, "indices") else 0}')
    print(f'Debug: imsum indices: {len(imsum.indices) if hasattr(imsum, "indices") else 0}')

    # Check column names and types
    print(f'Debug: ztab columns: {ztab.colnames}')
    imsum_subset = imsum['file1', 'Image', 'FILTER', 'EXPTIME', 'XEXPTIME']
    print(f'Debug: imsum subset columns: {imsum_subset.colnames}')

    # Check for common columns
    common_cols = set(ztab.colnames) & set(imsum_subset.colnames)
    print(f'Debug: Common columns for join: {common_cols}')

    if not common_cols:
        print('ERROR: No common columns found - join will fail!')
    else:
        for col in common_cols:
            print(f'Debug: {col}: ztab dtype={ztab[col].dtype}, imsum dtype={imsum_subset[col].dtype}')

    # Check memory usage
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024**2
    print(f'Debug: Current memory usage: {mem_mb:.1f} MB')

    gc.collect()

    ztab.info()
    imsum.info()
    print('Debug: before join', flush=True)

    ztab_result = join(ztab, imsum['file1', 'Image', 'FILTER', 'EXPTIME', 'XEXPTIME'], 
                      join_type='left')
    del ztab
    ztab = ztab_result
    del ztab_result
    gc.collect()

    print('Debug: joined ztab', flush=True)

    out_name = 'Summary/%s_%s_xxx.txt' % (field, tile)
    ztab.write(out_name, format='ascii.fixed_width_two_line', overwrite=True)
    print('Debug: Wrote %s with %d lines' % (out_name, len(ztab)))
    return ztab


def steer(argv):
    """
    Parse command-line arguments and execute background statistics.
    
    Main entry point for command-line execution. Handles argument parsing
    and orchestrates tile processing.

    Parameters
    ----------
    argv : list
        Command-line arguments (typically sys.argv).

    Returns
    -------
    None
        Results written to files and logged.
    
    Notes
    -----
    **Processing:**
    
    1. Parse arguments
    2. Expand -all to tiles T01-T16 if specified
    3. Open log file
    4. Process each tile with timing
    5. Optionally remove BackPrep files on success
    6. Close log
    
    **Error Handling:**
    
    If global ierror flag is set, does not remove BackPrep files even
    with -rm option.
    
    Examples
    --------
    From command line::
    
        python BackStats.py LMC_c42 T07 T08
        python BackStats.py -all -np 8 LMC_c42
        python BackStats.py -rm -sigma_clipped LMC_c42 T07
    """
    field = ''
    tiles = []
    xall = False
    nproc = 0
    xremove = False
    calc_sigma_clipped = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-all':
            xall = True
        elif argv[i] == '-rm':
            xremove = True
        elif argv[i] == '-sigma_clipped':
            calc_sigma_clipped = True
        elif argv[i] == '-np':
            i += 1
            nproc = int(argv[i])
        elif argv[i][0] == '-':
            print('BackStats: Error: Unknown switch %s' % argv[i])
            return
        elif field == '':
            field = argv[i]
        else:
            tiles.append(argv[i])
        i += 1

    if xall:
        tiles = []
        i = 1
        while i < 17:
            tiles.append('T%02d' % i)
            i += 1

    if len(tiles) == 0:
        print('BackStats: The tiles to be processed must be listed after the field, unless -all is invoked')
        
    open_log('%s.log' % field, reinitialize=False)
    
    for one in tiles:
        if calc_sigma_clipped:
            log_message('BackStats: Starting %s %s with %d processors (w/ sigma_clipped stats)' % 
                       (field, one, nproc))
        else:
            log_message('BackStats: Starting %s %s with %d processors (w/o sigma_clipped_stats)' % 
                       (field, one, nproc))

        start_time = timeit.default_timer()
        try:
            do_one_tile(field, one, nproc, calc_sigma_clipped)
        except IOError:
            log_message('BackStats: Failed on %s %s' % (field, one))
            return

        current_time = timeit.default_timer()
        log_message('BackStats: Finished %s %s in %.1f s' % (field, one, current_time-start_time))

    if ierror:
        log_message('Error: This run of BackStats.py had errors')
        log_message('Error: BackPrep files not removed')
    elif xremove:
        for one in tiles:
            xdir = '%s/%s/%s' % (BACKDIR, field, one)
            print('BackStats: There were no errors reported, removing %s' % (xdir))
            shutil.rmtree(xdir)

    close_log()
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
