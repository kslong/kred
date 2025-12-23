#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Carry out better photometry of an image
given a model psf and a group of objects
on which to carry out the photmetry


Command line usage (if any):

    Usage: PsfPhot.py -out whatever image.fits psf.fits stars.fits

Description:  

Primary routines:

    do_croweded

Notes:
                                       
History:

251217 ksl Coding begun
251223 ksl Addit to kred

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt


#!/usr/bin/env python
# coding: utf-8



import numpy as np
from astropy.io import fits
from astropy.table import Table
from photutils.psf import PSFPhotometry, SourceGrouper, ImagePSF
from photutils.background import LocalBackground, MMMBackground

from astropy.nddata import NDData, StdDevUncertainty 
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip

from astropy.table import vstack
from tqdm import tqdm  # Optional: for progress bar
import matplotlib.pyplot as plt
from photutils.psf import SourceGrouper


def do_crowded(image='artificial.fits',psf='artificial_moffat_psf.fits',source_table='artificial_all_stars.fits',outroot='psf_phot'):
    # Load your data
    image_data = fits.getdata(image)
    psf_data = fits.getdata(psf)
    print(psf_data.shape)
    # Check PSF normalization
    psf_sum = np.sum(psf_data)
    print(f"PSF sum before normalization: {psf_sum}")

    psf_data = psf_data / psf_sum

    # Create normalized PSF model
    # normalize=True ensures the PSF integrates to 1
    psf_model = ImagePSF(psf_data, oversampling=1)

    # Load your DAOStarFinder results table
    sources = Table.read(source_table)  # or however you have it stored
    # Create an NDData object for the image
    nddata = NDData(data=image_data)




    # Calculate background and RMS
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(image_data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    # Use background RMS as uncertainty
    uncertainty = bkg.background_rms

    print(f"Median background RMS: {np.median(uncertainty):.2f}")

    # Create NDData with this uncertainty
    nddata = NDData(data=image_data, uncertainty=StdDevUncertainty(uncertainty))


    sources.info()




    # Optional: Set up local background estimation
    bkg_estimator = LocalBackground(5, 10, MMMBackground())


    # Set up grouping - stars closer than this will be fit together
    # Use ~2-3× the PSF FWHM in pixels
    pixel_scale = 0.27  # arcsec/pixel (adjust to your value)
    seeing_pixels = 1.5 / pixel_scale  # Your 1.5 arcsec seeing
    min_separation = 2.5 * seeing_pixels  # ~2.5× FWHM
    min_separation = 2.0 * seeing_pixels  # ~2.5× FWHM

    grouper = SourceGrouper(min_separation=min_separation)

    # Now include the grouper in PSFPhotometry
    fitter = LevMarLSQFitter()

    psf_phot = PSFPhotometry(
        psf_model=psf_model,
        fit_shape=(25, 25),
        finder=None,
        grouper=grouper,  # ← ADD THIS
        fitter=fitter,
        fitter_maxiters=500,
        localbkg_estimator=bkg_estimator,
        aperture_radius=10
    )

    # Process in chunks of 10,000 sources
    chunk_size = 10000
    n_sources = len(sources)
    n_chunks = int(np.ceil(n_sources / chunk_size))

    all_results = []

    print(f"Processing {n_sources} sources in {n_chunks} chunks of {chunk_size}")


    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, n_sources)

        print(f"\nChunk {i+1}/{n_chunks}: sources {start_idx} to {end_idx}")

        # Create init_params for this chunk
        init_params_chunk = Table()
        init_params_chunk['x_init'] = sources['xcentroid'][start_idx:end_idx]
        init_params_chunk['y_init'] = sources['ycentroid'][start_idx:end_idx]

        # Perform photometry on this chunk
        try:
            phot_chunk = psf_phot(nddata, init_params=init_params_chunk)  # ← CHANGED HERE
            all_results.append(phot_chunk)
            print(f"  Successfully fitted {len(phot_chunk)} sources")
        except MemoryError:
            print(f"  Memory error in chunk {i+1}, trying smaller sub-chunks...")
            # If still too large, split this chunk further
            sub_chunk_size = chunk_size // 2
            for j in range(2):
                sub_start = start_idx + j * sub_chunk_size
                sub_end = min(sub_start + sub_chunk_size, end_idx)
                init_params_sub = Table()
                init_params_sub['x_init'] = sources['xcentroid'][sub_start:sub_end]
                init_params_sub['y_init'] = sources['ycentroid'][sub_start:sub_end]
                phot_sub = psf_phot(nddata, init_params=init_params_sub)  # ← AND HERE
                all_results.append(phot_sub)
        except Exception as e:
            print(f"  Error in chunk {i+1}: {e}")
            continue

    # Combine all results
    if all_results:
        phot_results = vstack(all_results)
        print(f"\n{'='*60}")
        print(f"COMPLETE: Total sources fitted: {len(phot_results)}")
        print(f"Columns: {phot_results.colnames}")

        # Quick quality check
        good_fits = np.sum(phot_results['qfit'] == 0)
        print(f"Good fits (qfit=0): {good_fits} ({100*good_fits/len(phot_results):.1f}%)")
    else:
        print("No results obtained")



    phot_results.info()

    if outroot=='':
        outroot='PhotPsf'

    phot_results.write('%s.fits' % outroot,format='fits',overwrite=True)
    print(len(sources),len(phot_results))

    # Make some plots
    plt.figure(1, (6, 6))
    plt.loglog(sources['Net'], phot_results['flux_fit'], '.', alpha=0.01)
    plt.xlim(10, 1e5)
    plt.ylim(10, 1e5)
    plt.xlabel('Aperture Net Flux')
    plt.ylabel('PSF Fit Flux')
    plt.savefig('%s_compare.png' % outroot) 
    plt.close()

    plt.figure(2, (6, 6))
    plt.semilogx(sources['Net'], phot_results['flux_fit']/sources['Net'], '.', alpha=0.01)
    plt.xlim(10, 1e5)
    plt.ylim(0, 2)
    plt.xlabel('Aperture Net Flux')
    plt.ylabel('PSF Flux / Aperture Flux')
    plt.axhline(1.0, color='r', linestyle='--', alpha=0.5)
    plt.savefig('%s_compare_ratio.png' % outroot)  
    plt.close()





def steer(argv):
    '''
    This is generally just a steering routine

    Usage: PsfPhot.py -out whatever image.fits psf.fits stars.fits
    '''

    image=''
    psf=''
    stars=''
    root=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            root=argv[i]
        elif argv[i][0]=='-':
            print('Error: Could not intepret commands: ',argv)
            return
        elif image=='':
            image=argv[i]
        elif psf=='':
            psf=argv[i]
        elif stars=='':
            stars=argv[i]
        else:
            print('Error: Improper commands: ',argv)
            return
        i+=1


    do_crowded(image=image,psf=psf,source_table=stars,outroot=root)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
