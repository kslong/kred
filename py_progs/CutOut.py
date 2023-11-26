#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create cutout images of objects contained in a "masterfile" 
in the images contained in a specific directory or set of directories


Command line usage (if any):

    usage: CutOut.py [-r] [-size 10] image_directory object_table

Description:  

    where 

        image_directory is the directory to be searched for fits images
        object_table is a list of objects with their positions

    and
        -h prints the documentatio and quits

        -r is an optional parameter that says not
        only to serach the directory names above, but also
        all of the subdirectories

        -size 10 is an optional parameter that gives the size
        of the cutouts in arc minutes


Primary routines:

    doit

Notes:

    The object_table should be an astropy table that contains (at least)
    the following columns:

    Source_name   - the name of the object (one word, with no spaces)
    RA 
    Dec
                                       
History:

230930 ksl Coding begun

'''




import os
from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
from glob import glob

def extract_region(source_name, ra, dec, size_arcmin, input_fits, outdir='test',default_value=0):
    # Open the input FITS file

    print('Starting %s  %f %f\n' % (source_name,ra,dec))
    hdul = fits.open(input_fits)
    
    # Extract WCS information
    wcs = WCS(hdul[0].header)
    
    # Convert RA and Dec to pixel coordinates
    coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
    x, y = wcs.all_world2pix(coords.ra.deg, coords.dec.deg, 0)
    
    # Convert size from arcminutes to pixels
    size_pixels = (size_arcmin / 60) / abs(wcs.wcs.cd[0, 0])
    
    # Define the region to extract
    xmin = int(max(0, x - size_pixels/2))
    xmax = int(min(hdul[0].data.shape[1], x + size_pixels/2))
    ymin = int(max(0, y - size_pixels/2))
    ymax = int(min(hdul[0].data.shape[0], y + size_pixels/2))

    # print('test1: ',xmin,xmax,ymin,ymax,xmax-xmin,ymax-ymin)

    if xmax-xmin>ymax-ymin:
        xmax=xmin+ymax-ymin
    elif ymax-ymin>xmax-xmin:
        ymax=ymin+xmax-xmin
    
    # Check if the position is within the image
    if not (0 <= x < hdul[0].data.shape[1] and 0 <= y < hdul[0].data.shape[0]):
        # print("RA and Dec not within the image.")
        hdul.close()
        return
    
    # Extract the region
    extracted_data = hdul[0].data[ymin:ymax, xmin:xmax]

    
    # Create a new array with the fixed size
    output_data = np.full((int(size_pixels), int(size_pixels)), default_value, dtype=extracted_data.dtype)

    # Insert the extracted data into the new array
    output_data[:extracted_data.shape[0], :extracted_data.shape[1]] = extracted_data
    
    # Update WCS information for the output image
    new_center_ra, new_center_dec = wcs.all_pix2world(xmin + size_pixels/2, ymin + size_pixels/2, 0)
    wcs_output = wcs.deepcopy()
    wcs_output.wcs.crval = [new_center_ra, new_center_dec]
    wcs_output.wcs.crpix = [size_pixels/2, size_pixels/2]  # Update reference pixel coordinates
    wcs_output.wcs.cd = wcs.wcs.cd  # Copy CD matrix for rotation
    
    # Update FITS header with the new WCS information
    header = wcs_output.to_header()
    
    # Create a new FITS file with the extracted data and updated WCS
    hdu = fits.PrimaryHDU(output_data, header=header)
    hdul_out = fits.HDUList([hdu])


    os.makedirs(outdir,exist_ok=True)

    word=input_fits.split('/')
    output_fits='%s/%s_%s' % (outdir,source_name,word[-1])

    
    print('Writing   %s ' % (output_fits))
    # Write to the output FITS file
    hdul_out.writeto(output_fits, overwrite=True)
    
    # Close both FITS files
    hdul.close()
    hdul_out.close()


def make_cutouts(object_table='lmc_snr.txt',size_arcmin=10,directory='T07',subdirectories=True,xoutdir=''):
    '''
    where 
        object table is an ascii table containing the names and positions of 
        a specific object

        size_in_arcmin  is the size of the region to be extracted

        director is the directory to search

        subdirectories indicates whether just to search the top level directory, or to
        search all subdirectories as well

    '''



    try:
        xtab=ascii.read(object_table)
    except:
        print('Could not read %s' % object_table)
        return

    if xoutdir=='':
        words=object_table.split('/')
        xoutdir=words[-1].replace('.txt','')
        xoutdir='%s_%s' % (xoutdir,'CutOuts')
    



    if subdirectories==True:
        pattern = '%s/**/*.fits*' % directory
        files=glob(pattern,recursive=True)
    else:
        pattern = '%s/*.fits*' % directory
        files=glob(pattern)

    if len(files)==0:
        print('Found no files in %s' % directory)
        return 
    else:
        print('Found %d files to search' % len(files))


    for one_file in files:
        for one_object in xtab:
            source_name=one_object['Source_name']
            ra=one_object['RA']
            dec=one_object['Dec']
            extract_region(source_name,ra,dec,size_arcmin, one_file,  outdir=xoutdir,default_value=-99)



def steer(argv):
    '''
    Just a steering routine
    '''

    recursive=False
    size=10
    xdir=''
    xobject_table=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-r':
            recursive=True
        elif argv[i]=='-size':
            i+=1
            size=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch  %s' % argv[i])
            return
        elif xdir=='':
            xdir=argv[i]
        elif xobject_table=='':
            xobject_table=argv[i]
        else:
            print(__doc__)
            print('Error: Too many arguments')
        i+=1

    if xdir=='' or xobject_table=='':
        print(__doc__)
        print('Not enough arguments')
        return


    make_cutouts(object_table=xobject_table,size_arcmin=size,directory=xdir,subdirectories=recursive)







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)




