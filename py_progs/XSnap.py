#!/usr/bin/env python
# coding: utf-8

"""XSnap - Create Snapshots of Astronomical Sources

Space Telescope Science Institute

Synopsis
--------

Create PNG images (and optionally FITS cutouts) of astronomical sources from
FITS images, with optional region overlays.

Command Line Usage
------------------

::

    XSnap.py [-size arcmin] [-type suffix] [-min vmin] [-max vmax] [-o outname] image.fits master_table
    XSnap.py [-size arcmin] [-type suffix] [-min vmin] [-max vmax] snapshot_table region_table

Modes
-----

The script operates in three modes depending on the inputs:

1. **Overview mode** (FITS file + master table, no -size):
   Creates a single PNG of the full FITS image with regions from the master
   table overlaid.

   Example::

       XSnap.py -o lmc_ha_overview DECam_SWARP/LMC_c42_T01.ha.fits config/lmc_snr.txt

   Output: ``lmc_ha_overview.png``

2. **Snapshot mode** (FITS file + master table + -size):
   Creates one PNG snapshot per source in the master table, all extracted from
   the same FITS image. Also creates FITS cutouts in ``xdata/``.

   Example::

       XSnap.py -size 10 -type ha -min -1 -max 20 DECam_SWARP/LMC_c42_T01.ha.fits config/lmc_snr.txt

   Output: ``ximage/{Source_name}.ha.png`` for each source, plus FITS cutouts
   in ``xdata/``

3. **Multi-file snapshot mode** (snapshot table + region table, no FITS file):
   The snapshot table must contain a ``filename`` column specifying a different
   FITS file for each source. Creates one PNG per row using the corresponding
   FITS file.

   Example::

       XSnap.py -size 10 -type ha snapshots.txt config/lmc_snr.txt

   Where ``snapshots.txt`` contains columns: Source_name, RA, Dec, filename

Options
-------

-size arcmin    Size of snapshot cutouts in arcminutes. Required for snapshot modes.
-type suffix    Suffix appended to output filenames (e.g., "ha" -> Source.ha.png).
                Useful for distinguishing filter/image types.
-o outname      Base name for output file in overview mode (produces outname.png).
-min vmin       Minimum value for image scaling (default: 5th percentile).
-max vmax       Maximum value for image scaling (default: 95th percentile).

Input Tables
------------

Master/region tables must contain at minimum: Source_name, RA, Dec

For region overlays, tables may also include:
- RegType: "circle" or "ellipse"
- Major, Minor: region sizes in arcseconds
- Theta: position angle for ellipses

Output
------

- PNG images are written to ``ximage/`` (snapshots) or current directory (overview)
- FITS cutouts are written to ``xdata/`` (snapshot modes only)

"""


# # Create  routine to prodces a Summary Overview of SNRS in MCELS



import os
from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.wcs as wcs
from astropy import units as u

from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.pyplot as plt
import matplotlib.patches as patches





def extract_region(source_name, ra, dec, size_arcmin, input_fits, outdir='test',default_value=0,frac_off=0.1):
    '''
    Extract a region of a given size, but do not write an image if a fraction of an
    image has not data exceeds frac_off


    '''
    # Open the input FITS file

    print('Starting %s  %f %f on %s' % (source_name,ra,dec,input_fits))
    hdul = fits.open(input_fits)
    
    # Extract WCS information
    wcs = WCS(hdul[0].header)
    
    # Convert RA and Dec to pixel coordinates
    coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
    x, y = wcs.all_world2pix(coords.ra.deg, coords.dec.deg, 0)
    
    # Convert size from arcminutes to pixels
    size_pixels = np.rint((size_arcmin / 60) / abs(wcs.wcs.cd[0, 0]))
    size_pixels=np.abs(size_pixels)
    if size_pixels % 2 != 0:
        size_pixels+=1
    print('XXX ',size_pixels)
    
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
    print('XXX',ymin,ymax,xmin,xmax)
    extracted_data = hdul[0].data[ymin:ymax, xmin:xmax]

    
    try:
        # Create a new array with the fixed size
        output_data = np.full((extracted_data.shape[0], extracted_data.shape[1]), default_value, dtype=extracted_data.dtype)
        # Insert the extracted data into the new array
        output_data[:extracted_data.shape[0], :extracted_data.shape[1]] = extracted_data
    except Exception as e:
        print(f"B An error occurred on {input_fits}: {e}")
        return

    num_default_value_pixels = np.sum(output_data == default_value)
    frac_default=num_default_value_pixels/output_data.size
    if frac_default>frac_off:
        print(f"Too much of the requested image at {ra} and {dec} appears to outside the region of the fits file: {input_fits}")
        print(f"This image had {frac_default} > {frac_off} so ignoring")
        return

    # Update WCS information for the output image.  Note that there may be some
    # issues if the original file does not contain a cs matrix
    offset=-1
    new_center_ra, new_center_dec = wcs.all_pix2world(xmin + size_pixels/2+offset, ymin + size_pixels/2+offset, 0)
    wcs_output = wcs.deepcopy()
    wcs_output.wcs.crval = [new_center_ra, new_center_dec]
    wcs_output.wcs.crpix = [size_pixels/2, size_pixels/2]  # Update reference pixel coordinates
    wcs_output.wcs.cd = wcs.wcs.cd  # Copy CD matrix for rotation
    if wcs_output.wcs.cd is None:
        raise ValueError("WCS does not contain a CD matrix, which will cause problems later")
    
    # Update FITS header with the new WCS information.  relax=Ture keeps wd approach.
    header = wcs_output.to_header(relax=True)
    
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
    return output_fits




def display_fits_image(image_file, scale='linear', ymin=None, ymax=None,invert=True,masterfile='',outfile=''):
    # Open the FITS file
    try:
        hdul = fits.open(image_file)
    except:
        print('Could not open: %s' % image_file)
        return None


    # Access the image data
    data = hdul[0].data

    # Access the WCS information
    wcs_info = wcs.WCS(hdul[0].header)

    # Close the FITS file
    hdul.close()

    flattened_data = data.flatten()
    bad_data_threshold=-5000
    flattened_data[flattened_data < bad_data_threshold] = np.nan
    good_data = flattened_data[~np.isnan(flattened_data)]
    
    xmed=np.median(good_data)
    xstd=np.std(good_data)
    
    vmax=xmed+xstd
    vmin=xmed-xstd


    # Calculate the 5th and 95th percentiles
    lower_percentile = np.percentile(good_data, 5)
    upper_percentile = np.percentile(good_data, 95)

    vmin=lower_percentile
    vmax=upper_percentile

    if ymin!=None:
        vmin=ymin
    if ymax!=None:
        vmax=ymax

    print('stats (med,std)      %8.2e  %8.2e' % (xmed,xstd))
    print('stats (5 per cent)   %8.2e  %8.2e' % (lower_percentile,upper_percentile))
    print('limits (min,max)     %8.2e  %8.2e' % (vmin,vmax))
    


    # Apply scaling to the image data
    if scale == 'linear':
        scaled_data = data
    elif scale == 'log':
        # Adjust vmin and vmax for logarithmic scaling
        if vmin is not None:
            vmin = np.max([vmin, np.min(data[data > 0])])
        if vmax is not None:
            vmax = np.min([vmax, np.max(data)])
        scaled_data = np.log10(data)
    elif scale == 'sqrt':
        scaled_data = np.sqrt(data)
    else:
        raise ValueError("Invalid scale. Available options are 'linear', 'log', and 'sqrt'.")

    # Invert the colors if invert is True

    # Create a figure and axes using wcsaxes
    plt.close(1)  # Close existing figure 1 if it exists
    fig = plt.figure(1, figsize=(10, 10))  # Adjust the figure size as needed
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs_info, aspect='equal')  # Set the aspect ratio to 'equal'
    fig.add_axes(ax)

    # Display the image with specified vmin and vmax, and origin at lower left
    if invert == True:
        im = ax.imshow(scaled_data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
    else:
        im = ax.imshow(scaled_data, cmap='gray', vmin=vmin, vmax=vmax, origin='lower')


    # Add RA and Dec axis labels
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    words=image_file.split('/')
    root=words[-1].replace('.fits','')
    root=root.replace('.gz','')
    ax.set_title('%s' % (root))

    # Create a separate axis for the colorbar
    cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Adjust the position and size of the colorbar

    # Add colorbar
    fig.colorbar(im, cax=cax)

    reg_info=False
    if masterfile!='':
        xmaster=ascii.read(masterfile)
        # Check that master file has region info:
        reg_info=False
        colnames=xmaster.colnames
        for one_col in colnames:
            if one_col=='RegType':
                reg_info=True

    ## Now add regions

    if reg_info==True:
        pixel_scale=np.abs(wcs_info.pixel_scale_matrix[0,0])*3600.

        for one in xmaster:
            # print('starting\n ',one)
            ra=float(one['RA'])
            dec=float(one['Dec'])
            # print(ra,dec)

            sky_coord=SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
            pix_coord=sky_coord.to_pixel(wcs_info)
            if one['RegType']=='circle':
                radius_pixels=one['Major']/pixel_scale
                # print('hello :',pix_coord[0],pix_coord[1],radius_pixels)
                circle=patches.Circle(pix_coord,radius_pixels,edgecolor='red',facecolor='none',linewidth=1)
                ax.add_patch(circle)

                # Add a label just above the circle
                label='%s' % one['Source_name']
                ax.text(pix_coord[0], pix_coord[1] + radius_pixels + 5, label, color='red', fontsize=12, ha='center')

            elif one['RegType']=='ellipse':
                rmajor=2*one['Major']/pixel_scale
                rminor=2*one['Minor']/pixel_scale
                pa=one['Theta']+90 # This seems required to get the correct orienation
                ellipse=patches.Ellipse(pix_coord,width=rminor,height=rmajor,angle=pa,edgecolor='red',facecolor='none',linewidth=2)
                ax.add_patch(ellipse)
                label='%s' % one['Source_name']
                ax.text(pix_coord[0], pix_coord[1] + 0.5* rmajor + 5, label, color='red', fontsize=12, ha='center')
            else:
                print('Unkown region type: ', one['RegType'])
    else:
        print('No region info to add')



    # plt.show()
    if outfile=='':
        plt.savefig('foo.png')
    else:
        try:
            plt.savefig(outfile)
        except:
            print('Could not save to %s' % outfile)
    return ax






def make_one_image(filename,master,ymin,ymax,outroot=''):
    '''
    Create a plot of an image and overlay the retions
    from a master file on it.
    '''
    print('Creating an image of one file with regions')
    xm=ascii.read(master)
    words=master.split('/')
    if outroot=='':
        outfile_name=words[-1].replace('.txt','')
        outfile_name='%s.png' % outfile_name
    else:
        outfile_name='%s.png' % outroot
    display_fits_image(image_file=filename, scale='linear', ymin=ymin,ymax=ymax,invert=True,masterfile=master,outfile=outfile_name)

    return

def make_many_images(filename,master,xtype,size,ymin,ymax,frac_off=0.1):
    '''
    Create cut-outs of an image, one for each source in a masterfile
    and overlay the regions from the master file on each sanpshot.
    '''

    xm=ascii.read(master)
    print(xm)
    if os.path.isdir('ximage')==False:
        os.mkdir('ximage')
    
    for one in xm:
        # print('starting\n',one)
        xsource_name='%s' % one['Source_name']
        xra=float(one['RA'])
        xdec=float(one['Dec'])
        if xtype==None:
            outfile_name='ximage/%s.png' % xsource_name
        else:
            outfile_name='ximage/%s.%s.png' % (xsource_name,xtype)

        stamp=extract_region(xsource_name, xra, xdec, size_arcmin=size, input_fits=filename, outdir='xdata',default_value=0,frac_off=frac_off)
        display_fits_image(image_file=stamp, scale='linear', ymin=ymin,ymax=ymax,invert=True,masterfile=master,outfile=outfile_name)
    
    
    print('Creating an image for each regions')
    return



def make_many_images2(master,reg,xtype,size,ymin,ymax):
    '''
    Create cut-outs of an image, one for each line in a mastefile  
    and overlay the regions from the master file on each sanpshot.
    '''

    print('Doing make_many_images2: %s %s' % (master,reg))

    xm=ascii.read(master)
    print(xm)
    if os.path.isdir('ximage')==False:
        os.mkdir('ximage')
    
    for one in xm:
        # print('starting\n',one)
        xsource_name='%s' % one['Source_name']
        xra=float(one['RA'])
        xdec=float(one['Dec'])
        filename=one['filename']
        if xtype==None:
            outfile_name='ximage/%s.png' % xsource_name
        else:
            outfile_name='ximage/%s.%s.png' % (xsource_name,xtype)

        stamp=extract_region(xsource_name, xra, xdec, size_arcmin=size, input_fits=filename, outdir='xdata',default_value=0,frac_off=0.01)
        print('XXX ',stamp)
        display_fits_image(image_file=stamp, scale='linear', ymin=ymin,ymax=ymax,invert=True,masterfile=reg,outfile=outfile_name)
    
    
    print('Creating an image for each regions')
    return


def steer(argv):
    '''

    XSnap.py [-size 10] [-type ha] [-min -1] [-max 20] -out ha [filename or table of snaps]   master_table_of_regions

    without - size we use the full image, and just display eveything
    with a size we make images of each of the source in the master tale


    '''
    filename=''
    master=''
    size=-1
    outroot=''
    ymin=None
    ymax=None
    xtype=None
    reg=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-size':
            i+=1
            size=int(argv[i])
        elif argv[i]=='-min':
            i+=1
            ymin=float(argv[i])
        elif argv[i]=='-max':
            i+=1
            ymax=float(argv[i])
        elif argv[i]=='-type':
            i+=1
            xtype=argv[i]
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        elif argv[i][0]=='-':
            print('Error: Cannot parse command line :',argv)
        elif filename=='' and argv[i].count('fits')>0:
            filename=argv[i]
        elif master=='':
            master=argv[i]
        elif reg=='':
            reg=argv[i]
        else:
            print('Error: Too many arguments :', argv)
        i+=1

    
    try:
        xtab=ascii.read(master)
        xmaster=True
        xfilenames=False
        for one_col in xtab.colnames:
            if one_col.count('filename'):
                xfilenames=True
    except:
        print('Error: Cannot find master file: %s' % master)
        return

    if os.path.isfile(filename):
        fits_exists=True
    else:
        fits_exists=False

    if xfilenames==True and fits_exists==False:
        # This is the case where we want to read multiple fits file from the master file
        make_many_images2(master,reg,xtype,size,ymin,ymax)
        return
    if fits_exists==True and size>0:
        # This is the case where we have a single fits file, but a master file with regions indecated.
        make_many_images(filename,master,xtype,size,ymin,ymax)
        return
    if fits_exists:
        make_one_image(filename,master,ymin,ymax,outroot)
        return

    # If we have reached this point something is return
    print('Error: The masterfile did not cantaing the fits files and none was provided, so exiting')




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
