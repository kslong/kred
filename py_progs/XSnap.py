#!/usr/bin/env python
# coding: utf-8
'''
Create one or more images and fits files of
snapshots of one or moe sets of souces in 
a masterfile

Usaage:
    XSnap.py [-size 10] [-type ha] [-min -1] [-max 20] image.fits  mastertable

    image.fits   a fitsfile, with the data in the PRIMARY header
    mastertable  a masterfile with postiions and sizes of objects in 
        a standard format

    The routine has two basic modes, if -size is not provided, an plot of
    the fist file will be made, and the regions will be ovelaid on the plot

    If -size is provided, cutouts of the image will be made, with a size 
    in arcmin given by the number that follows the size command

    The  fits images associated with each cutout will be placed in the ximage 
    directory

    The plots will be in the ximage directory

    -type ha  is just used to help name the plots

    In the absence of -min or -max, the images are autoscaled, if -xmin or -xmax
    are provided then one or the other of these values will replace
    what the autoscled values would have been 




'''

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





def extract_region(source_name, ra, dec, size_arcmin, input_fits, outdir='test',default_value=0,frac_off=0.01):
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
    if size_pixels % 2 != 0:
        size_pixels+=1
    
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
        print(f"This image had {frac_default} > {frac_off} so ignoring")
        return

    # Update WCS information for the output image
    offset=-1
    new_center_ra, new_center_dec = wcs.all_pix2world(xmin + size_pixels/2+offset, ymin + size_pixels/2+offset, 0)
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
    return output_fits




def display_fits_image(image_file, scale='linear', ymin=None, ymax=None,invert=True,masterfile='',outfile=''):
    # Open the FITS file
    try:
        hdul = fits.open(image_file)
    except:
        print('Could not open: ' % image_file)


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
    fig = plt.figure(1,figsize=(10, 10))  # Adjust the figure size as needed
    fig.clf()
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


    ## Now add regions

    if masterfile!='':
        pixel_scale=np.abs(wcs_info.pixel_scale_matrix[0,0])*3600.
        xmaster=ascii.read(masterfile)
        for one in xmaster:
            print('starting\n ',one)
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
                pa=-one['Theta']
                ellipse=patches.Ellipse(pix_coord,rmajor,rminor,angle=pa,edgecolor='red',facecolor='none',linewidth=1)
                ax.add_patch(ellipse)
                label='%s' % one['Source_name']
                ax.text(pix_coord[0], pix_coord[1] + 0.5* rmajor + 5, label, color='red', fontsize=12, ha='center')
            else:
                print('Unkown region type: ', one['RegType'])



    # plt.show()
    if outfile=='':
        plt.savefig('foo.png')
    else:
        try:
            plt.savefig(outfile)
        except:
            print('Could not save to %s' % outfile)
    return ax






def make_one_image(filename,master,ymin,ymax):
    '''
    Create a plot of an image and overlay the retions
    from a master file on it.
    '''
    print('Creating an image of one file with regions')
    xm=ascii.read(master)
    words=master.split('/')
    outfile_name=words[-1].replace('.txt','')
    outfile_name='%s.png' % outfile_name
    display_fits_image(image_file=filename, scale='linear', ymin=ymin,ymax=ymax,invert=True,masterfile=master,outfile=outfile_name)

    return

def make_many_images(filename,xtype,master,size,ymin,ymax):
    '''
    Create cut-outs of an image, one for each source in a masterfile
    and overlay the regions from the master file on each sanpshot.
    '''

    xm=ascii.read(master)
    print(xm)
    if os.path.isdir('ximage')==False:
        os.mkdir('ximage')
    
    for one in xm:
        print('starting\n',one)
        xsource_name='%s' % one['Source_name']
        xra=float(one['RA'])
        xdec=float(one['Dec'])
        if xtype==None:
            outfile_name='ximage/%s.png' % xsource_name
        else:
            outfile_name='ximage/%s.%s.png' % (xsource_name,xtype)

        stamp=extract_region(xsource_name, xra, xdec, size_arcmin=size, input_fits=filename, outdir='xdata',default_value=0,frac_off=0.01)
        display_fits_image(image_file=stamp, scale='linear', ymin=ymin,ymax=ymax,invert=True,masterfile=master,outfile=outfile_name)
    
    
    print('Creating an image for each regions')
    return


def steer(argv):
    '''

    XSnap.py [-size 10] [-type ha] [-min -1] [-max 20] filename  mastertable

    without - size we use the full image, and just display eveything
    with a size we make images of each of the source in the master tale


    '''
    filename=''
    master=''
    size=-1
    ymin=None
    ymax=None
    xtype=None

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
        elif argv[i][0]=='-':
            print('Error: Cannot parse command line :',argv)
        elif filename=='':
            filename=argv[i]
        elif master=='':
            master=argv[i]
        else:
            print('Error: Too many arguments :', argv)
        i+=1

    if os.path.isfile(filename)==False:
        print('Error: Cannot find image file: %s' % filename)
        return
    
    if os.path.isfile(master)==False:
        print('Error: Cannot find master file: %s' % master)
        return

    if size<0:
        make_one_image(filename,master,ymin,ymax)
    else:
        make_many_images(filename,xtype,master,size,ymin,ymax)
    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
