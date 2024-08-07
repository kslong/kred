#!/usr/bin/env python
# coding: utf-8
'''

Create plots of the images that are in the DECam_Swarp  after running Swarp.py, CleanStars.py 
etc 


Usage:   SwarpEval.py [-all] xdir field [tiles]


where -all will cause swarp to be run on all 16 tiles.  With these inputs, the routine will
use the files ending in _sw.tab to set up run files

and
    xdir is the top level directory containing subdirectories with fields and tiles

and 
    filed is a field

If one wants to run only 1 or a few tiles then the command will be something like

SwarpEval.py DECam_SWARP2  LMC_c42  T01 T03 


The evaluation imags will be stored in a directory eval, located in this case in DECam_SWARP2/LMC_c42


Note: 231009 - This routine and the command line interface to it functions 
differently from most other routines, in that one specifies a top level
directory.   
'''




from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.wcs as wcs
import astropy.units as u
from astropy.visualization.wcsaxes import WCSAxes
from glob import glob
from astropy.table import Table
import os
import numpy as np


exec_dir=os.getcwd()



def display_fits_image(image_file, scale='linear', invert=False, vmin=None, vmax=None,outfile=''):
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
    if invert:
        scaled_data = -scaled_data

    # Create a figure and axes using wcsaxes
    fig = plt.figure(1,figsize=(10, 10))  # Adjust the figure size as needed
    fig.clf()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs_info, aspect='equal')  # Set the aspect ratio to 'equal'
    fig.add_axes(ax)

    # Display the image with specified vmin and vmax, and origin at lower left
    im = ax.imshow(scaled_data, cmap='gray', vmin=vmin, vmax=vmax, origin='lower')

    # Add RA and Dec axis labels
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    if image_file.count('_b'): 
        xword = 'Match'
    else:
        xword= 'Orig'
    words=image_file.split('/')
    root=words[-1].replace('.fits','')
    ax.set_title('%s - Bkg = %s' % (root,xword))

    # Create a separate axis for the colorbar
    cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Adjust the position and size of the colorbar

    # Add colorbar
    fig.colorbar(im, cax=cax)


    # plt.show()
    if outfile=='':
        plt.savefig('foo.png')
    else:
        try:
            plt.savefig(outfile)
        except:
            print('Could not save to %s' % outfile)

# Example usage

def get_images(xdir='DECam_SWARP2',field='LMC_c01',tile='T01'):
    
    swarpdir='%s/%s/%s/' % (xdir,field,tile)
    xfiles=glob('%s/*.fits' % swarpdir)
    files=[]
    for one in xfiles:
        if one.count('weight')==0:
            files.append(one)
                     
    return(files)
    
    




def get_stats(xfiles):
    '''
    Characterize the stats of the images, assuing
    that 0 represents a pixel should be masked.
    '''
    root=[]
    med=[]
    mean=[]
    std=[]
    for one in xfiles:
        words=one.split('/')
        root.append(words[-1].replace('.fits',''))
        try:
            f=fits.open(one)
        except:
            print('Error: get stats could not open: ',one)
            continue
        zmask=np.select([f[0].data!=0],[0],1)
        
        z=np.ma.masked_array(f[0].data,mask=zmask)
        zz=z.compressed()

        xmed=np.median(zz)
        xmean=np.average(zz)
        xstd=np.std(zz)
        med.append(xmed)
        mean.append(xmean)
        std.append(xstd)
    
    xtab=Table([xfiles,root,med,mean,std],names=['Filename','Root','Median','Average','STD'])
    return xtab
        
        

    


def make_plots(xdirectory='DECam_SWARP2',field='LMC_c01',tile='T01'):

    
    print('Beginning %s %s' % (field,tile))
    files=get_images(xdirectory,field,tile)
    if len(files)==0:
        print('No images were found for field %s tile %s in %s' % (field,tile,xdirectory))
        return

    xtab=get_stats(files)
    xdir='Figs/%s/%s' % (xdirectory,field)
    
    if os.path.isdir(xdir)==False:
        os.makedirs(xdir)
    for one in xtab:
        xmin=one['Median']-0.1*one['STD']
        xmax=one['Median']+0.1*one['STD']
        xout='%s/%s_%s.png' % (xdir,one['Root'],tile)
        display_fits_image(one['Filename'], scale='linear', invert=False, vmin=xmin, vmax=xmax,outfile=xout)
    print('Finished %s %s' % (field,tile))
        

def steer(argv):
    '''
    This is just a steering routine for running swarp on one or more
    tiles from the command line

    '''

    xdir=''
    field=''
    tiles=[]
    xall=False
    redo=True
    bsub=False
    star_sub=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s' % argv[i])
            return
        elif xdir=='':
            xdir=argv[i]
        elif field=='':
            field=argv[i]
        else:
            tiles.append(argv[i])
        i+=1

    if xall:
        tiles=[]
        i=1
        while i<17:
            tiles.append('T%02d' % i)
            i+=1


    if os.path.isdir('%s/%s' % (xdir,field))==False:
        print(__doc__)
        print('Check inputs: %s/%s does not appear to be a directory' % (xdir,field))
        return
    if len(tiles)==0:
        print(__doc__)
        print('Check inputs: No tiles requested',argv)
        return




    for one in tiles:
        make_plots(xdir, field,one)


    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)



