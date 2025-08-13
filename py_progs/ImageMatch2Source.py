#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Find which processed images contain sources in a 
masterfile compatable table 


Command line usage (if any):

    usage: ImageMactch2Source.py [-h] [-out whatever] [-sep 33.]  [-n_closest 3]  file_tab source_tab immage_type

    where:
        -h prints this documentation and returns
        -out whatever defines the name of the output file; if this is not given 
            the name is based on the source_tab and image_tyype
        -sep 33 sets the maximum separation between a source and image (center) from a default value of 33 arcmin to
            whatever
        -n_closest cause the routine to return a match for a soruce to up to in this case 3 files
        file_tab - An astropy table produced by ImageSum.py that indicates what files are avalilabel to match
        source_tab - A table containing columns Source_name, RA, Dec (at least).  The program identifies files that
            match these postions
        image_type - The type of image that one wants to match to such as s2_sub_r

Description:  

    The purpose of this routine is to find out what fits files in a kred-like directory structure match multiple 
    sources.  The routien produces an astropy table that contains the Source_name, the matching file, the separation of
    each source from the center of the file.  Normally, one just searches for the best match, but one can also find
    the n files that are the closest match.

    One first runs ImageSum.py to produce a file that identifies the files of various types that are in the
    directory structure, and then runs this routine to produce the match file.

    The table that is returned contains the  Source_name, RA and Dec (of the source), the filename, the separation
    from the source to the center of the image, and the rank of the match (1 is closest)

    At that point one could use the output to crate snap shots that are centered on each soruce if one wished.

Primary routines:

    doit

Notes:
                                       
History:

250813 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
from astropy.table import Table, vstack, join
from astropy.coordinates import SkyCoord
from astropy import units as u

def find_closest(source_table, file_table, max_sep_arcmin, n):
    """
    Find the closest n files to each source within a maximum separation.
    
    Parameters:
    -----------
    source_table : astropy.table.Table
        Table containing sources with columns 'Source_name', 'RA', 'Dec'
    file_table : astropy.table.Table
        Table containing files with columns 'filename', 'RA', 'Dec'
    max_sep_arcmin : float
        Maximum separation allowed in arc_min
    n : int
        Maximum number of files to return for each source
    
    Returns:
    --------
    astropy.table.Table
        Table with columns: Source_name, filename, separation_deg, rank
    """

    max_sep=max_sep_arcmin/60.
    
    # Create SkyCoord objects for both tables
    source_coords = SkyCoord(
        ra=source_table['RA'] * u.degree, 
        dec=source_table['Dec'] * u.degree
    )
    
    file_coords = SkyCoord(
        ra=file_table['RA'] * u.degree, 
        dec=file_table['Dec'] * u.degree
    )
    
    # Lists to store results
    source_names = []
    filenames = []
    separations = []
    ranks = []
    
    # For each source, find the closest files
    for i, source_coord in enumerate(source_coords):
        source_name = source_table['Source_name'][i]
        
        # Calculate separations to all files
        separations_to_files = source_coord.separation(file_coords)
        
        # Find files within max_sep
        within_sep = separations_to_files <= max_sep * u.degree
        
        if np.any(within_sep):
            # Get indices of files within separation, sorted by distance
            valid_indices = np.where(within_sep)[0]
            valid_separations = separations_to_files[within_sep]
            
            # Sort by separation (closest first)
            sort_indices = np.argsort(valid_separations)
            sorted_file_indices = valid_indices[sort_indices]
            sorted_separations = valid_separations[sort_indices]
            
            # Take up to n closest files
            n_matches = min(n, len(sorted_file_indices))
            
            for j in range(n_matches):
                file_idx = sorted_file_indices[j]
                sep_deg = sorted_separations[j].degree
                
                source_names.append(source_name)
                filenames.append(file_table['filename'][file_idx])
                separations.append(sep_deg)
                ranks.append(j + 1)  # rank starts from 1
    
    # Create the result table
    result_table = Table()
    result_table['Source_name'] = source_names
    result_table['filename'] = filenames
    result_table['separation_arcmin'] = np.array(separations)*60.
    result_table['rank'] = ranks
    result_table['separation_arcmin'].format='.3f'
    
    return result_table


def doit(source_tab='smc_snr_cotton24.txt',file_tab='Image_Sum_DECam_SUB2.txt',image_type='s2_sub_r',nmatch=1,max_sep=20.):
    
    try:
        ftab=ascii.read(file_tab)
    except:
        print('Error: Could not read ',file_tab)
        return []
    fftab=ftab[ftab['Image_type']==image_type]
    if len(fftab)==0:
        print('Error: No images of type %s in %s' % (image_type,file_tab))
        xtypes=np.unique(ftab['Image_type'])
        print('The image types are:')
        for one in xtypes:
            print(one)
        return []
    try:
        stab=ascii.read(source_tab)  
    except:
        print('Error: Could not read ',source_tab)
        return []
    
    ztab=find_closest(source_table=stab, file_table=fftab, max_sep_arcmin=max_sep, n=nmatch)
    if len(ztab)==0:
        print('Error: Found no sources withing separation of %f arcmin' % max_sep)
        return []

    ztab=join(ztab,stab['Source_name','RA','Dec'],join_type='left')

    ztab=ztab['Source_name','RA','Dec','filename','separation_arcmin','rank']

    ztab['RA'].format='.5f'
    ztab['Dec'].format='.5f'

    found=np.unique(ztab['Source_name'])
    print('Of %d sources in %s, matched %d in images' % (len(stab),source_tab, len(found)))
    return ztab



def steer(argv):
    '''
    ImageMactch2Source.py [-h] [-out whatever] [-sep 33.]  [-n_closest 3]  file_tab source_tab immage_type
    '''
    file_tab=''
    source_tab=''
    image_type=''
    sep=33.
    n_closest=1
    outname=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            outname=argv[i]
        elif argv[i][:4]=='-sep':
            i+=1
            sep=eval(argv[i])
        elif argv[i][:2]=='-n':
            i+=1
            n_closest=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: unknwon of %s in comand line' % argv[i],argv)
            return
        elif file_tab=='':
            file_tab=argv[i]
        elif source_tab=='':
            source_tab=argv[i]
        elif image_type=='':
            image_type=argv[i]
        else:
            print('Error: Could not parse command line: ',argv)
            return

        i+=1

    if image_type == '':
        print('Error: Not enough argmuments in command line: ', argv)
        return

    print('Searching %s for images of type %s for matches to sources in %s' % (file_tab,image_type,source_tab))

    ztab=doit(source_tab,file_tab,image_type,nmatch=n_closest,max_sep=sep)

    if len(ztab)==0:
        print('No matches found')
        return

    if outname=='':
        outname='XX_%s.%s' % (image_type,source_tab)

    print('Writing results to %s' % (outname))

    ztab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    return









# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)

