#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Setup the directories and identify the files that are required to carry
out the initial processing steps to run PrepFiles on a set of data


Command line usage (if any):

    usage: SetupTile.py [-h] [-all] [-xtab myfile.txt] field tile1 tile2 ...

    where:
        field is that name of a field, e.g LMC_c42
        tile1, tiles are the names of tiles as specified in the MC_tiles.txt file 
    and the optional prameters 
        -h prints this help 
        -all  causes the program to create talbles containing the relevant
            infomation for all of the tiles in a field
        -xtab myfile.txt causes the program to use a different tile defination file

Description:  

    This program reads a tile defination file with the folowing format and 
    to search for CCD images that should prepped to create fits files 
    to construct images in tiles of a field.  It writes the results to
    a table in the Summary directory

    The program searches for the tile definition file in the local direcotry
    and in the $KRED/config directory.

    The program works by reading the det.tab files created by MefSum.py



Primary routines:

    doit

Notes:
                                       
History:

230513 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
import os
from astropy.table import vstack, join
from glob import glob

from astropy.io import ascii
from log import *

CWD=os.getcwd()
DATADIR='%s/DECam_CCD/' % (CWD)
PREPDIR='%s/DECam_PREP/' % (CWD)


def find_files_to_prep(field='LMC_c45',tile='T07',ra_center=81.1,dec_center=-66.17,size_deg=.1):
    
    tab_file='Summary/%s_det.tab' % (field)
    if os.path.isfile(tab_file):
        det=ascii.read(tab_file)
    else:
        print('Cannot find %s' % tab_file)
        return

    mef_file='Summary/%s_mef.tab' % (field)


    if os.path.isfile(mef_file):
        mef=ascii.read(mef_file)
    else:
        print('Cannot find %s' % mef_file)
        return


    x=join(det,mef,join_type='left')

    

    print('Looking for CCD images with RA and DEC of %.5f %.5f and size of %.2f deg' % (ra_center,dec_center,size_deg))


    dec_min=dec_center-0.5*size_deg
    dec_max=dec_center+0.5*size_deg
    delta_ra=size_deg/np.cos(dec_center/57.29578)
    ra_min=ra_center-0.5*delta_ra
    ra_max=ra_center+0.5*delta_ra
    
    
    x['Dec_max'] = np.max([x['COR1DEC1'],x['COR2DEC1'],x['COR3DEC1'],x['COR4DEC1']],axis=0)
    x['Dec_min'] = np.min([x['COR1DEC1'],x['COR2DEC1'],x['COR3DEC1'],x['COR4DEC1']],axis=0)
    
    x=x[x['Dec_min']<dec_max]

    if len(x)==0:
            print('Failed to find any images with Dec less than %.5f ' % dec_max)
            return []

    x=x[x['Dec_max']>dec_min]
    if len(x)==0:
            print('Failed to find any images with Dec greater than %.5f ' % dec_min)
            return []


    
    x['RA_max'] = np.max([x['COR1RA1'],x['COR2RA1'],x['COR3RA1'],x['COR4RA1']],axis=0)
    x['RA_min'] = np.min([x['COR1RA1'],x['COR2RA1'],x['COR3RA1'],x['COR4RA1']],axis=0)
    
    x=x[x['RA_min']<ra_max]
    if len(x)==0:
            print('Failed to find any images with RA  less than %.5f ' % ra_max)
            return []


    # print(len(x))
    x=x[x['RA_max']>ra_min]
    if len(x)==0:
            print('Failed to find any images with RA  greater than %.5f ' % ra_min)
            return []


    
    x['Delta_Dec']=x['CENDEC1']-dec_center
    x['Delta_RA']=(x['CENRA1']-ra_center)/np.cos(dec_center/57.29578)
    
    x['Delta_Dec'].format='.2f'
    x['Delta_RA'].format='.2f'

    print('Found %d CCD extensions that satisfy the input criteria ' % len(x))

    xfiles=[]
    for one in x:
        one_name='%s_%s.fits' % (one['Root'],one['EXTNAME'])
        xfiles.append(one_name)

    x['Filename']=xfiles
    x['Field']=field
        
    
    return x['Field','Filename','Root','EXTNAME','CENRA1','CENDEC1','COR1RA1','COR1DEC1','COR2RA1','COR2DEC1','COR3RA1','COR3DEC1','COR4RA1','COR4DEC1','Delta_Dec','Delta_RA','FILTER','EXPTIME']



def setup_one_tile(field='LMC_c42',tile='T07',ra=81.108313,dec=-66.177280,size_deg=0.67):
    

    
        
    x=find_files_to_prep(field,tile,ra,dec,size_deg)

    if len(x)>0:
        outfile='Summary/%s_%s.txt' % (field,tile)
        x.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('Failed to set up %s %s' % (field,tile))

    data_dir='%s/%s/data' % (DATADIR,field)

    if os.path.isdir(data_dir) == False:
        print('Error: %s does not appear to exist' % data_dir)
        print('Run PrepMef on this field first')
        raise IOError

    tile_dir='%s/%s/%s/' % (PREPDIR,field,tile)

    if os.path.isdir(tile_dir)==False:
        print('Creating the Prep directory as %s' % tile_dir)
        os.makedirs(tile_dir)
    else:
        print('The Prep directory %s already exists' % tile_dir )
        xfiles=glob('%s/*fits*' % tile_dir)
        if len(xfiles):
            print('Removing existing fits files or links and reinitializing')
            for one in xfiles:
                os.remove(one)

    nerrors=0
    for one in x:
        one_file=one['Filename']
        data_file='%s/%s' % (data_dir,one_file)
        tile_file='%s/%s' % (tile_dir,one_file)
        if os.path.isfile(data_file)==False:
            print('Failed to find %s in %s ' % (one_file,data_dir))
            nerrors+=1
        else:
            os.symlink(data_file,tile_file)
    if nerrors:
        print('Failed to find %d of %d files in %s' % (nerrors,len(x),data_dir))
        raise IOError
    else:
        print('Successfully linked files from %s to %s' % (data_dir,tile_dir))
            





        
    return

def setup_tiles(xtab):
    '''
    Run setup_one_tile multiple times, with the possibility of parallel
    processing
    '''

    for one in xtab:
        try:
            setup_one_tile(one['Field'],one['Tile'],one['RA'],one['Dec'],one['Size'])
        except IOError:
            print('Error: Something is wrong, and must be sorted before continuing')
            print('Try rerunning MefPrep.py -finish, and then repeat this step')
            return

    return


def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    table=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-xtab':
            i+=1
            table=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif field=='':
            field=argv[i]
        else:
            if field=='':
                field=argv[i]
            else:
                tiles.append(argv[i])
        i+=1

    if xall==False and len(tiles)==0:
        print('Sorry: there seems to be nothing to do')
        print('-all not set and no tiles to set up provided' )
        return


    if table=='':
        table='MC_tiles.txt'

    kred=os.getenv('KRED')
    if os.path.isfile(table):
        xtab=ascii.read(table)
    elif os.path.isfile(kred+'/config/'+table):
        xtab=ascii.read(kred+'/config/'+table)
    else:
        print('Error: Could not find %s in local director or in kred/config' % table)
        return

    xtab=xtab[xtab['Field']==field]

    if len(xtab)==0:
        print('Could not find field %s in %s' % (field,table))
        return 0

    if xall==False:

        i=0         
        for one_tile in tiles:
            q=xtab[xtab['Tile']==one_tile]
            if i==0:
                xtiles=q.copy()
                i+=1
            else:
                xtiles=vstack([xtiles,q])
    else:
        xtiles=xtab

    if len(xtiles)==0:
        print('Sorry: there seems to be nothing to do')
        print('Looked for the following tiles: ',tiles)
        return

    setup_tiles(xtiles)

    fields=np.unique(xtiles['Field'])
    for one in fields:
        open_log('%s.log' % one)
        foo=xtiles[xtiles['Field']==one]
        for one_tile in foo:
            log_message('SetupTile %s %s' % (one_tile['Field'],one_tile['Tile']))
        close_log()


    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)

