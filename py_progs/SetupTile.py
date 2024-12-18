#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Setup the directories and identify the files that are required to carry
out the initial processing steps to run PrepFiles on a set of data


Command line usage (if any):

    usage: SetupTile.py [-h] [-all] [-S7] [-seeing_max 1.3[ [-xtab myfile.txt] [-use_all_data] field tile1 tile2 ...

    where:
        field is that name of a field, e.g LMC_c42
        tile1, tiles are the names of tiles as specified in the MC_tiles.txt file 
    and the optional prameters 
        -h prints this help 
        -all  causes the program to create talbles containing the relevant
            infomation for all of the tiles in a field
        -xtab myfile.txt causes the program to use a different tile defination file
        -S7 causes the S7 chip to be ignored
        -seeing_max causes any images with seing greater than the value given to be ignored.
        -use_all_data causes all of the data regardless of field to be included (as long
            as that field has been Prepped.

Description:  

    This program reads a tile definition file and 
    searches for CCD images that should prepped to create fits files 
    to construct images in tiles of a field.  It writes the results to
    a table in the Summary directory

    The program searches for the tile definition file in the local directory
    and in the $KRED/config directory.

    The program works by reading the det.tab files created by MefSum.py



Primary routines:

    doit

Notes:
    If the directory exists where one is going to place links to the
    files that are to be prepped, sll the file(link)s are removed so
    that only the files including in this run should be in the direcory

    The routine reads tables that are contained in the Summary directory
    to decide what images to include.
                                       
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

    

    print('Looking in %s for CCD images with RA and DEC of %.5f %.5f and size of %.2f deg' % (field,ra_center,dec_center,size_deg))


    dec_min=dec_center-0.5*size_deg
    dec_max=dec_center+0.5*size_deg
    delta_ra=size_deg/np.cos(dec_center/57.29578)
    ra_min=ra_center-0.5*delta_ra
    ra_max=ra_center+0.5*delta_ra
    
    
    x['Dec_max'] = np.max([x['COR1DEC1'],x['COR2DEC1'],x['COR3DEC1'],x['COR4DEC1']],axis=0)
    x['Dec_min'] = np.min([x['COR1DEC1'],x['COR2DEC1'],x['COR3DEC1'],x['COR4DEC1']],axis=0)
    
    x=x[x['Dec_min']<dec_max]

    if len(x)==0:
            print('  Failed to find any images with Dec less than %.5f ' % dec_max)
            return []

    x=x[x['Dec_max']>dec_min]
    if len(x)==0:
            print('  Failed to find any images with Dec greater than %.5f ' % dec_min)
            return []


    
    x['RA_max'] = np.max([x['COR1RA1'],x['COR2RA1'],x['COR3RA1'],x['COR4RA1']],axis=0)
    x['RA_min'] = np.min([x['COR1RA1'],x['COR2RA1'],x['COR3RA1'],x['COR4RA1']],axis=0)
    
    x=x[x['RA_min']<ra_max]
    if len(x)==0:
            print('  Failed to find any images with RA  less than %.5f ' % ra_max)
            return []


    # print(len(x))
    x=x[x['RA_max']>ra_min]
    if len(x)==0:
            print('  Failed to find any images with RA  greater than %.5f ' % ra_min)
            return []


    
    x['Delta_Dec']=x['CENDEC1']-dec_center
    x['Delta_RA']=(x['CENRA1']-ra_center)/np.cos(dec_center/57.29578)
    
    x['Delta_Dec'].format='.2f'
    x['Delta_RA'].format='.2f'

    # print('Found %d CCD extensions that satisfy the input criteria ' % len(x))

    xfiles=[]
    for one in x:
        one_name='%s_%s.fits' % (one['Root'],one['EXTNAME'])
        xfiles.append(one_name)

    x['Filename']=xfiles
    x['Field']=field

    print('  Found %d relevant files in %s' % (len(x),field))
        
    
    return x['Field','Filename','Root','EXTNAME','CENRA1','CENDEC1','COR1RA1','COR1DEC1','COR2RA1','COR2DEC1','COR3RA1','COR3DEC1','COR4RA1','COR4DEC1','Delta_Dec','Delta_RA',
            'FILTER','EXPTIME','MAGZERO','SEEING']



def setup_one_tile(field='LMC_c42',tile='T07',ra=81.108313,dec=-66.177280,size_deg=0.67,s7=True,seeing_max=1000.,use_all_data=False):
    '''
    find the chip images that should be used for a specific tile with a specific size.  

    if s7 is False elimated the S7 chip
    if the seeing max is less than a large value, choose only files that have good seeing
    '''

        
    if use_all_data==True:
        words=field.split('_')
        galaxy=words[0]
        xfiles=glob('Summary/%s*det.tab' % galaxy)
        print('There are %d fields to survey' % len(xfiles))
        fields=[]
        for one_file in xfiles:
            xx=one_file.replace('Summary/','')
            xx=xx.replace('_det.tab','')
            fields.append(xx)
        print('hello Knox ',fields)
        z=[]
        for one in fields:
            one_result=find_files_to_prep(one,tile,ra,dec,size_deg)
            if len(one_result):
                z.append(one_result)
        print('hello Knox ',len(z))
        x=vstack(z)
            
    else:
        x=find_files_to_prep(field,tile,ra,dec,size_deg)

    x.write('foo.txt',format='ascii.fixed_width_two_line',overwrite=True)

    original_length=len(x)
    delta_s7=0
    delta_seeing=0
    if s7==False:
        foo=x[x['EXTNAME']=='S7']
        delta_s7=len(foo)
        x=x[x['EXTNAME']!='S7']
    if seeing_max<100.:
        print('Toasty inside:',len(x),seeing_max)
        foo=x[x['SEEING']>seeing_max]
        delta_seeing=len(foo)
        x=x[x['SEEING']<=seeing_max]
    if len(x)<original_length:
        print('Eliminated %d s7 images and %d bad seeing images of %d that were possible' % (delta_s7,delta_seeing,original_length))
    else:
        print('Using all %d images possible for this tile' % original_length)


    x.meta['comments']=['RA %f' % ra, 'DEC %f' % dec]


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
        one_dir='%s/%s/data' % (DATADIR,one['Field'])
        one_file=one['Filename']
        data_file='%s/%s' % (one_dir,one_file)
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


def setup_tiles(xtab,s7,seeing_max,use_all_data): 
    '''
    Run setup_one_tile multiple times, with the possibility of parallel
    processing

    where 
        if s7 is True, we will include chip s7, but if false we will not
        if seeeing_max is less than some large value this will be used to
            limit thes size
    '''

    for one in xtab:
        try:
            setup_one_tile(one['Field'],one['Tile'],one['RA'],one['Dec'],one['Size'],s7,seeing_max,use_all_data)
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
    seeing_max=1000.
    use_s7=True
    use_all_data=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-S7' or argv[i]=='-s7':
            use_s7=False
        elif argv[i]=='-use_all_data':
            use_all_data=True 
        elif argv[i]=='-xtab':
            i+=1
            table=argv[i]
        elif argv[i]=='-seeing_max':
            i+=1
            seeing_max=eval(argv[i])
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

    setup_tiles(xtiles,use_s7,seeing_max,use_all_data)

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

