#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Setup the directories and identify the files that are required to carry
out the initial processing steps to run PrepFiles on a set of data


Command line usage (if any):

    usage: SetupSpecial  -h [-all]  [-field LMC_c40] [-field LMC_c42] source_table.txt  [source_name1] [Source_name]

    where:

        source_table.txt is an astrpy table with the source names and positions

    and the optional prameters 

        -field is that name of a field, e.g LMC_c42.  For multiple files one can add 
        tile1, tiles are the names of tiles as specified in the MC_tiles.txt file 
        If no fields are provided all of the data will be searched.

        -h prints this help 
        -all  causes the program to create 'tiles'  containing the relevant
            infomation for all of the sources in the table

Description:  

    This program reads an astropy table containing at least the following columns

    Source_name  - a string (one_word no space) naming the object
    RA, Dec       - the right ascension and declination of the object in degrees

    and constructs 'tiles' of one or more of the objects from some of all of the
    available data.

    The astropy table can be local, or it can be stored in the config directory

    If the table also contains a column 

    Size   - a size in degrees

    then the that size will be used in selecting the images that can contribute 
    to the total.  If it does not then a default is assumed.

    At present the program does not use the size for anything except what images
    to include in the processing of that tile.


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


def find_files_to_prep(field='LMC_c45',ra_center=81.1,dec_center=-66.17,size_deg=.1):


    
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



def setup_one_tile(field=['LMC_c42'],xmaster='snr',tile='foo',ra=81.108313,dec=-66.177280,size_deg=0.67):
    

    print('Field a',field)

    if len(field)==0:
        det_files=glob('Summary/*det.tab')
        print('Det files ',det_files)
        xfield=[]
        for one in det_files:
            words=one.split('/')
            xfield.append(words[-1].replace('_det.tab',''))
    else:
        xfield=field

    print('Fields to try', xfield)

    i=0
    for one in xfield:
        if i==0:
           x=find_files_to_prep(one,ra,dec,size_deg)
           if len(x)>0:
               print(x)
               i+=1
        else:
            xx=find_files_to_prep(one,ra,dec,size_deg)
            x=vstack([x,xx])
            print(x)
        

    if i>0:
        outfile='Summary/%s_%s.txt' % (xmaster,tile)
        x.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('Failed to set up %s %s  %s' % (field,xmaster,tile))
        return

    tile_dir='%s/%s/%s/' % (PREPDIR,xmaster,tile)

    print('HELLO',tile_dir)
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

        data_dir='%s/%s/data' % (DATADIR,one['Field'])
        print('Toasty ',data_dir)

        if os.path.isdir(data_dir) == False:
            print('Error: %s does not appear to exist' % data_dir)
            print('Run PrepMef on this field first')
            raise IOError


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
        setup_one_tile(one['Field'],one['Tile'],one['RA'],one['Dec'],one['Size'])


def setup_objects(source_names,table_name,fields,xdir):
    '''
    '''


    kred=os.getenv('KRED')
    if os.path.isfile(table_name):
        x=ascii.read(table_name)
    elif os.path.isfile(kred+'/config/'+table_name):
        x=ascii.read(kred+'/config/'+table_name)
    else:
        print('Could not find ',table_name)
    print('got here')
    # x=ascii.read(table_name)

    for one in source_names:
        ztab=x[x['Source_name']==one]
        print('ok',ztab)
        name=ztab['Source_name'][0]
        xra=ztab['RA'][0]
        xdec=ztab['Dec'][0]
        try:
            xsize=ztab['Size'][0]
        except:
            xsize=0.3

        setup_one_tile(field=fields,xmaster=xdir,tile=name,ra=xra,dec=xdec,size_deg=xsize)





def steer(argv):
    '''
    This is just a steering routine

    The command line it interprets is

    SetupSpecial  [-all]  [-field LMC_c40] [-field LMC_c42] source_table source_name1 [Source_name]

    where -all says to do all of the sources in the xtab file
    '''
    field=''
    tiles=[]
    xall=False
    table=''
    xfields=[]
    xobjects=[]
    top_dir=''

    i=1

    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-field':
            i+=1
            xfields.append(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif table=='':
            table= argv[i]
        else:
            xobjects.append(argv[i])
        i+=1

    if table=='':
        print('Sorry: there seems to be nothing to do')
        print('No table of sources provided' )
        return

    if xall==False and len(xobjects)==0:
        print('Sorry: there seems to be nothing to do')
        print('-all not set and no tiles to set up provided' )
        return

    kred=os.getenv('KRED')
    if os.path.isfile(table):
        xtab=ascii.read(table)
    elif os.path.isfile(kred+'/config/'+table):
        xtab=ascii.read(kred+'/config/'+table)
    else:
        print('Error: Could not find %s in local director or in kred/config' % table)
        return




    if xall==True:
        xobjects=xtab['Source_name']
    else:
        # Verify the objects are in the table
        xfound=[]
        objects_in_table=xtab['Source_name']
        for one in xobjects:
            good=False
            igood=[]
            i=0
            for one_object in objects_in_table:
                if one==one_object:
                    xfound.append(one)
                    print('Found %s in table' % one)
                    good=True
                    break
                i+=1
            if good==False:
                print('Could not find %s in table' % one)

            xobjects=xfound

    if len(xobjects)==0:
        print('Sorrry: there seems to be nothing to do')
        return
    else:
        print('There are %d objects to setup' % len(xobjects))

    # At this point the objects we want to analyze are contained in xobjects

    if top_dir=='':
        words=table.split('/')
        top_dir=words[-1].replace('.txt','')



    setup_objects(source_names=xobjects,table_name=table,fields=xfields,xdir=top_dir)

    # fields=np.unique(xtiles['Field'])
    # for one in fields:
    #     open_log('%s.log' % one)
    #     foo=xtiles[xtiles['Field']==one]
    #     for one_tile in foo:
    #         log_message('SetupTile %s %s' % (one_tile['Field'],one_tile['Tile']))
    #     close_log()


    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)

