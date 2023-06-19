#!/usr/bin/env python
# coding: utf-8



import sys
from astropy.io import ascii,fits
import numpy as np
import os
from astropy.table import vstack, join
from glob import glob
import timeit
from multiprocessing import Pool




CWD=os.getcwd()
PREPDIR='%s/DECam_PREP/' % (CWD)

def do_one_file(filein,fileout,b):
    # print(filein,fileout,b)
    f=fits.open(filein)
    f[1].data-=b
    f.writeto(fileout,overwrite=True)
    f.close()
    return


def subtract_one_tile(field='LMC_c42',tile='T07', np=1):
    
    
    table_name='Summary/%s_%s_bbb.txt' % (field,tile)
    try:
        xtab=ascii.read(table_name)
    except:
        print('Failed to find %s' % table_name)
        raise IOError
    
    data_dir='%s/%s/%s/' % (PREPDIR,field,tile)
    out_dir='%s/%s/%s_b/' % (PREPDIR,field,tile)
    if os.path.isdir(data_dir) == False:
        print('Error: %s does not appear to exist' % data_dir)
        raise IOError
    
    if os.path.isdir(out_dir)==False:
        print('Creating directory %s for results' % out_dir)
        os.makedirs(out_dir)
    else:
        print('The directory %s for results already exists' % out_dir)
        
    start_time=timeit.default_timer()
        
    if np<2:
        for one in xtab:
            filein='%s/%s' % (data_dir,one['file'])
            fileout='%s/%s' % (out_dir,one['file'])
            b=one['b']
            do_one_file(filein,fileout,b)
    
    else:
        all_inputs=[]
        for one in xtab:
            filein='%s/%s' % (data_dir,one['file'])
            fileout='%s/%s' % (out_dir,one['file'])
            b=one['b']
            all_inputs.append([filein,fileout,b])
        with Pool(np) as p:
            z=p.starmap(do_one_file,all_inputs)



        
    current_time=timeit.default_timer()
    print('\n***Completely done for field %s tile %s in %.1f s\n' % (field,tile,current_time-start_time))

    
    return


def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True
    nproc=1

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
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

    if len(tiles)==0:
        print('The tiles to be processed must be listed after the field, unless -all is invoked')
        return

    for one in tiles:
        subtract_one_tile(field,one,nproc)

    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)



