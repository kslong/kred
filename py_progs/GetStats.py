#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create a table (or tables) that contains information about 
how overlap regions of various images compare


Command line usage (if any):

    usage: GetStats.py  [-np 8] [-all] field 1 2 3

    where field is a field name, and the following numbers are one 
    or more tiles

    -all implies to do all tiles

    -np 8  implies to perform calculations in parallel, with in
    this case 8 threads

Description:  



Primary routines:


Notes:

    This routine can take a while to run since individual tiles
    involve a large number of files and one must seach for
    essentially  all x-matches.  If there are 100 files
    one needs to look at potentially 100 x 99 choices, and 
    in some cases we have close to 1000 files

    Running it with as many threads as possile is helpful.

    This routine could likely be sped up using an input
    table that identified the corners of a tile in RA and
    Dec.
                                       
History:

230505 ksl Coding begun and routine parallelized

'''

import sys
import os
from glob import glob
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats
from multiprocessing import Pool


def calculate_one(i,files,indir):
    j=i+1
    one=fits.open('%s/%s' % (indir,files[i]))
    records=[]
        
    while j<len(files):
        one_record=[]
        two=fits.open('%s/%s' % (indir,files[j]))
        overlap=np.select([one[0].data*two[0].data!=0],[1],default=0)
        nonzero=np.count_nonzero(overlap)

        if nonzero>0:
            one_record=[files[i],files[j],nonzero]

            # We now need a mask, but we now need the good values 0 and the masked
            # values to be something else, so we subtract 1
            xmask=overlap-1
            one_masked=np.ma.array(one[0].data,mask=xmask)
            two_masked=np.ma.array(two[0].data,mask=xmask)
            mean,median,std=sigma_clipped_stats(one_masked,sigma_lower=3,sigma_upper=2,grow=3)
            one_record=one_record+[mean,median,std]
            mean,median,std=sigma_clipped_stats(two_masked,sigma_lower=3,sigma_upper=2,grow=3)
            one_record=one_record+[mean,median,std]  
            # print(i,j,one_record)
                
            records.append(one_record)

                

        two.close()
        j+=1
    one.close()
    return records


        

def calculate_stats(files,indir,out_name,nproc=1):
    records=[]

    i=0
    all_inputs=[]
    while i<len(files):
        all_inputs.append([i,files,indir])
        i+=1

    if nproc<2:
        i=0
        while i<len(files):
            xrecords=calculate_one(all_inputs[i][0],all_inputs[i][1],all_inputs[i][2])
            records=records+xrecords
            i+=1
    else:
        with Pool(nproc) as p:
            zrecords=p.starmap(calculate_one,all_inputs)
        records=[]
        for one in zrecords:
            records=records+one

    xrecords=np.array(records)
    
    xtab=Table(xrecords,names=['file1','file2','npix','mean1','med1','std1','mean2','med2','std2'])

    colnames=xtab.colnames
    i=0
    while i<len(colnames):
        z=xtab[colnames[i]]
        xxx=[]
        xfloat=True
        for one_value in z:
            try:
                xxx.append(eval(one_value))
            except:
                xfloat=False
        if xfloat==True:
            xtab.remove_column(colnames[i])
            xtab[colnames[i]]=xxx
        i+=1
        

    ztab=xtab[colnames]
    # ztab.info()
    # ztab['mean1'].format='.3f'
    # ztab['mean2'].format='.3f'
    # ztab['med1'].format='.3f'
    # ztab['med2'].format='.3f'
    # ztab['std1'].format='.3f'
    # ztab['std2'].format='.3f'
    ztab.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote %s with %d lines' % (out_name,len(ztab)))
    return ztab
        
def prep_one_tile(field='LMC_c42',tile=7,nproc=1):
    try:
        sumfile='../workspace/Summary/%s_%d_imsum.txt' % (field,tile)
        x=ascii.read(sumfile)
    except:
        print('Could not find %s' % sumfile)
        return



    indir='../workspace/%s_b/%d' % (field,tile)
    if os.path.exists(indir)==False:
        print('Could not find the directory %s with the data' % indir)
        return



    filters=np.unique(x['Filter'])

    for one_filter in filters:
        xx=x[x['Filter']==one_filter]
        print('Starting %s with %d files' % (one_filter,len(xx)))
        outname='../workspace/Summary/%s_%d_%s_stats.txt' % (field,tile,one_filter)
        # outname='%s_%d_%s_stats.txt' % (field,tile,one_filter)
        calculate_stats(xx['Filename'],indir,outname,nproc)



    return


def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    nproc=0

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknwon switch' % argv[i])
            return
        elif field=='':
            field=argv[i]
        else:
            tiles.append(int(argv[i]))
        i+=1

    if xall:
        tiles=[]
        i=1
        while i<17:
            tiles.append(i)
            i+=1

    if len(tiles)==0:
        print('The tiles to be processed must be listed after the field, unless -all is invoked')
    for one in tiles:
        prep_one_tile(field,one,nproc)

    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


