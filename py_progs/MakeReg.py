#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create region files for fits files in a directory,
assuming that the images have the corners of the images
store in the headers.


Command line usage (if any):

    usage: MakeReg.py xdir     

    where xdir is a directory that contains fits files with
    the corners of the images given

Description:  

    The routine produces a set of region files, one
    for each unique filter and exposure time for
    the files that are in the directory

Primary routines:

    gen_regions

Notes:
                                       
History:

231017 ksl Coding begun

'''





from glob import glob
from astropy.io import fits
from astropy.table import Table
import numpy as np


reg_header='''# Region file format: DS9 version 4.1
global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''


def gen_regions(xdir='.',out_root='test'):
    '''
    Locate the fits files in a directory and
    create region files for each unique
    filter and exposure time.
    '''
    
    files=glob('%s/*fit*' % (xdir))
    if len(files)==0:
        print('Error: no fits files in %s' % xdir)
        return
    
    xfilter=[]
    xtime=[]
    name=[]
    ra=[]
    dec=[]
    d_ra=[]
    d_dec=[]
    for one in files:
        word=one.split('/')
        word=word[-1]
        word=word.split('.')
        xname=word[0]
        
        # try:
        x=fits.open(one)
        h=x[0].header
        qfilter=h['FILTER']
        qtime=h['EXPTIME']
        word=qfilter.split()
        xfilter.append(word[0])
        xtime.append(qtime)
        
        h=x[1].header
        xra=h['CENRA1']
        xdec=h['CENDEC1']
        edge_ra=[h['COR1RA1'],h['COR2RA1'],h['COR3RA1'],h['COR4RA1']]
        edge_de=[h['COR1DEC1'],h['COR2DEC1'],h['COR3DEC1'],h['COR4DEC1']]
        dec_min=np.min(edge_de)
        dec_max=np.max(edge_de)
        ra_min=np.min(edge_ra)
        ra_max=np.max(edge_ra)
            
        delta_dec=(dec_max-dec_min)
        delta_ra=(ra_max-ra_min)*np.cos(xdec/57.29578)
            
        name.append(xname)
        ra.append(xra)
        dec.append(xdec)
        d_ra.append(delta_ra)
        d_dec.append(delta_dec)
            
            
        #except:
        #print('Error: Problem with %s' % one)
    
    xtab=Table([name,ra,dec,d_ra,d_dec,xfilter,xtime],names=['Filename','RA','Dec','Delta_RA','Delta_Dec','Filter','Exptime'])
    
    
    
    filters=np.unique(xtab['Filter'])
    for one_filter in filters:
        xxtab=xtab[xtab['Filter']==one_filter]
        
        times=np.unique(xxtab['Exptime'])
        
        for one_time in times:
            xxxtab=xxtab[xxtab['Exptime']==one_time]
        
            outname='%s_%s_%03d.reg' % (out_root,one_filter,one_time)
        
            f=open(outname,'w')
            f.write(reg_header)
    
    
            for one in xxxtab:
                f.write('box(%.6f,%.6f,%.3f,%.3f,360.0) # text={%s}\n' % (one['RA'],one['Dec'],one['Delta_RA'],one['Delta_Dec'],one['Filename']))
            
            f.close()
    
    return xtab


def steer(argv):
    '''
    Setup the inputs for this routine
    '''

    xdir=argv[1]

    out_root=xdir.replace('/','_')
    if out_root[-1]=='_':
        out_root=out_root[0:len(out_root)-1]


    gen_regions(xdir,out_root) 


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)


