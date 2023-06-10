#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Calculate backgrounds for images for which
the statistics have been gathered


Command line usage (if any):

    usage: CalcBack.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

230609 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np

from astropy.table import Table,join,vstack
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt



def xeval(xall):
    '''
    xeval returns a metric that descibes now close we actually are to the
    original offsets.  Note that we do not have access to this information
    for real observations; it is only known here because we are using fake
    data.
    '''
    dd=(xall['b']-xall['back'])-(np.average(xall['b']-xall['back']))
    xd=np.sqrt(np.sum(dd*dd))/len(xall)
    return xd



def diff_eval(xdelta,files,return_diff=False):
    '''
    Given a table that contains the measured differences between 
    the flux in different images, and another the table that 
    gives the current version the backgrounds for the images,
    calculated the 'effective chi*2'
    
    Note that this is the only statistic we have access to
    for a real case, since we do not know what the actual offsets
    are.
    '''
    j=0
    xsum=0
    delta=[]
    for one in xdelta:
        diff=one['Delta']-(files['b'][one['j']]-files['b'][one['i']])
        delta.append(diff)
    delta=np.array(delta)
    xsum=np.dot(delta,delta)
    xsum=np.sqrt(xsum)/len(delta)
    if return_diff==True:
        return xsum,delta
    else:
        return xsum




def xpath(xdelta,files,nstart=2):
    '''
    Do a path search to find the offsets
    
    where xdelta is a table containing the differeces
    beween measurments of the "background" in overlapping
    regions, and files is a table listing the files,
    and the offsets for the files.  nstart is
    an index for the file that is used as the start 
    of the search.
    
    The routine returns the files table, with
    the offsets determined
    
    What is returned depends on the order of the
    files in the xdelta table
    
    '''
    xfiles=files.copy()
    xfiles['b']=0.0
    xfiles['Offset_known']=False
    xfiles['Offset_known'][nstart]=True
    # xfiles['b'][nstart]=xfiles['back'][nstart]
    i=0
    while i<len(xfiles):
        known=xfiles[xfiles['Offset_known']==True]
        if len(known)==len(xfiles):
            break
        unknown=xfiles[xfiles['Offset_known']==False]
        xmatch=join(known,xdelta,join_type='left')
        for one_unknown in unknown:
            for one_match in xmatch:
                if one_unknown['i']==one_match['j']:
                    xfiles['b'][one_match['j']]=one_match['b']+one_match['Delta']
                    xfiles['Offset_known'][one_match['j']]=True
                       
        i+=1
    return xfiles
            
        

def monte(xdelta,files):
    i=0
    chi_best=1e10
    xbest=files.copy()
    xchi=[]
    while i<len(files):
        z=xpath(xdelta,files,nstart=i)
        chi=diff_eval(xdelta,z)
        xchi.append(chi)
        # print(chi)
        if chi<chi_best:
            xbest=z.copy()
            chi_best=chi
        i+=1
    print('chi: min %.2f  max %.2f ave %.2f std %.2f' % 
          (np.min(xchi),np.max(xchi),np.average(xchi),np.std(xchi)))
    return xbest
        
def svdskyfit(xcross,threshold=0.1, verbose=True):

    i=xcross['i']
    j=xcross['j']
    npix=xcross['npix']
    delta=xcross['Delta']



    """Solve for sky values using singular value decomposition"""

    i = np.asarray(i)
    j = np.asarray(j)
    npix = np.asarray(npix)
    delta = np.asarray(delta)
    assert len(i)==len(j) and len(npix)==len(j) and len(delta)==len(j)
    assert j.max() == i.max() and j.min() == 0 and i.min() == 0
    n = i.max() + 1

    # compute weights
    # divide weight by mean sum of columns

    wt = np.zeros((n,n),dtype=float)
    wt[i,j] = npix
    wt = wt/(wt.sum(axis=1).mean())

    # stick data from the table into the array

    delta = np.zeros((n,n),dtype=float)
    delta[xcross['i'],xcross['j']] = xcross['Delta']

    a = wt - np.diag(wt.sum(axis=1))
    b = (wt*delta).sum(axis=1)

    if verbose:
        print(f"determinant of A {np.linalg.det(a):.4}")
        print(f"product of diagonal {np.prod(np.diag(a)):.4}")
    assert (a == np.transpose(a)).all()

    # Compute the SVD
    u, w, vt = np.linalg.svd(a, full_matrices=False, compute_uv=True, hermitian=True)
    if verbose:
        print(f"Largest SV {w[0]}")
        print("10 smallest SVs:", w[-10:])

    # Edit the singular values and solve for sky offsets
    winvert = (w>=threshold)/(w + (w < threshold))
    if verbose:
        print(f"{(w<threshold).sum()} SVs removed by threshold {threshold}")
    bkgd = (vt.T @ np.diag(winvert) @ u.T) @ b
    return bkgd


def do_svd(xcross,bfile,threshold=0.1):

    xbfile=bfile.copy()

    background=svdskyfit(xcross,threshold=threshold,verbose=True)

    xbfile['b']=background

    return xbfile


def create_inputs(infile='data/LMC_c45_T07_xxx.txt',xfilter='Ha',exptime=800):
    
    try:
        x=ascii.read(infile)
    except:
        print('Could not read %s' % infile)
        raise IOError
    # print(x.info)
    
    y=x[x['Filter']==xfilter]
    print('Found %d xmatches with %s' % (len(y),xfilter))
    
    y=y[y['Exptime']==exptime]
    
    # print(y.info)
    
    # eliminate any bad values
    
    y['npix'].fill_value=-999
    y=y.filled()
    y=y[y['npix']!=-999]
    
    if len(y)==0:
        print('There are no xmatches to work with')
        raise IOError
    else:
        print('Found %d xmatches with exposure %f' % (len(y),exptime))
    
    x1=np.array(y['file1'])
    x2=np.array(y['file2'])
    both=np.append(x1,x2)
    files=np.unique(both)
    print(len(files))
    
    # Get initial estimate of the background to use
    
    y.sort('npix',reverse=True)
    b=[]
    for one_file in files:
        for one_line in y:
            if one_file==one_line['file1']:
                b.append(-one_line['med1'])
                break
            elif one_file==one_line['file2']:
                b.append(-one_line['med2'])
                break
    # print(len(b))
    
    number=np.linspace(0,len(b)-1,len(b),dtype=int)
    b=np.array(b)

    
    btab=Table([files,number,-b],names=['file','i','b_init'])
    
    btab['b']=0.0
    btab['b_init'].format=btab['b'].format='.3f'
    
    # Now prepare the difference array, which is just a symmetric version
    # of y
    y['Delta']=y['med2']-y['med1']
    y['Delta'].format='.3f'
    
    y1=y['file1','file2','npix','Delta']
    y2=Table([y['file2'],y['file1'],y['npix'],-y['Delta']],
             names=['file1','file2','npix','Delta'])
    y_all=vstack([y1,y2])
    
    # Finally add i and j
    
    bb=btab['file','i']
    bb.rename_column('file','file1')
    
    y_all=join(y_all,bb,join_type='left')
    bb.rename_columns(['file1','i'],['file2','j'])
    y_all=join(y_all,bb,join_type='left')
    
    y_all.sort(['i','j'])
    y_all.sort('npix',reverse=True)
    
    print(np.max(y_all['j']),np.max(y_all['i']))
            
    return btab,y_all


def do_one_tile(field,tile):
    '''
    Geneerate a set of backgrounds to be used for maaching images
    in a tile.
    '''

    infile='Summary/%s_%s_xxx.txt' % (field,tile)
    try:
        x=ascii.read(infile)
    except:
        print('Could not read %s with all of the stats' % infile)
        return 


    xfilt=np.unique(x['Filter'])

    records=[]

    for one_filter in xfilt:
        xx=x[x['Filter']==one_filter]
        xexp=np.unique(xx['Exptime'])
        for one_exp in xexp:
            one_record=[one_filter,one_exp]
            btab,y_all=create_inputs(infile,one_filter,one_exp)
            bfile=infile.replace('xxx.txt','bbb_%s_%d.txt' % (one_filter,one_exp))
            btab.write(bfile,format='ascii.fixed_width_two_line',overwrite=True)
            xfile=infile.replace('xxx.txt','xxx_%s_%d.txt' % (one_filter,one_exp))
            y_all.write(xfile,format='ascii.fixed_width_two_line',overwrite=True)

            print('Field %s  Tile %s Filter %s Exptime %d ' % (field,tile,one_filter,one_exp))

            orig=diff_eval(y_all,btab)

            btab['b']=btab['b_init']
            simple=diff_eval(y_all,btab)

            best_offsets=btab.copy()
            best_val=simple

            xbest=monte(y_all,btab)
            xmonte=diff_eval(y_all,xbest)

            if xmonte<best_val:
                best_offsets=xbest.copy()
                best_val=xmonte


            svd=do_svd(y_all,btab)
            xsvd=diff_eval(y_all,svd)

            if xsvd<best_val:
                best_offset=xsvd.copy()
                best_val=xsvd


            svd=do_svd(y_all,btab,0.01)
            xsvd1=diff_eval(y_all,svd)
            if xsvd1<best_val:
                best_offset=svd.copy()
                best_val=xsvd1


            print('Original b=0:     %.2f' % orig)
            print('Simple  b=b_init:  %.2f' % simple)
            print('Monte           :  %.2f' % xmonte)
            print('SVD             :  %.2f' % xsvd)
            print('SVD1            :  %.2f' % xsvd1)

            one_record.append('%.2f' % orig)
            one_record.append('%.2f' % simple)
            one_record.append('%.2f' % xmonte)
            one_record.append('%.2f' % xsvd)
            one_record.append('%.2f' % xsvd1)


            records.append(one_record)

    records=np.array(records)
    xout=Table(records,names=['Filter','Exptime','None','Simple','Path','SVD','SVD1'])
    xout.write('test.txt',format='ascii.fixed_width_two_line',overwrite=True)

    best_offset.write('b_test.txt',format='ascii.fixed_width_two_line',overwrite=True)


    return



    
    
def steer(argv):
    '''
    This is just a steering routine
    '''
    field=''
    tiles=[]
    xall=False
    redo=True
    # nproc=1
    xback=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-finish':
            redo=False
        elif argv[i]=='sub_back':
            xback=True
        # elif argv[i]=='-np':
        #     i+=1
        #     nproc=int(argv[i])
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
        do_one_tile(field,one)

    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)


