#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

The MCELS images (both for the LMC and SMC) have
some header info that is depracated, notably:

WARNING: FITSFixedWarning: RADECSYS= 'FK5 ' 
the RADECSYS keyword is deprecated, use RADESYSa. [astropy.wcs.wcs]
WARNING: FITSFixedWarning: 'datfix' made the change 'Set MJD-OBS to 51146.000000 from DATE-OBS.
Changed DATE-OBS from '29/11/1998' to '1998-11-29''. [astropy.wcs.wcs]

These generate warnings if the WCS is accessed.

This routine fixes these warnings. The original files 
are copied to a subdirectory old


Command line usage (if any):

    usage: head_fix.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250731 ksl Coding begun

'''


import sys
import os
import shutil
import re
from astropy.io import fits
from astropy.time import Time

def parse_non_iso_date(date_str):
    """Try to convert 'DD/MM/YY' → 'YYYY-MM-DD'"""
    m = re.match(r'(\d{2})/(\d{2})/(\d{2})$', date_str.strip())
    if m:
        day, month, year = m.groups()
        year = int(year)
        year += 1900 if year >= 70 else 2000  # Simple Y2K-ish pivot
        return f"{year:04d}-{int(month):02d}-{int(day):02d}"
    return None

def fix_fits_header(filename):
    updated = False

    # Open the FITS file in update mode
    with fits.open(filename, mode='update') as hdul:
        hdr = hdul[0].header

        # Fix RADECSYS → RADESYS
        if 'RADECSYS' in hdr:
            radecsys = hdr['RADECSYS']
            del hdr['RADECSYS']
            hdr['RADESYS'] = radecsys
            print(f"Replaced RADECSYS with RADESYS = '{radecsys}'")
            updated = True
            if radecsys.strip().upper() == 'FK5' and 'EQUINOX' not in hdr:
                hdr['EQUINOX'] = 2000.0
                print("Added EQUINOX = 2000.0")

        # Fix DATE-OBS
        if 'DATE-OBS' in hdr:
            original = hdr['DATE-OBS'].strip()
            iso_pattern = re.compile(r'\d{4}-\d{2}-\d{2}')
            if iso_pattern.fullmatch(original):
                pass  # Already correct
            else:
                converted = parse_non_iso_date(original)
                if converted:
                    try:
                        t = Time(converted, format='iso', scale='utc')
                        hdr['DATE-OBS'] = t.to_value('iso', subfmt='date')
                        hdr['MJD-OBS'] = t.mjd
                        print(f"Converted DATE-OBS: '{original}' → '{t.value}', MJD-OBS = {t.mjd:.6f}")
                        updated = True
                    except Exception as e:
                        print(f"Warning: Failed to parse converted DATE-OBS='{converted}': {e}")
                else:
                    print(f"Warning: Skipping unrecognized DATE-OBS='{original}'")

    # Backup and save if updated
    if updated:
        backup_dir = 'old'
        os.makedirs(backup_dir, exist_ok=True)
        backup_path = os.path.join(backup_dir, os.path.basename(filename))
        shutil.copy2(filename, backup_path)
        print(f"Backed up original file to: {backup_path}")
        print(f"Header updated: {filename}")
    else:
        print(f"No updates necessary: {filename}")


def steer(argv):
    '''
    Just steer
    '''

    i=1
    files=[]
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        else:
            files.append(argv[i])
        i+=1

    for filename in files:
        if not os.path.isfile(filename):
            print(f"Error: File not found: {filename}")
            sys.exit(1)
        print('filename',filename)
        fix_fits_header(filename)


if __name__ == "__main__":
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
