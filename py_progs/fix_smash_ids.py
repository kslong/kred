#!/usr/bin/env python
"""Fix id and random_id columns from float to string in SMASH FITS tables."""

import gc
import os
import sys
from astropy.table import Table

def _drop_cache(filename):
    """Advise the kernel to evict this file from the page cache."""
    try:
        with open(filename, 'rb') as f:
            os.posix_fadvise(f.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    except (AttributeError, OSError):
        pass

def fix_one(filename):
    t = Table.read(filename)
    changed = []
    for col in ('id', 'random_id'):
        if col in t.colnames and t[col].dtype.kind == 'f':
            t[col] = t[col].astype(str)
            changed.append(col)
    if changed:
        t.write(filename, format='fits', overwrite=True)
        print(f'{filename}: fixed {changed}')
    else:
        print(f'{filename}: nothing to fix')
    del t
    gc.collect()
    _drop_cache(filename)

for f in sys.argv[1:]:
    fix_one(f)
