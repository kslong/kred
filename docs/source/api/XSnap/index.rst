XSnap
=====

.. py:module:: XSnap

.. autoapi-nested-parse::

   XSnap - Create Snapshots of Astronomical Sources

   Space Telescope Science Institute

   Synopsis
   --------

   Create PNG images (and optionally FITS cutouts) of astronomical sources from
   FITS images, with optional region overlays.

   Command Line Usage
   ------------------

   ::

       XSnap.py [-size arcmin] [-type suffix] [-min vmin] [-max vmax] [-o outname] image.fits master_table
       XSnap.py [-size arcmin] [-type suffix] [-min vmin] [-max vmax] snapshot_table region_table

   Modes
   -----

   The script operates in three modes depending on the inputs:

   1. **Overview mode** (FITS file + master table, no -size):
      Creates a single PNG of the full FITS image with regions from the master
      table overlaid.

      Example::

          XSnap.py -o lmc_ha_overview DECam_SWARP/LMC_c42_T01.ha.fits config/lmc_snr.txt

      Output: ``lmc_ha_overview.png``

   2. **Snapshot mode** (FITS file + master table + -size):
      Creates one PNG snapshot per source in the master table, all extracted from
      the same FITS image. Also creates FITS cutouts in ``xdata/``.

      Example::

          XSnap.py -size 10 -type ha -min -1 -max 20 DECam_SWARP/LMC_c42_T01.ha.fits config/lmc_snr.txt

      Output: ``ximage/{Source_name}.ha.png`` for each source, plus FITS cutouts
      in ``xdata/``

   3. **Multi-file snapshot mode** (snapshot table + region table, no FITS file):
      The snapshot table must contain a ``filename`` column specifying a different
      FITS file for each source. Creates one PNG per row using the corresponding
      FITS file.

      Example::

          XSnap.py -size 10 -type ha snapshots.txt config/lmc_snr.txt

      Where ``snapshots.txt`` contains columns: Source_name, RA, Dec, filename

   **Optional Arguments:**

   -size arcmin
       Size of snapshot cutouts in arcminutes. Required for snapshot modes.

   -type suffix
       Suffix appended to output filenames (e.g., "ha" produces Source.ha.png).
       Useful for distinguishing filter/image types.

   -o outname
       Base name for output file in overview mode (produces outname.png).

   -min vmin
       Minimum value for image scaling. Default: 5th percentile.

   -max vmax
       Maximum value for image scaling. Default: 95th percentile.

   Input Tables
   ------------

   Master/region tables must contain at minimum: Source_name, RA, Dec

   For region overlays, tables may also include:
   - RegType: "circle" or "ellipse"
   - Major, Minor: region sizes in arcseconds
   - Theta: position angle for ellipses

   Output
   ------

   - PNG images are written to ``ximage/`` (snapshots) or current directory (overview)
   - FITS cutouts are written to ``xdata/`` (snapshot modes only)



Functions
---------

.. autoapisummary::

   XSnap.display_fits_image
   XSnap.extract_region
   XSnap.make_many_images
   XSnap.make_many_images2
   XSnap.make_one_image
   XSnap.steer


Module Contents
---------------

.. py:function:: display_fits_image(image_file, scale='linear', ymin=None, ymax=None, invert=True, masterfile='', outfile='')

.. py:function:: extract_region(source_name, ra, dec, size_arcmin, input_fits, outdir='test', default_value=0, frac_off=0.1)

   Extract a region of a given size, but do not write an image if a fraction of an
   image has not data exceeds frac_off




.. py:function:: make_many_images(filename, master, xtype, size, ymin, ymax, frac_off=0.1)

   Create cut-outs of an image, one for each source in a masterfile
   and overlay the regions from the master file on each sanpshot.


.. py:function:: make_many_images2(master, reg, xtype, size, ymin, ymax)

   Create cut-outs of an image, one for each line in a mastefile  
   and overlay the regions from the master file on each sanpshot.


.. py:function:: make_one_image(filename, master, ymin, ymax, outroot='')

   Create a plot of an image and overlay the retions
   from a master file on it.


.. py:function:: steer(argv)

   XSnap.py [-size 10] [-type ha] [-min -1] [-max 20] -out ha [filename or table of snaps]   master_table_of_regions

   without - size we use the full image, and just display eveything
   with a size we make images of each of the source in the master tale




