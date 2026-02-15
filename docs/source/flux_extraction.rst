====================================
Fluxes Estimates of extended sources
====================================

Kred contains a tool for estimating fluxes, or more properly, DN, for extended sources
such as SNRs.  The basic idea is that one defines regions for the source and background
and uses the routine :doc:`GetImageFlux.py <api/GetImageFlux/index>` to extract the
fluxes.

The basic input, aside from the fits images for which one wishes to extract the flux,
are 'masterfiles', which minimally must have at least the following columns::

    Source_name         RA        Dec RegType    Major    Minor    Theta      Color
    ----------- ---------- ---------- ------- -------- -------- -------- ----------
     J0041-7336  10.257080 -73.608440 ellipse   123.00   110.00    15.00 yellow
     J0046-7308  11.669170 -73.137470 ellipse    92.50    70.00   -55.00 yellow
     J0047-7308  11.819170 -73.143470 ellipse    55.00    50.00   -45.00 yellow
     J0047-7309  11.902080 -73.155560 ellipse    90.00    60.00   -15.00 yellow
     J0048-7319  12.081670 -73.327670 ellipse    82.50    65.00     0.00 yellow

to define the area contained by the source.  The supported RegTypes are ``ellipse``,
``circle``, and ``annulus``.  For ``circle`` regions, only the ``Major`` column is used
(as the radius) and ``Theta`` is ignored.  Extra columns are allowed.  Generally, the
routines (read and) write out astropy tables in 'ascii.fixed_width_two_line' format.

If you have defined region files by hand on ds9, and saved the region files in a format
similar to::

    # Region file format: DS9 version 4.1
    # Filename: smc_snr_cotton24.txt.reg
    global color=yellow width=3 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    fk5
    ellipse(10.257080,-73.608440,123.00",110.00",15.0) # color=yellow text={J0041-7336}
    ellipse(11.669170,-73.137470,92.50",70.00",-55.0) # color=yellow text={J0046-7308}
    ellipse(11.819170,-73.143470,55.00",50.00",-45.0) # color=yellow text={J0047-7308}
    ellipse(11.902080,-73.155560,90.00",60.00",-15.0) # color=yellow text={J0047-7309}
    ellipse(12.081670,-73.327670,82.50",65.00",0.0) # color=yellow text={J0048-7319}

then you can use the routine  :doc:`reg2master.py <api/reg2master/index>` to comvert the
region file to a master file.

Running GetImageFlux
--------------------

The typical use case is to process multiple images at once with automatic background
generation::

    GetImageFlux.py -auto_back masterfile image1.fits image2.fits ...

or using a wildcard::

    GetImageFlux.py -auto_back smc_snr_cotton24.txt xdata/*fits

The ``-auto_back`` flag constructs a circular annulus background region outside of each
source region and writes out a new region file (``<input>_with_back.txt``) with this
kind of format::

    No. Source_name    RA     Dec RegType  Major  Minor Theta Color  SourceBack
    --- ----------- -----  ------ ------ ------ -----_ ------ ----------
      1  J0041-7336 10.26 -73.61 ellipse 123.00 110.00    15 yellow     Source
      1  J0041-7336 10.26 -73.61 annulus 171.48 126.00    15 yellow       Back
      2  J0046-7308 11.67 -73.14 ellipse  92.50  70.00   -55 yellow     Source
      2  J0046-7308 11.67 -73.14 annulus 124.88  95.50   -55 yellow       Back

Some editing may be requred, but then one can convert this into a region file, with the
routine  :doc:`master2reg.py <api/master2reg/index>` which can be displayed in ds9, and
modified as desired.  Once one is satisfied with the choices of background, one can save
the region file, convert it back to a masterfile, and run::

    GetImageFlux.py masterfile_with_back image1.fits image2.fits ...

without ``-auto_back``.

The routine :doc:`GetImageFlux.py <api/GetImageFlux/index>` has various options, which
can be explored with::

    GetImageFlux.py -h

If one includes the ``-viz`` option, figures are produced that show the source and
background regions in a local ``Figs_Flux/`` directory.


Output
------

When processing multiple images, a single consolidated output file is written::

    Flux_<region_file>.txt

This file contains three types of rows for each source in each image, identified by
the ``SourceBack`` column:

* **Source** -- photometry of the source region
* **Back** -- photometry of the background region
* **Net** -- background-subtracted values

The rows are grouped as Source, Back, Net for each source in each image.  An ``Image``
column identifies which FITS file each measurement came from.

The key columns in the output are:

* ``flux`` -- total flux (DN) in the aperture; for Net rows this is the background-subtracted flux
* ``mean``, ``median`` -- mean and median pixel values; for Net rows these are source minus background
* ``surface_brightness_per_arcsec2`` -- flux per square arcsecond
* ``num_pixels_used`` -- number of unmasked pixels in the aperture
* ``frac_in_image`` -- fraction of the aperture that falls within the image (1.0 = fully contained)

Example output (abbreviated)::

    Source_name SourceBack      flux   mean median frac_in_image                           Image
    ----------- ---------- --------- ------ ------ ------------- -------------------------------
     J0041-7336     Source 515630.51  48.52  39.88          1.00 SMC_c06_T06.ha_sub_r
     J0041-7336       Back 321620.23  14.75  11.60          1.00 SMC_c06_T06.ha_sub_r
     J0041-7336        Net 457170.07  33.77  28.28          1.00 SMC_c06_T06.ha_sub_r

The Net flux is computed as: ``source_flux - num_pixels_used * back_median``.  The Net
mean and median are simply the source value minus the background value.

Sources that fall entirely outside an image are silently skipped.  For sources that are
only partially within an image, ``frac_in_image`` will be less than 1.0, and the
photometry will be based only on the pixels that fall on the image.  Users should
examine this column to decide whether to trust partial measurements.

If no regions in the region file overlap with a given image, that image is skipped
with a message, and it does not contribute to the output file.
