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

to define the area contained by the source.  Various RegTypes are allowed, as detailed
in the routine.  Extra columns are allowed.  Generally, the routines (read and) write out 
astropy tables in 'ascii.fixed_width_two_line' format. 

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


Now if you run::
    
    GetImageFlux.py -autoback whatever.fits masterfile

the routine will constuct a circular region file outside of the source region and use
this for creating a background, and calculate fluxes.  It will also write out a new
master file, with this kind of format::

    No. Source_name    RA     Dec RegType  Major  Minor Theta Color  SourceBack
    --- ----------- -----  ------ ------ ------ -----_ ------ ----------
      1  J0041-7336 10.26 -73.61 ellipse 123.00 110.00    15 yellow     Source
      1  J0041-7336 10.26 -73.61 annulus 171.48 126.00    15 yellow       Back
      1  J0041-7336 10.26 -73.61 annulus 281.00 158.00    15  green     Source
      1  J0041-7336 10.26 -73.61 annulus 353.63 284.00    15  green       Back
      2  J0046-7308 11.67 -73.14 ellipse  92.50  70.00   -55 yellow     Source
      2  J0046-7308 11.67 -73.14 annulus 124.88  95.50   -55 yellow       Back


Some editing may be requred, but then one can convert this into a region file, with the 
routine  :doc:`master2reg.py <api/master2reg/index>` which can be displayed in ds9, and
modified as desired.  Once one is satisfied with the choices of background, one can save
the region file, convert it back to a masterfile, and run::

    GetImageFlux.py whatever.fits new_masterfile

without -autoback.  

The routine :doc:`GetImageFlux.py <api/GetImageFlux/index>` has various options, which can be explored with::

    GetImageFlux.py -h


If one includes the -viz option figures are produced that show the soruce and background regions in a local Figs\_flux directory

The output that is produced has the following format::

    Source_name   Src_flux Src_num_pixels_used Src_mean Src_median Src_mode  Src_min Src_max Back_num_pixels_used Back_mean Back_median Back_mode Back_min Back_max   Net_flux   Med_flux
    ----------- ---------- ------------------- -------- ---------- -------- -------- ------- -------------------- --------- ----------- --------- -------- -------- ---------- ----------
     J0041-7336  515630.51                5041    48.52      39.88    22.60  -196.59  422.03                21992     14.75       11.60      5.30 -1056.70  4373.91  457170.07  142568.93
     J0046-7308 1637613.09                2389   322.97     264.74   148.28  -146.48 4063.22                15016    457.00      258.00   -140.01 -1228.32  4402.29 1021251.81   16100.82
     J0047-7308 2578933.54                1037  1194.03    1138.03  1026.03    58.42 3840.28                 3852    637.29      485.58    182.16 -1462.81  4332.35 2075383.92  676590.43
     J0047-7309 4544013.48                2247  1071.19    1051.90  1013.32   488.10 3687.06                10122    583.66      390.00      2.68 -1042.10  4332.35 3667689.86 1487290.79
     J0048-7133  499304.77               22417    11.00      10.60     9.80  -486.80 3660.66               104763      7.35        7.18      6.84 -2094.51  2970.70  338334.11   76649.03
     J0048-7319 2055453.84                1977   488.06     278.01  -142.09  -536.73 2240.81                 8645    319.37      144.11   -206.40 -1376.09  4399.86 1770546.53  264713.95
     J0049-7314 1843292.78                2824   260.07     231.82   175.33 -1305.88 1644.70                19154    266.97      170.37    -22.84 -1771.55  3942.05 1362180.20  173553.41


which need to be examined in detail.  Src refers to measurements of the source region, while Back of course refers to the Background region. Net\_flux and Med\_flux are two ways of obtaining
the total background subtracted flux.  Net\_flux is based on the median value of the background and is the more accurate.  Med\flux uses the median of the source and background regions 
instead. 
