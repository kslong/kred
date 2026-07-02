medianrange
===========

.. py:module:: medianrange

.. autoapi-nested-parse::

   Compute median plus error bounds using binomial probabilities
   from eqn 1 of Gott et al. 2001, ApJ, 549, 1

   This computes the percentile range about the median (using a simple
   model fit for N>10 points and a direct summation for N<10) and then
   uses the np.percentile() function to calculate the range.

   R. White, 2023 June 21



Attributes
----------

.. autoapisummary::

   medianrange.rng


Functions
---------

.. autoapisummary::

   medianrange.getbounds
   medianrange.medianrange


Module Contents
---------------

.. py:function:: getbounds(ndata, ncut=10, a=1.2, verbose=False)

   Compute percentile bounds using binomial probabilities.

   Based on equation 1 of Gott et al. 2001, ApJ, 549, 1.
   https://ui.adsabs.harvard.edu/abs/2001ApJ...549....1G/abstract

   Parameters
   ----------
   ndata : int
       Number of input values
   ncut : int, optional
       Use direct calculation for ndata <= ncut points, and model
       fit for ndata > ncut points. Default is 10.
   a : float, optional
       Parameter in model fit. Default is 1.2.
   verbose : bool, optional
       If True, prints the percentile range. Default is False.

   Returns
   -------
   pct_med : float
       Percentile for median of array (always 50.0)
   pct_medlo : float
       Percentile for lower bound of confidence interval
   pct_medhi : float
       Percentile for upper bound of confidence interval
       


.. py:function:: medianrange(data, **kw)

   Compute median and error bounds using binomial probabilities.

   Based on equation 1 of Gott et al. 2001, ApJ, 549, 1.
   https://ui.adsabs.harvard.edu/abs/2001ApJ...549....1G/abstract

   Parameters
   ----------
   data : array-like
       Array of input values
   kw : dict
       Additional keyword arguments passed to getbounds()

   Returns
   -------
   med : float
       Median of array
   medlo : float
       Lower bound of confidence interval
   medhi : float
       Upper bound of confidence interval

   Notes
   -----
   Requires numpy >= 1.23.5 for the percentile method.



.. py:data:: rng

