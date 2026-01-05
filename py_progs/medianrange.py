"""Compute median plus error bounds using binomial probabilities
from eqn 1 of Gott et al. 2001, ApJ, 549, 1

This computes the percentile range about the median (using a simple
model fit for N>10 points and a direct summation for N<10) and then
uses the np.percentile() function to calculate the range.

R. White, 2023 June 21
"""

import numpy as np
from scipy.special import gammaln

def medianrange(data, **kw):
    """Compute median and error bounds using binomial probabilities.
    
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
    
    """

    data = np.asarray(data)
    ndata = data.size
    pvals = getbounds(ndata, **kw)
    dvals = np.percentile(data, pvals, method="weibull")
    return dvals[0], dvals[1], dvals[2]

# cache to save directly computed values in getbounds()
_bound_cache = {}

def getbounds(ndata, ncut=10, a=1.2, verbose=False):
    """Compute percentile bounds using binomial probabilities.
    
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
        
    """

    percentile = 50.0
    if ndata > ncut:
        prange = percentile*np.sqrt((ndata - a)/ndata**2)
        pvals = np.array([percentile-prange, percentile+prange])
    else:
        try:
            pvals = _bound_cache[ndata]
        except KeyError:
            # calculate probabilities
            confidence = 0.682690 # 1-sigma confidence range
            fpercentile = percentile/100.0
            ii = np.arange(ndata+1, dtype=float)
            p = np.exp(ii*np.log(fpercentile)+(ndata-ii)*np.log(1-fpercentile) +
                       gammaln(ndata+1) - gammaln(ii+1) - gammaln(ndata-ii+1))

            # find percentile bounds that give confidence interval
            ivals = np.interp([(1-confidence)/2, (1+confidence)/2],
                              p.cumsum(), ii)
            pvals = (ivals+1)*100.0/(ndata+1)
            _bound_cache[ndata] = pvals
    if verbose:
        print(f"Percentiles for {ndata}: {pvals}")
    return (percentile, pvals[0], pvals[1])


if __name__ == "__main__":
    # create random test data with 10000 samples
    # normal distribution, sigma=10, mean=5
    # specify a seed to make sample repeatable
    rng = np.random.default_rng(seed=1234579)
    v = rng.normal(loc=5.0, scale=10.0, size=10000)

    # smaller samples with various numbers of points
    vsamples = [v[:3], v[:7], v[:10], v[:15], v[:100], v[:1000], v]

    for v in vsamples:
        med, medlo, medhi = medianrange(v)
        print(f"{len(v):5d} model {med:7.3f} +{medhi-med:7.3f} -{med-medlo:7.3f}")
        # set ncut to high value to force full calculation rather than model fit
        med, medlo, medhi = medianrange(v, ncut=len(v)+1)
        print(f"{len(v):5d} exact {med:7.3f} +{medhi-med:7.3f} -{med-medlo:7.3f}")
        mean = v.mean()
        sigma = v.std()/np.sqrt(len(v))
        print(f"{len(v):5d} mean  {mean:7.3f} +-{sigma:6.3f}")
    print("Simple medians (should be same as above):")
    for v in vsamples:
        print(f"{len(v):5d} {np.median(v):7.3f}")
