import numpy as np
import pandas as pd
import os

from scipy.signal import convolve
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.integrate import quad


FunY0 = None


def get_funY_approx():
    """ Approximate the integral of gaussian cdf as a piecewise linear function.
    This is solely to speed up the evaluation because the function is called repeatedly in MCMC simulation.

    Returns
    -------
    A linearly-interpolated function
    """
    fn = f'{os.path.dirname(os.path.abspath(__file__))}{os.sep}FunY0_xy.xz'
    if os.path.exists(fn):
        dat = pd.read_pickle(fn)
        k1, k2 = dat['x'], dat['y']
    else:
        k1 = np.concatenate([np.arange(-75, -2), np.arange(-2, 2, 0.0002), np.arange(2, 75)])
        k2 = [quad(norm.cdf, 0, x)[0] for x in k1]
    return lambda x: np.interp(x, k1, k2)


def num_xNorm(wid, fCut, locErrX, locErrY, dilF, MBper, approx: bool = False):
    """ Model based xNorm histogram

    Parameters
    ----------
    wid, fCut, locErrX, locErrY, dilF, MBper: model parameters
    approx: whether to use an approximated version of the integral func, see get_funY_approx
    
    Returns
    -------
    newX : x-bins
    comb: probability mass
    """

    # Parameters
    # binWidth is 2/nBin, x ranges from [-1,1]
    nBin = 100
    binEdges = np.linspace(-1, 1, nBin + 1)
    binCenters = np.mean(np.vstack((binEdges[:-1], binEdges[1:])), axis=0)

    # [Membrane]     calculate the surface density (arc length) in each bin, normalized later
    surf = -np.diff(np.arccos(binEdges))  # Equals to arc length for memb
    
    # [Cytoplasmic]  calculate the area in each bin, normalized later
    # here we used term y to indicate z dimension in the paper.
    y = np.sqrt(1 - binCenters ** 2)  # z-coordinate of each bin for memb

    # Rescale the model cell to match real inner membrane,  r=1 --> R/dilF
    R = wid / 2  # wid = cell width (experimental, from bright-field imaging)
    scaleF = dilF / R  # dilF = factor to reflect that inner membrane is off the cell width 
    
    # fCut is determined from the bottom of the cell ~ 0.3 um
    # fCutHere is relative to the cell midline, with 0 at the middle of cell
    fCutHere = (R - fCut) * scaleF # how much cut in z from z = 0, normalizaed by R, dilF multiplied here and divided later

    # Gaussian blue vector [1,51] in x axis, the locErr was scaled by scaleF
    blurX = norm.pdf(binEdges, 0, locErrX * scaleF / 1000) 

    sigma = locErrY * scaleF / 1000 #loc error in z dimension, but written as y. Assumed to be the same as loc error in x dimension

    # to calculate cdf of the blur function in y axis
    def funY(x):
        return norm.cdf(x, 0, 1)

    # [Membrane Part Blur]
    # arc of z>0 and arc of z<0 for each x. considered blur in z
    surfCut = surf * (funY((fCutHere + y) / sigma) + funY((fCutHere - y) / sigma)) / 2
    surfFinal = convolve(surfCut, blurX)  # apply (gausssian) blur in x now.

    
    # [Cytoplasmic Part Blur]
    cyto = y

    if approx:
        global FunY0
        if FunY0 is None:
            FunY0 = get_funY_approx() # define interpolation function

        re = FunY0((fCutHere + y) / sigma)
        le = FunY0((fCutHere - y) / sigma)

        cytoPart = (re - le) / 2 / y * sigma # change of variables
    else:
        # performe integral for each bin
        cytoPart = np.zeros_like(y)

        for i in range(len(y)//2):
            # note: quad is integral along x at a given x.
            integral_result, _ = quad(funY, (fCutHere - y[i]) / sigma, (fCutHere + y[i]) / sigma)
            cytoPart[i] = integral_result / (2 * y[i]) * sigma
            cytoPart[len(y)-i-1] = cytoPart[i]  # by symmetry

    cytoCut = cyto * cytoPart  # apply fCut with locErrY
    cytoFinal = convolve(cytoCut, blurX) # apply locErr in x axis

    # conv: [1,50]*[1,51]-->[1,100], new binEdges:[-2,2] with 2*nBin    
    conX = np.linspace(binCenters[0] - 1, binCenters[-1] + 1, nBin * 2) #x coordinate after convolution

    # normalizing constants
    surf_sum = np.sum(surf)
    cyto_sum = np.sum(cyto)

    # Combine two parts based on the ratio before applying fCut [surf,cyto]
    comb = (surfFinal / surf_sum * MBper + cytoFinal / cyto_sum * (100 - MBper)) / 100
    newX = conX / dilF # real bin centers after scaling by dilF (membrane/outline)

    return newX, comb


def fitxNorm(x, wid, centers, xNormInput, approx: bool = False):
    """ A function for fitting xNorm data with the model

    Parameters
    ----------
    x: model parameters to estimate
    centers: x-values where the histogram bins are centered at
    wid, xNormInput: xNorm histogram data
    approx: whether to use an approximated version of the integral func, see get_funY_approx

    Returns
    -------
    Root mean-squared error
    """

    fCut, locErrX, dilF, MBper = x
    locErrY = locErrX

    newX, Comb = num_xNorm(wid, fCut, locErrX, locErrY, dilF, MBper, approx=approx)
    interp_func = interp1d(newX, Comb, kind='cubic', fill_value='extrapolate')
    vq = interp_func(centers)
    normF = np.nansum(vq)
    xNormNum = vq / normF

    rmse = np.sqrt(np.sum((xNormNum - xNormInput) ** 2) / len(xNormInput))
    return rmse
