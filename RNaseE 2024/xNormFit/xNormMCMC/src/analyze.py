import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

from pymcmcstat import MCMC

from xnormfit.src.fminsearchbnd import fminsearchbnd
from xnormfit.src.num_xNorm import num_xNorm, fitxNorm


DOTTED_LINE = '-' * 100


def xNormFunc(centers, params):
    """
    Estimate an xNorm histogram

    Parameters
    ----------
    centers : histogram centers
    params : args for num_xNorm.num_xNorm func

    Returns
    -------
    Normalized xNorm histogram
    """

    newX, comb = num_xNorm(**params)
    interp_func = interp1d(newX, comb, kind='cubic', fill_value='extrapolate')
    vq = interp_func(centers)
    normF = np.nansum(vq)
    return vq / normF


def run_fmin(wid, edges, xNormInput, lb_, ub_, pnames, use_approx=True):
    """
    Run an fmin search

    Parameters
    ----------
    lb_ : the lower bounds
    ub_ : the upper bounds
    pnames : parameter names
    use_approx: whether to use an approximated version of the integral func, see num_xNorm.get_funY_approx

    Returns
    -------
    Normalized xNorm histogram
    """
    print("")
    print("Executing fmin search")
    x0 = (lb_ + ub_) / 2

    def fun(x):
        centers = (edges[1:] + edges[:-1]) / 2
        return fitxNorm(x, wid, centers, xNormInput, approx=use_approx)

    options = ''
    opt_res, fval, exitfag, output = fminsearchbnd(fun, x0, lb_, ub_, options)

    opt_res_df = pd.Series(opt_res, index=pnames)

    print("")
    print("Estimated fmin parameters")
    from IPython.display import display
    display(opt_res_df)
    return opt_res_df


def mcmc_setup(wid, edges, xNormInput, lb_: list, ub_: list, pnames: list,
               use_approx: bool = True, nsimu: int = 2000, x0=None):
    """
    Set up an MCMC instance using the xNorm information

    Parameters
    ----------
    lb_ : the lower bounds
    ub_ : the upper bounds
    pnames : parameter names
    nsimu : the number of simulations
    x0 : the initial values

    Returns
    -------
    An mcstat instance
    """
    mcstat = MCMC.MCMC()

    # Add data : we add trivial data, only used for determining the size of the histogram with which we calculate the sos
    mcstat.data.add_data_set(np.full(len(xNormInput), 0), np.full(len(xNormInput), 0))

    # the sum of squares function
    def sos_fun(x, data):
        centers = (edges[1:] + edges[:-1]) / 2
        fCut, locErrX, dilF, MBper = x
        params = dict(wid=wid, fCut=fCut, locErrX=locErrX, locErrY=locErrX,
                      dilF=dilF, MBper=MBper, approx=use_approx)
        xNormNum = xNormFunc(centers, params)
        return np.sum((xNormNum - xNormInput) ** 2)

    # Define model settings
    mcstat.model_settings.define_model_settings(sos_function=sos_fun)

    # Add model parameters
    if x0 is None:
        x0 = (lb_ + ub_) / 2

    # Define simulation options
    mcstat.simulation_options.define_simulation_options(nsimu=nsimu, updatesigma=True)

    for i in range(len(pnames)):
        mcstat.parameters.add_model_parameter(
            name=pnames[i],
            theta0=x0[i],  # initial value
            minimum=lb_[i],  # lower limit
            maximum=ub_[i],  # upper limit
        )
    return mcstat


def run_analysis(wid, edges, xNormInput, lb_, ub_, pnames, nsimu: int = 5000, burnin: int = None):
    """
    Perform an MCMC analysis

    Parameters
    ----------
    wid, edges, xNormInput: model data
    lb_ : the lower bounds
    ub_ : the upper bounds
    pnames : parameter names
    nsimu : the number of simulations
    burnin : burn-in period

    Returns
    -------
    An mcstat instance after a simulation is performed
    """

    # run fmin : this is mainly for a comparison
    print(DOTTED_LINE)
    run_fmin(wid, edges, xNormInput, lb_, ub_, pnames)

    # run MCMC
    print(DOTTED_LINE)
    print()
    print("running MCMC")
    mcstat = mcmc_setup(wid, edges, xNormInput, lb_, ub_, pnames, use_approx=True, nsimu=nsimu, x0=None)
    mcstat.run_simulation()

    # Extract results
    results = mcstat.simulation_results.results
    chain = results['chain']
    mcpl = mcstat.mcmcplot
    _ = mcpl.plot_chain_panel(chain, pnames, figsizeinches=(12, 6))
    plt.show()

    if burnin is None:
        # you can choose this parameter based on the previous chain plot
        burnin = input("Burn-in period :")
    burnin = int(burnin)

    print("Estimated MCMC parameters")
    # display chain statistics
    mcstat.chainstats(chain[burnin:, :], results)
    plt.show()

    # display density plots for estimated parameters
    _ = mcpl.plot_density_panel(chain[burnin:, :], pnames, figsizeinches=(12, 6))
    plt.show()
    print(DOTTED_LINE)

    return mcstat
