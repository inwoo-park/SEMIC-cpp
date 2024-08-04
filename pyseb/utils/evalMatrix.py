#!/usr/bin/env python3
import numpy as np

__all__ = ['nCRMSE','nRMSE']
def nCRMSE(Xobs, Ymod,omitnan:bool=1): # {{{
    '''nCRMSE - normalized centered root mean square error.

    Usage
    -----
    E = nCRMSE(X, Y)

    Parameters
    ---------
    Xobs - observed variable.
    Ymod - modeled or predicted variable.

    Returns
    -------
    E - noramlized centered root mean square error.
    '''

    if omitnan:
        pos = np.where((~np.isnan(Xobs)) & (~np.isnan(Ymod)))
        Xobs = Xobs[pos]
        Ymod = Ymod[pos]

        if not np.any(pos):
            return 0

    # initialize default variable
    Xstd  = np.std(Xobs)
    Xmean = np.mean(Xobs)
    Ymean = np.mean(Ymod)

    output = np.mean(((Xobs-Xmean) - (Ymod-Ymean))**2)**(0.5)/Xstd

    return output
    # }}}

def nRMSE(Xobs, Ymod, ismethod='std', omitnan:bool=True): # {{{
    '''nRMSE - normalized root mean square error.

    Usage
    -----
    E = nRMSE(Xobs, Ymod)

    Parameters
    ---------
    Xobs - observed variable.
    Ymod - modeled or predicted variable.
    ismethod - which do you want to use for normlize the misfit value? (default: std)
            "std", "quantile", "minmax" are available.

    Returns
    -------
    E - noramlized root mean square error.
    '''

    if omitnan:
        pos = np.where((~np.isnan(Xobs)) & (~np.isnan(Ymod)))
        Xobs = Xobs[pos]
        Ymod = Ymod[pos]

    #calculate interquartile range
    if ismethod == 'quantile':
        q3, q1 = np.percentile(Xobs, [75 ,25])
        Xstd = q3 - q1
    elif ismethod == 'std':
        Xstd = np.std(Xobs)
    elif ismethod == 'minmax':
        Xstd = np.amax(Xobs) - np.amin(Xobs)
    else:
        raise Exception('ERROR: Given method(=%s) is not available. Use "std","quantile", and "minmax" are available.')

    #if Xstd < 0.000001:
    #    raise Exception('ERROR: Xstd encounters the zero value.')

    return np.sqrt(np.mean((Xobs-Ymod)**2))/Xstd
    # }}}
