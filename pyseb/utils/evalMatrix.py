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

    output = (np.mean(((Xobs-Xmean) - (Ymod-Ymean))**2)/Xstd**2 + ((Xmean-Ymean)/Xstd)**2)**(0.5)

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

    Returns
    -------
    E - noramlized root mean square error.
    '''

    if omitnan:
        pos = np.where((~np.isnan(Xobs)) & (~np.isnan(Ymod)))
        Xobs = Xobs[pos]
        Ymod = Ymod[pos]

    #Xstd = np.std(Xobs)
    #Xstd = np.amax(Xobs) - np.amin(Xobs)

    #calculate interquartile range 
    q3, q1 = np.percentile(Xobs, [75 ,25])
    Xstd = q3 - q1

    if Xstd < 0.000001:
        raise Exception('ERROR: Xstd encounters the zero value.')

    return np.sqrt(np.mean((Xobs-Ymod)**2))/Xstd
    # }}}
