from .AirDensityFromSpecificHumidity import *
from .Dewpoint2SpecificHumidity import *
from .VaporPressure import *
from .evalMatrix import *
from .TransientTotalTimeSeries import *
from .isnotebook import *
from .iswindow import iswindow
from .gethostname import gethostname
from .cftime2datetime import cftime2datetime

def is_djf(month):
    '''Find specific month gith given month

    Example
    -------
    code-block::python
        
        era5_djf = era5.isel(time=pyseb.utils.is_djf(era5['time.year']))

    Returns
    -------
    output: bool in np.ndarray
    '''
    return (month == 12) | (month == 1) | (month == 2)

def args_str2bool(default=True, type=None,help=None, dest=None, option_strings=None):
    '''Return bool depending on input argument...
    '''
    v = default
    import argparse
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
