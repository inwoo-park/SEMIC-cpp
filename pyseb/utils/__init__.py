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
