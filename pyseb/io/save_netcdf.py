#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import xarray
import warnings

__all__ = ['save_netcdf']
def save_netcdf(semic, time:list, fname:str=None, data_type:type=np.float32): # {{{
    '''
    Explain
    -------
    Export results in SEMIC module in netcdf (.nc) format.


    Example
    -------
    >>> import pyseb
    >>> pyseb.io.save_netcdf(semic, 'results.nc')

    Parameters
    ----------
    semic: pyseb.SEMIC
        class SEMIC contains results in format....

    fname: str
        export file name...

    time: list
        datetime array.

    data_type: str (default: np.float32)
        data type for "xarray.DataArray".

    Ouput
    -----
    nc: xarray.Dataset
        Result array defined in "xarray.Dataset".
    '''
    #raise Exception('ERROR: This is under development.')

    nx    = semic.nx # load number of grid.
    ntime = len(time) # time stamp for semic simulation.

    # Initialize Dataset array.
    nc = xarray.Dataset(coords={'time':time,'ncells':np.arange(nx)})

    # Check given data type
    if not data_type in (np.float32, np.float64, float):
        raise Exception(f"ERROR: Given data type (={data_type}) is not available.")

    # check requested output
    output_request = semic.output_request
    for varname in output_request:
        if varname in 'smb':
            value = semic.Result.smb.get_value()
            attrs = {'units':'m w.e. s-1',
                     'long_name':'surface mass balance'}
        elif varname in 'melt':
            value = semic.Result.melt.get_value()
            attrs = {'units':'m w.e. s-1',
                     'long_name':'surface melting'}
        elif varname in 'tsurf':
            value = semic.Result.tsurf.get_value()
            attrs = {'units':'K',
                     'long_name':'surface temperature'}
        elif varname in 'alb':
            value = semic.Result.alb.get_value()
            attrs = {'units':'-',
                     'long_name':'integrated surface albedo considering snow, ice, and land.'}
        elif varname in 'alb_snow':
            value = semic.Result.alb_snow.get_value()
            attrs = {'units':'-',
                     'long_name':'snow albedo'}
        elif varname in 'hsnow':
            value = semic.Result.hsnow.get_value()
            attrs = {'units':'m',
                     'long_name':'snow height'}
        elif varname in 'subl':
            value = semic.Result.subl.get_value()
            attrs = {'units':'m w.e. s-1',
                     'long_name':'sublimation rate'}
        elif varname in 'evap':
            value = semic.Result.evap.get_value()
            attrs = {'units':'m w.e s-1',
                     'long_name':'evapotranspiration'}
        else:
            raise Exception('ERROR: Given variable name (=%s) is not supported.'%(varname))

        nc[varname] = xarray.DataArray(value.astype(data_type).T, dims=('time','ncells'))

    print(f'Save {fname}')
    if not fname is None:
        nc.to_netcdf(fname)
    else:
        warnings.warn('WARNING: filename for exporting SEMIC result is not defined! Skip saving result in netcdf')

    return nc
    # }}}
