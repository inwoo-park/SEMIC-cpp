#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import xarray

__all__ = ['save_netcdf']
def save_netcdf(semic, fname:str, time:list): # {{{
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
    '''
    #raise Exception('ERROR: This is under development.')

    nx    = semic.nx # load number of grid.
    ntime = len(time) # time stamp for semic simulation.

    # Initialize Dataset array.
    nc = xarray.Dataset(coords={'time':time,'nx':np.arange(nx)})

    # check requested output
    output_request = semic.output_request
    for varname in output_request:
        if varname in 'smb':
            value = semic.Result.smb.get_value()
        elif varname in 'melt':
            value = semic.Result.melt.get_value()
        elif varname in 'tsurf':
            value = semic.Result.tsurf.get_value()
        elif varname in 'alb':
            value = semic.Result.alb.get_value()
        elif varname in 'alb_snow':
            value = semic.Result.alb_snow.get_value()
        elif varname in 'hsnow':
            value = semic.Result.hsnow.get_value()
        elif varname in 'subl':
            value = semic.Result.subl.get_value()
        elif varname in 'evap':
            value = semic.Result.evap.get_value()
        else:
            raise Exception('ERROR: Given variable name (=%s) is not supported.'%(varname))

        nc[varname] = xarray.DataArray(value, dims=('nx','ntime'))

    print(f'Save {fname}')
    nc.to_netcdf(fname)
    # }}}
