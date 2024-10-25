#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import xarray
import warnings

__all__ = ['save_netcdf']
def save_netcdf(semic, time:list, fname:str=None, data_type:type=np.float32,
                elements=None, x=None, y=None,
                # output results
                smb=[], melt=[], tsurf=[], alb=[], alb_snow=[], subl=[], evap=[], hsnow=[], swsn=[]): # {{{
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

    elements: list (default: None)

    Ouput
    -----
    nc: xarray.Dataset
        Result array defined in "xarray.Dataset".
    '''

    nx    = semic.nx # load number of grid.
    ntime = len(time) # time stamp for semic simulation.
    #print(f'nx    = {nx}')
    #print(f'ntime = {ntime}')
    if ntime == 1:
        raise Exception('ERROR!')

    # Initialize Dataset array.
    nc = xarray.Dataset(coords={'time':(('time'),np.arange(ntime)),
                                'ncells':(('ncells'),np.arange(nx))})
    #nc['time']   = xarray.DataArray(time, dims='time')
    #nc['ncells'] = xarray.DataArray(np.arange(nx), dims='ncells')

    if np.any(elements):
        nc.coords.update({'elements':(('nelem','ntri'), elements)})
    if np.any(x):
        nc.coords.update({'x':(('ncells'), x)})
    if np.any(y):
        nc.coords.update({'y':(('ncells'), y)})

    # Check given data type
    if not data_type in (np.float32, np.float64, float):
        raise Exception(f"ERROR: Given data type (={data_type}) is not available.")

    # check requested output
    output_request = semic.output_request
    for varname in ['smb','melt','tsurf','alb','alb_snow','subl','evap','swsn']:
        isfind = 0
        if (varname == 'smb') & ((varname in output_request) | np.any(smb)):
            isfind = 1
            if ~np.any(smb):
                value = semic.Result.smb.get_value()
            else:
                value = smb
            attrs = {'units':'m w.e. s-1',
                     'long_name':'surface mass balance'}
        elif (varname == 'melt') & ((varname in output_request) | np.any(melt)):
            isfind = 1
            if ~np.any(melt):
                value = semic.Result.melt.get_value()
            else:
                value = melt
            attrs = {'units':'m w.e. s-1',
                     'long_name':'surface melting'}
        elif (varname == 'tsurf') & ((varname in output_request) | np.any(tsurf)):
            isfind = 1
            if ~np.any(tsurf):
                value = semic.Result.tsurf.get_value()
            else:
                value = tsurf
            attrs = {'units':'K',
                     'long_name':'surface temperature'}
        elif (varname == 'alb') & ((varname in output_request) | np.any(alb)):
            isfind = 1
            if ~np.any(alb):
                value = semic.Result.alb.get_value()
            else:
                value = alb
            attrs = {'units':'-',
                     'long_name':'integrated surface albedo considering snow, ice, and land.'}
        elif (varname == 'alb_snow') & ((varname in output_request) | np.any(alb_snow)):
            isfind = 1
            if ~np.any(alb_snow):
                value = semic.Result.alb_snow.get_value()
            else:
                value = alb_snow
            attrs = {'units':'-',
                     'long_name':'snow albedo'}
        elif (varname == 'hsnow') & ((varname in output_request) | np.any(hsnow)):
            isfind = 1
            if ~np.any(hsnow):
                value = semic.Result.hsnow.get_value()
            else:
                value = hsnow
            attrs = {'units':'m',
                     'long_name':'snow height'}
        elif (varname == 'subl') & ((varname in output_request) | np.any(subl)):
            isfind = 1
            if ~np.any(subl):
                value = semic.Result.subl.get_value()
            else:
                value = subl
            attrs = {'units':'m w.e. s-1',
                     'long_name':'sublimation rate'}
        elif (varname == 'evap') & ((varname in output_request) | np.any(evap)):
            isfind = 1
            if ~np.any(evap):
                value = semic.Result.evap.get_value()
            else:
                value = evap
            attrs = {'units':'m w.e s-1',
                     'long_name':'evapotranspiration'}
        elif (varname == 'swsn') & np.any(swsn):
            isfind = 1
            if ~np.any(swsn):
                value = semic.Result.swsn.get_value()
            else:
                value = swsn
            attrs = {'units':'W m-2',
                     'long_name':'net downward shortwave radiation = (1-alpha)*swd.'}

        if isfind:
            nc[varname] = xarray.DataArray(value.astype(data_type).T, dims=('time','ncells'),
                                           attrs=attrs)

    print(f'Save {fname}')
    if not fname is None:
        nc.to_netcdf(fname)
    else:
        warnings.warn('WARNING: filename for exporting SEMIC result is not defined! Skip saving result in netcdf')

    return nc
    # }}}
