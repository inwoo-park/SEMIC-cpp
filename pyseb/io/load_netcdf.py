#!/usr/bin/env python3
import xarray.backends
from ..utils import iswindow
import netCDF4, xarray

def load_netcdf(fname:str):
    '''Wrap-up loading netcdf file with using netCDF4

    Example
    -------
    code-block:: python

        import pyseb
        nc = pyset.io.load_netcdf('test.nc')

    Parameters
    ----------
    fname: str
        Netcdf filename

    Returns
    -------
    nc: xarray.Dataset
    '''
    if iswindow():
        nc = netCDF4.Dataset(fname)
        return xarray.load_dataset(xarray.backends.NetCDF4DataStore(nc))
    else:
        return xarray.load_dataset(fname)