#!/usr/bin/env python3
import numpy as np
from ..classes.Model import Model, Mesh2d
import netCDF4

__all__ = ['loadmodel_netcdf']
def loadmodel_netcdf(fname:str): # {{{
    '''return simple model written in ISSM model.

    Parameters:
    -----------
    fname: str - input file names
    '''

    # load netcdf file
    nc = netCDF4.Dataset(fname)

    # just load mesh information
    md = Model()
    setattr(md.mesh, 'x', nc['/mesh/x'][:])
    setattr(md.mesh, 'y', nc['/mesh/y'][:])
    setattr(md.mesh, 'elements', nc['/mesh/elements'][:])
    setattr(md.mesh, 'lat', nc['/mesh/lat'][:])
    setattr(md.mesh, 'long', nc['/mesh/long'][:])

    setattr(md.mesh, 'numberofvertices', nc['/mesh/numberofvertices'][:])
    setattr(md.mesh, 'numberofelements', nc['/mesh/numberofelements'][:])

    setattr(md.mesh, 'vertexonboundary', nc['/mesh/vertexonboundary'][:])

    setattr(md.mesh, 'extractedvertices', nc['/mesh/extractedvertices'][:])
    setattr(md.mesh, 'extractedelements', nc['/mesh/extractedelements'][:])    

    return md
    # }}}
