#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import scipy.interpolate
import scipy.spatial

def InterpFromMesh2dToGridScipy(elements:list, x:list, y:list, data:list, xg, yg, method:str='linear'):# {{{
    '''InterpFromMesh2dToGridScipy

    Example
    -------
    >>> data = InterpFromMesh2dToGridScipy(md.mesh.elements-1, md.mesh.x, md.mesh.y, md.geometry.surface, xg,yg)

    Parameters
    ----------
    elements: np.ndarray((int,3)) or 2-D array
        Element in triangle meht

    x, y: np.ndarray((float,)) or 1-D array
        Coordinate of triangle meth
    
    data: np.ndarray((float,)) or 1-D array
        Data in triangle

    xg, yg: np.ndarray((float,)) or 2-D array
        x, y coordinates defined in grid. 

    Returns
    -------
    mdata: np.ndarray((float,))
        Interpolated grid data from triangle mesh.
    '''
    # check length of x, y
    if (len(x) != len(y)) | (len(x) != len(data)):
        raise Exception('ERROR: len(x) ~= len(y)')

    # Initialize interpolator
    xg, yg = np.meshgrid(xg,yg)
    points = np.vstack([xg.ravel(), yg.ravel()]).T
    mdata = scipy.interpolate.griddata((x,y), data, points,
                                        method=method) 

    # Fill nan value at outside of elements
    if 0:
        tri = scipy.spatial.Delaunay(np.vstack((x,y)).T)
        tri.simplices = np.array(elements, dtype=int) # update triangle information

        # find points in triangle
        #indices = (tri.find_simplex(points) >= 1)
        indices = tri.find_simplex(points)
    else:
        import matplotlib.tri
        tri = matplotlib.tri.Triangulation(x,y, triangles=elements).get_trifinder()
        indices = (tri(xg.ravel(), yg.ravel()) >= 0)

    mdata[~indices] = np.nan # assign nan value outside of triangles
    mdata = np.reshape(mdata, np.shape(xg))

    return mdata
    # }}}
