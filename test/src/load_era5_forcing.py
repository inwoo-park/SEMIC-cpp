#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import pyseb
import scipy

def load_era5_forcing(nx:int=None, ntime:int=None):
    '''load era5 forcing variable

    Returns
    ------
    forcing: pyseb.SemicForcings
    '''
    dirname = os.path.dirname(__file__)

    force = scipy.io.loadmat(
            os.path.join(dirname, '../../data/Prepare/ANT_InterpERA5_Day_1980.mat')
            )
    if ntime == None:
        ntime = force['t2m'].shape[0]
    if nx == None: # define number of grid
        nx = force['t2m'].shape[1]

    print(f'   shape of t2m = ({nx}, {ntime})')

    print(f'   Initialize forcing variables.')
    sf  = force['msr'].T[:nx,:]
    rf  = force['mtpr'].T[:nx,:] - force['msr'].T[:nx,:]
    t2m = force['t2m'].T[:nx,:]
    lwd = force['msdwlwrf'].T[:nx,:]
    swd = force['msdwswrf'].T[:nx,:]
    sp  = force['sp'].T[:nx,:]
    wind= force['wind2m'].T[:nx,:]
    
    # set air density
    qq = pyseb.utils.Dewpoint2SpecificHumidity(force['d2m'].T[:nx,:], force['sp'].T[:nx,:])
    rhoa = pyseb.utils.AirDensityFromSpecificHumidity(force['sp'].T[:nx,:], force['t2m'].T[:nx,:], qq)
        
    print(f'   Set SemicForcings')
    f = pyseb.SemicForcings()
    f.nx = nx
    f.ntime = ntime
    f.sf.set_value(sf)
    f.rf.set_value(rf)
    f.sp.set_value(sp)
    f.t2m.set_value(t2m)
    f.lwd.set_value(lwd)
    f.swd.set_value(swd)
    f.wind.set_value(wind)
    f.qq.set_value(qq)
    f.rhoa.set_value(rhoa)

    return f
