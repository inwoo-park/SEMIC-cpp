#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import pyseb

def test_load(): # {{{
    nx = 10
    semic = pyseb.libpysemic.SEMIC()
    semic.Initialize(nx)
    semic.Display()
    print(semic.sf)

    ONES = np.ones((nx,))
    semic.sf = 10*ONES
    print(semic.sf)
    # }}}

def test_load_parameters(): # {{{
    semic = pyseb.SEMIC()
    param = semic.Param
    const = semic.Const
    print(semic.Param.ceff)
    print(semic.Const.clv)

    print('Show Const classs')
    semic.Const.Display()
    # }}}

def test_LongwaveRadiation(): # {{{
    nx = 10
    semic = pyseb.SEMIC()
    semic.Initialize(nx)

    ONES = np.ones((nx,))
    tsurf = 273.15*ONES
    for i in range(nx):
        tsurf[i] = tsurf[i] - i
    semic.tsurf = tsurf

    semic.LongwaveRadiationUp()
    print(semic.lwup)
    # }}}

def test_RunSemic(): # {{{
    import matplotlib.pyplot as plt

    semic = pyseb.SEMIC()
    nx = 1
    ONES = np.ones((nx,))

    # initialize semic module
    semic.Initialize(nx)

    # load dataset
    inputs = np.loadtxt('../data/c01_input.txt')
    sf  = inputs[:,0]
    rf  = inputs[:,1]
    swd = inputs[:,2]
    lwd = inputs[:,3]
    wind = inputs[:,4]
    sp   = inputs[:,5]
    rhoa = inputs[:,6]
    qq   = inputs[:,7]
    t2m  = inputs[:,8]

    semic.mask = [1]
    semic.verbose = 0 

    # initialize albedo
    semic.alb = 0.85*ONES[:]
    semic.tsurf = (273.15-10)*ONES[:]

    # initial parameters
    semic.Param.amp = 3*ONES[:]

    ntime = len(t2m)
    print(f'   t2m: {semic.t2m}')
    print(f'   ntime: {ntime}')

    # initialize results array
    smb = np.zeros((nx,ntime))
    melt= np.zeros((nx,ntime))

    nloop = 20
    for loop in range(nloop):
        print(f'loop = {loop}')
        for i in range(365):
            semic.sf = [sf[i]]
            semic.rf = [rf[i]]
            semic.swd = [swd[i]]
            semic.lwd = [lwd[i]]
            semic.wind = [wind[i]]
            semic.sp = [sp[i]]
            semic.rhoa = [rhoa[i]]
            semic.qq = [qq[i]]
            semic.t2m = [t2m[i]]
            semic.RunEnergyAndMassBalance()
            if loop == nloop-1:
                smb[:,i] = semic.smb
                melt[:,i] = semic.melt
    print('done')

    fig, ax = plt.subplots()
    ax.plot(np.arange(ntime), smb.ravel())
    plt.show()
    # }}}

def test_tqdm(): # {{{
    '''check tqdm show progress bar
    '''
    import tqdm

    istqdm = 1

    # disable option for not showing tqdm.
    for i in tqdm.tqdm(range(10), disable=(not istqdm)):
        print(i)
    # }}}

if __name__ == '__main__':
    print('   Do main')
    #test_load()
    #test_load_parameters()
    #test_LongwaveRadiation()
    # test_RunSemic()
    test_tqdm()