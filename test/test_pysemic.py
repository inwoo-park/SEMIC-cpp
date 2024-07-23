#!/usr/bin/env python3
import numpy as np
import pytest
import os, sys, platform
import pyseb

@pytest.fixture
def test_load(): # {{{
    nx = 10
    semic = pyseb.libpysemic.SEMIC()
    semic.Initialize(nx)
    semic.Display()
    # print(semic.sf)

    ONES = np.ones((nx,))
    semic.sf = 10*ONES
    # print(semic.sf)
    # }}}

@pytest.mark.skip(reason='Skip loading 2d array')
def test_load2d(): # {{{
    nx = 10
    ntime = 365
    semic = pyseb.libpysemic.SEMIC()
    semic.Initialize(nx)
    semic.Display()
    print(semic.sf)

    ONES = np.ones((nx,ntime))
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

def test_openmp(): # {{{
    semic = pyseb.SEMIC()
    print(f'Without settings...')
    print(f'num threads = {semic.GetOpenmpThreads()}')

    print(f'With setting nproc = 8.')
    semic.num_threads = 8
    semic.SetOpenmpThreads(8)
    print(f'num threads = {semic.GetOpenmpThreads()}')
    # }}}

def test_SemicForcings(): #{{{
    a = 100*np.ones((10,365))
    forcings = pyseb.SemicForcings(10,365)

    forcings.sf[:] = a[:]

    # print(forcings.sf[:])

    # return shape
    print(np.shape(forcings.sf[:]))
    # }}}

def test_DoubleMatrix(): # {{{
    # initialize 2d array.
    b = 10*np.ones((10,2))
    c = 100*np.ones((10,))
    a = pyseb.DoubleMatrix(b)

    # get values.
    print(f'a[:,0] = {a[:,0]}')
    print(f'a[0,:] = {a[0,:]}')
    print(f'a[:] = {a[:]}')
    
    # set values
    a[:,0] = c[:]
    print(a[:,0], a[:,1])

    # replace all values.
    a[:] = b[:]
    # }}}

@pytest.mark.skip(reason='Skip loading 2d array')
def test_SemicForcings_ERA5(): # {{{
    import tqdm
    import xarray, netCDF4
    import scipy.io
    import datetime, pandas

    # move to specific directory
    _dirname = os.path.dirname(__file__)
    os.chdir(_dirname)

    force = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980.mat')
    ntime, nx = force['t2m'].shape
    # nx = 1000 # only 100 nodes are required

    print(f'   shape of t2m = ({nx}, {ntime})')

    semic_f = pyseb.SemicForcings(nx, ntime)
    sf  = force['msr'].T[:nx,:]
    rf  = force['mtpr'].T[:nx,:] - force['msr'].T[:nx,:]
    t2m = force['t2m'].T[:nx,:]
    lwd = force['msdwlwrf'].T[:nx,:]
    swd = force['msdwswrf'].T[:nx,:]
    sp  = force['sp'].T[:nx,:]
    wind= force['wind2m'].T[:nx,:]
    
    semic_f.rf = pyseb.DoubleMatrix(rf)
    semic_f.sf = pyseb.DoubleMatrix(sf)
    semic_f.t2m = pyseb.DoubleMatrix(t2m)
    semic_f.lwd = pyseb.DoubleMatrix(lwd)
    semic_f.swd = pyseb.DoubleMatrix(swd)
    semic_f.sp = pyseb.DoubleMatrix(sp)
    semic_f.wind = pyseb.DoubleMatrix(wind)

    # set air density
    qq = pyseb.utils.Dewpoint2SpecificHumidity(force['d2m'].T[:nx,:], force['sp'].T[:nx,:])
    rhoa = pyseb.utils.AirDensityFromSpecificHumidity(force['sp'].T[:nx,:], force['t2m'].T[:nx,:], qq)
    qq = qq
    rhoa = rhoa
    semic_f.qq  = pyseb.DoubleMatrix(qq)
    semic_f.rhoa= pyseb.DoubleMatrix(rhoa)

    print(f"Show information.")
    print(f"-- Shape of qq = {np.shape(semic_f.qq[:])}")
    # }}}

def test_semic_openmp_ERA5(): # {{{
    '''test run semic with openmp
    '''
    import tqdm
    import xarray, netCDF4
    import scipy.io
    import datetime, pandas

    # move to specific directory
    _dirname = os.path.dirname(__file__)
    os.chdir(_dirname)

    force = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980.mat')
    ntime, nx = force['t2m'].shape
    # nx = 1000 # only 100 nodes are required

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
    qq = qq
    rhoa = rhoa

    # now calculate surface and enery balance
    print(f"   -- Load SEMIC module.")
    semic = pyseb.SEMIC()
    semic.Initialize(nx)
    
    print(f'   -- Prepare initialization and parameters')
    ONES = np.ones((nx,))
    semic.tsurf = (273.15-10)*ONES.copy()
    semic.mask  = 2*np.ones((nx,),dtype=int)
    semic.hsnow = 5*ONES.copy()
    semic.hice  = 10*ONES.copy()
    semic.alb   = 0.8*ONES.copy()
    semic.qmr   = 0.0*ONES.copy()

    # set initial parameter
    semic.Param.amp = 2.0*np.ones((nx,))

    # initialize logger
    df = pandas.DataFrame(columns=['num_thread','time'])
    
    semic.verbose = False

    for num_threads in [1, 2, 3, 4, 5, 6, 8, 10]:
        semic.num_threads = num_threads
        semic.SetOpenmpThreads()
        print(f'   Num threads in SEMIC = {semic.GetOpenmpThreads()}')
        print(f'   Run energy balance!')

        tstart = datetime.datetime.now()
        for _ in range(1):
            for i in range(ntime):
                semic.sf    = sf[:,i]
                semic.rf    = rf[:,i]
                semic.t2m   = t2m[:,i]
                semic.lwd   = lwd[:,i]
                semic.swd   = swd[:,i]
                semic.sp    = sp[:,i]
                semic.wind  = wind[:,i]
                semic.rhoa  = rhoa[:,i]
                semic.qq    = qq[:,i]
                
                semic.RunEnergyAndMassBalance()
        elp_time = datetime.datetime.now()-tstart
        print(f'Elapsed time = {elp_time}')

        df.loc[len(df)] = {'num_thread':num_threads,
                           'time':elp_time}

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(df['num_thread'], df['time'])
    plt.show()
    # }}}

if __name__ == '__main__':
    print('   Do main')
    #test_load()
    # test_load2d()
    #test_load_parameters()
    #test_LongwaveRadiation()
    # test_RunSemic()
    # test_tqdm()
    # test_SemicForcings()
    # test_openmp()
    # test_DoubleMatrix()
    #test_SemicForcings_ERA5()
    test_semic_openmp_ERA5()
