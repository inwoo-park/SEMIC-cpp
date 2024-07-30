#!/usr/bin/env python3
import numpy as np
import pytest
import os, sys, platform
import pyseb
import socket
hostname = socket.gethostname().lower().replace('-','')

def InitializeSEMIC(nx): # {{{
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
    semic.verbose = False

    return semic
# }}}

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
    '''test RunSemic with AWS dataset from SEMIC-Krapp2017 repository.
    '''
    import matplotlib.pyplot as plt
    dirname = os.path.dirname(__file__)

    semic = pyseb.SEMIC()
    nx = 1
    ONES = np.ones((nx,))

    # initialize semic module
    semic.Initialize(nx)

    # load dataset
    inputs = np.loadtxt(os.path.join(dirname, '../data/SEMIC-Krapp2017/c01_input.txt'))
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

    print('Load ERA5 interpolate dataset.')
    force = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980.mat')
    ntime, nx = force['t2m'].shape
    nx = 1000 # only 100 nodes are required

    print(f'Shape of t2m = ({nx}, {ntime})')

    print(f'Constrcut SemicForcings class.')
    semic_f = pyseb.SemicForcings()

    print(f'... Okay, now initialize variables.')
    sf  = force['msr'].T[:nx,:]
    rf  = force['mtpr'].T[:nx,:] - force['msr'].T[:nx,:]
    t2m = force['t2m'].T[:nx,:]
    lwd = force['msdwlwrf'].T[:nx,:]
    swd = force['msdwswrf'].T[:nx,:]
    sp  = force['sp'].T[:nx,:]
    wind= force['wind2m'].T[:nx,:]
    
    print('Set rainfall flux')
    semic_f.rf.set_value(rf)
    semic_f.sf.set_value(sf)
    semic_f.t2m.set_value(t2m)
    semic_f.lwd.set_value(lwd) 
    semic_f.swd.set_value(swd)
    semic_f.sp.set_value(sp)
    semic_f.wind.set_value(wind)
    
    print("Calculate specific humidity and air densit.")
    qq = pyseb.utils.Dewpoint2SpecificHumidity(force['d2m'].T[:nx,:], force['sp'].T[:nx,:])
    rhoa = pyseb.utils.AirDensityFromSpecificHumidity(force['sp'].T[:nx,:], force['t2m'].T[:nx,:], qq)

    print('-- Set specific humidity.')
    semic_f.qq.set_value(qq)
    print('-- Set air density.')
    semic_f.rhoa.set_value(rhoa)

    del sf, rf, t2m, lwd, swd, sp, wind
    del rhoa, qq
    del force

    # check variable information.
    print(f'SemicForcings.t2m information')
    print(f'nrow = {semic_f.t2m.nrow}')
    print(f'ncol = {semic_f.t2m.ncol}')
    # }}}

@pytest.mark.skip(reason='Skip Running SEMIC with OpenMP')
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
    nx = 10000 # only 100 nodes are required

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

    # initialize logger
    df1 = pandas.DataFrame(columns=['num_thread','time'])
    df2 = df1.copy(deep=True)
    
    # Initialize num threads depending on system.
    if hostname in ['simba00','simba20']:
        NUM_THREADS = [1, 4, 8, 16, 24]
    else:
        NUM_THREADS = [1, 2, 4, 8]
    
    # Use python and cpp interactive
    for num_threads in NUM_THREADS:
        semic = InitializeSEMIC(nx)
        semic.num_threads = num_threads
        semic.SetOpenmpThreads()
        print(f'   Num threads in SEMIC = {semic.GetOpenmpThreads()}')
        print(f'   Run energy balance!')

        tstart = datetime.datetime.now()
        for _ in range(2):
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

        df1.loc[len(df1)] = {'num_thread':num_threads,
                        'time':elp_time}
    
    # Use cpp only!
    smb = []
    for num_threads in NUM_THREADS:
        semic = InitializeSEMIC(nx)
        semic.num_threads = num_threads
        semic.SetOpenmpThreads()
        print(f'   Num threads in SEMIC = {semic.GetOpenmpThreads()}')

        # Before run semic initialize default value.
        ONES = np.ones((nx,))
        semic.tsurf = (273.15-10)*ONES.copy()
        semic.mask  = 2*np.ones((nx,),dtype=int)
        semic.hsnow = 5*ONES.copy()
        semic.hice  = 10*ONES.copy()
        semic.alb   = 0.8*ONES.copy()
        semic.qmr   = 0.0*ONES.copy()

        print(f'   Run energy balance!')
        # 5-times loop
        tstart = datetime.datetime.now()
        semic.RunEnergyAndMassBalance(f, 2)
        elp_time = datetime.datetime.now()-tstart
        print(f'Elapsed time = {elp_time}')

        df2.loc[len(df2)] = {'num_thread':num_threads,
                        'time':elp_time}
        
        # Okay, check difference between threads
        smb.append(np.array(semic.Result.smb.get_value()))

    # Check final value!
    print(smb[0][0,0], smb[1][0,0])
    assert(np.all(smb[0] == smb[1]))

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(df1['num_thread'], df1['time'], label='C++-Python')
    ax.plot(df2['num_thread'], df2['time'], label='C++only')
    ax.legend()
    plt.show()
    # }}}

@pytest.mark.skip(reason='Skip testing output_request')
def test_OutputRequest(): # {{{
    import tqdm
    import xarray, netCDF4
    import scipy.io
    import datetime, pandas

    # move to specific directory
    _dirname = os.path.dirname(__file__)
    os.chdir(_dirname)

    force = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980.mat')
    ntime, nx = force['t2m'].shape
    nx = 10000 # only 100 nodes are required

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

    # initialize semic!
    semic = InitializeSEMIC(nx)
    semic.num_threads = 4
    semic.SetOpenmpThreads()

    # show default request output
    print(f'Request output = {semic.output_request}')

    # now, modify
    semic.output_request = ['smb']

    # SemicForcing class = f
    # nloop = 2
    semic.RunEnergyAndMassBalance(f, 1)

    try:
        print(semic.Result.hsnow.get_value())
    except:
        print('semic.Result.hsnow is not defined!')
    print(np.shape(semic.Result.smb.get_value()))
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
    # test_SemicForcings_ERA5()
    #test_semic_openmp_ERA5()
    test_OutputRequest()
