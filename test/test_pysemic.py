#!/usr/bin/env python3
import numpy as np
import pytest
import os, sys, platform
import pyseb
import socket
hostname = socket.gethostname().lower().replace('-','')

# check memory usage
from memory_profiler import profile
import psutil

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
    semic.qmr_res= 0.0*ONES.copy()

    # set initial parameter
    semic.Param.amp = 2.0*np.ones((nx,),dtype=float)
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

@profile
def test_memory_location():# {{{
    nx = int(1e+6)
    s1 = pyseb.SEMIC()
    s2 = pyseb.SEMIC()

    s1.Initialize(nx)
    s2.Initialize(nx)

    print(hex(id(s1)))
    print(hex(id(s2)))

    del s1
    del s2
    # }}}

@profile
def test_SemicForcings(): #{{{
    '''test SemicForcings class.
    '''
    nx = 10000
    ntime = 365
    #a = 100*np.ones((nx, ntime),dtype=float)
    a = np.random.random((nx, ntime))

    if 1: # this process cannot clear the memory
        forcings = pyseb.SemicForcings(nx, ntime)
        forcings.t2m.set_value(1.0)
        forcings.sf.set_value(2.0)
        forcings.rf.set_value(3.0)

        forcings.t2m.set_value(a)
        forcings.wind.set_value(a)

        tmp = forcings.wind.get_value()
        del tmp
        del forcings

    # initialize empty array
    forcings = pyseb.SemicForcings()
    assert(forcings.t2m.nrow == 0)
    assert(forcings.t2m.ncol == 0)
    print('Shape of t2m')
    print(forcings.t2m.nrow, forcings.t2m.ncol)

    # assign variable. 
    forcings.t2m.set_value(a)
    forcings.sf.set_value(a)
    forcings.rf.set_value(a)

    print('Shape of t2m')
    print(forcings.t2m.nrow, forcings.t2m.ncol)

    t2m = forcings.t2m.get_value()
    sf  = forcings.sf.get_value()
    rf  = forcings.rf.get_value()
    assert(np.all(a == t2m))
    del a
    del t2m, rf, sf
    del forcings

    if 0:
        # return shape
        k = forcings.sf.get_value()
        assert(np.shape(k) == (nx, ntime))

        print('clear memory')
        print('   clear a')
        del a
        del k
        print('   clear forcings')
        del forcings
    # }}}

@profile # check memory usage
def test_Semic(): # {{{
    nx = 10000
    ntime = 365*2

    a = np.random.random((nx, ntime))

    # load semic module
    s = pyseb.SEMIC()
    s.Initialize(nx)
    s.output_request = ['smb','melt','alb','tsurf','hsnow']

    s.num_threads = 4
    s.SetOpenmpThreads()

    # set initial variable
    ONES = np.ones((nx,), dtype=float)
    s.tsurf = (273.15-10)*ONES.copy()
    s.hsnow = 5*ONES.copy()
    s.hice  = 10*ONES.copy()
    s.alb   = 0.8*ONES.copy()
    s.qmr   = 0.0*ONES.copy()
    #s.InitializeSemicResult(nx, ntime)

    s.Result.smb.set_value(a)
    s.Result.melt.set_value(a)
    s.Result.alb.set_value(a)
    s.Result.tsurf.set_value(a)
    s.Result.hsnow.set_value(a)

    tmp1 = s.Result.smb.get_value()
    tmp2 = s.Result.melt.get_value()

    print(tmp1[0,0])
    print(tmp2[0,0])

    del tmp1
    del tmp2
    del a
    del s
    # }}}

@profile # check memory usage
def test_DoubleMatrix(): # {{{
    # initialize 2d array.
    nrow = 1000
    ncol = 5000

    b = np.random.random((nrow, ncol))
    a = pyseb.DoubleMatrix()
    a.set_value(b)
    tmp = a.get_value()
    del tmp
    del a
    del b

    # step1
    b = 10*np.ones((nrow, ncol))
    a = pyseb.DoubleMatrix(nrow, ncol)
    a.set_value(b)
    tmp = a.get_value()

    del tmp
    del a
    del b

    # step2
    b = 10*np.ones((nrow, ncol))
    a = pyseb.DoubleMatrix(b)
    tmp = a.get_value()

    del tmp
    del a
    del b

    # step3
    a = pyseb.DoubleMatrix(nrow, ncol)
    a.set_value(1.0)
    tmp = a.get_value()
    del tmp
    del a

    # step4
    b = 10*np.ones((nrow, ncol),dtype=float)
    a = pyseb.DoubleMatrix(b)
    tmp = a.get_value()
    del tmp
    del a
    del b
    # }}}

@profile
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

@profile # chek memory usage
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
        #NUM_THREADS = [1, 4, 8, 16, 24]
        NUM_THREADS = [1, 2, 4]
    else:
        NUM_THREADS = [1, 2, 4, 8]
    
    if 0:
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
        print(semic.output_request)
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
        smb.append(semic.Result.smb.get_value())

    # Check final value!
    print(smb[0][0,0], smb[1][0,0])
    #assert(smb[0][0,0] == smb[1][0,0])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(df1['num_thread'], df1['time'], label='C++-Python')
    ax.plot(df2['num_thread'], df2['time'], label='C++only')
    ax.legend()
    plt.show()

    del semic
    del f
    del smb
    del force
    # }}}

@profile
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
    #nx   = 1 # only 100 nodes are required
    #ntime=10

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
    f = pyseb.SemicForcings(nx, ntime)
    #f.nx = nx
    #f.ntime = ntime
    f.sf.set_value(sf)
    f.rf.set_value(rf)
    f.sp.set_value(sp)
    f.t2m.set_value(t2m)
    f.lwd.set_value(lwd)
    f.swd.set_value(swd)
    f.wind.set_value(wind)
    f.qq.set_value(qq)
    f.rhoa.set_value(rhoa)

    t2m = f.t2m.get_value()
    print(f't2m[0,0] = {t2m[0,0]}')

    # initialize semic!
    semic = InitializeSEMIC(nx)
    semic.num_threads = 4
    semic.SetOpenmpThreads()
    semic.verbose = False

    # show default request output
    print(f'Request output = {semic.output_request}')

    # now, modify
    semic.output_request = ['smb','melt','tsurf']

    # SemicForcing class = f
    # nloop = 2
    print('Run semic!')
    semic.RunEnergyAndMassBalance(f, 1)

    smb = semic.Result.smb.get_value()
    melt = semic.Result.melt.get_value()
    tsurf = semic.Result.tsurf.get_value()
    # smb2 = semic.Result.smb.get_value_list()
    print(smb[0,0])
    print(melt[0,0])
    print(tsurf[0,0])
    print(np.shape(smb))

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # ax.plot(np.arange(365), tsurf[0,:])
    # plt.show()
    #try:
    #    print(semic.Result.hsnow.get_value())
    #except:
    #    print('semic.Result.hsnow is not defined!')
    #smb = semic.Result.smb.get_value()
    #del smb
    #del sf, rf, sp
    #del swd, lwd, wind ,rhoa, qq ,t2m
    del f
    del force
    del semic
    del smb, melt, tsurf
    # }}}

if __name__ == '__main__':
    print('   Do main')
    #test_load()
    # test_load2d()
    # test_load_parameters()
    # test_LongwaveRadiation()
    # test_RunSemic()
    # test_tqdm()
    # test_memory_location()

    if 0:
        proc = psutil.Process(os.getpid())
        mem_before = proc.memory_info().rss / 1024**2
        test_SemicForcings()
        mem_after  = proc.memory_info().rss / 1024**2
        print('Memory usage')
        print(f'   before = {mem_before}')
        print(f'   after  = {mem_after}')

    if 0:
        proc = psutil.Process(os.getpid())
        mem_before = proc.memory_info().rss / 1024**2
        test_Semic()

        mem_after  = proc.memory_info().rss / 1024**2
        print('Memory usage')
        print(f'   before = {mem_before}')
        print(f'   after  = {mem_after}')

    # test_openmp()
    #test_DoubleMatrix()
    if 0:
        proc = psutil.Process(os.getpid())
        mem_before = proc.memory_info().rss / 1024**2
        test_SemicForcings_ERA5()

        mem_after  = proc.memory_info().rss / 1024**2
        print('Memory usage')
        print(f'   before = {mem_before}')
        print(f'   after  = {mem_after}')
    if 1:
        proc = psutil.Process(os.getpid())
        mem_before = proc.memory_info().rss / 1024**2
        test_semic_openmp_ERA5()

        mem_after  = proc.memory_info().rss / 1024**2
        print('Memory usage')
        print(f'   before = {mem_before}')
        print(f'   after  = {mem_after}')

    if 0:
        proc = psutil.Process(os.getpid())
        mem_before = proc.memory_info().rss / 1024**2
        test_OutputRequest()

        mem_after  = proc.memory_info().rss / 1024**2
        print('Memory usage')
        print(f'   before = {mem_before}')
        print(f'   after  = {mem_after}')
