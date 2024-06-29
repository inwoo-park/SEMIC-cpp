#!/usr/bin/env python
# coding: utf-8

# # Explains
# 
# * To convert ipynb file to python (*.py) and run optimization
# ```
# jupyter-nbconvert --to script run_OptimizationERA5.ipynb
# python run_OptimizationERA5.py
# ```

# # Load modules

# In[ ]:


def isnotebook():
   try:
      from IPython import get_ipython
      if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
         return False
   except ImportError:
      return False
   except AttributeError:
      return False
   return True
if isnotebook():
   get_ipython().run_line_magic('matplotlib', 'inline')
   get_ipython().run_line_magic('load_ext', 'autoreload')
   get_ipython().run_line_magic('autoreload', '2')
import os, sys
import numpy as np
import pandas, json
import matplotlib.pyplot as plt
sys.path.append('../install/lib')
sys.path.append('../src')
import pyseb
from pyseb.utils import Dewpoint2SpecificHumidity, VaporPressure, AirDensityFromSpecificHumidity
from pyseb.utils.evalMatrix import *

import copy
from issmmodules import *
import scipy
import scipy.io
import tqdm
import datetime # checking elapsed time.

# parallel computing
import multiprocessing

from interpRACMO23p2_ERA5 import *
from interpXarrayGridToMesh import *
import xarray
import shutil
import argparse


# # Initialize arguments

# In[ ]:


parser = argparse.ArgumentParser(prog='Optimization',
                            description='Estimate the SEMIC parameters with PSO method')
parser.add_argument('-ngen',type=int,
                   help='Set number of generation',
                   default=2)
parser.add_argument('-npop',type=int,
                   help='Set number of population',
                   default=10)
parser.add_argument('-ncpu',default=30,type=int,
                   help='Set set number of cpu for parallel computing',)
parser.add_argument('-omega',default=0.5, type=float,
                   help='Inertia weight (high -> global optial, low -> local optimal)',
                   )
parser.add_argument('-f',
                   help='dummy part for ipython')
parser.add_argument('-datadir',type=str,
                   default='./data/PSO')
parser.add_argument('-debug',type=int,
                   default=0)
parser.add_argument('-freq',type=str,default='day',
                   help='frequnecy of forcing variable. "day" or "mon" is available.')
parser.add_argument('-amp_max',type=float, default=5,
                    help='Available maximum temperature amplitude for "monthly" experiment". (default: 5)',
                    )
parser.add_argument('-tqdm',type=int, default=0,
                   help='Show tqdm progress bar. (default: 0)')
parser.add_argument('-opt_method',type=int,default=1,
                   help='''optimization method for monthly dataset.
                   1: use homogeneous melting
                   2: use heterogeneous temperature amplitude at each basin.
                   3: use some multiplication at each baseh''')
args = parser.parse_args()

if isnotebook():
    args.freq      = 'mon'
    args.ngen  = 2
    args.npop  = 10
    args.nloop = 2
    args.opt_method=2
    args.debug = 1 # force to change
    args.amp_max  = 10.
    args.tqdm = True

print(f'Information with given argument.')
print(f'datadir: {args.datadir}')
print(f'npop:    {args.npop}')
print(f'ngen:    {args.ngen}')
print(f'freq:    {args.freq}')
print(f'-- amp_max: {args.amp_max}')
print(f'debug:   {args.debug}')
print(f'opt_method: {args.opt_method}')
print(f'===== system information')
print(f'number of cpu: {args.ncpu}')

if os.path.isdir(args.datadir):
    print('   Remove existing directory')
    shutil.rmtree(args.datadir)
os.mkdir(args.datadir)


# # Load ANT model

# In[ ]:


print(f'Load ANT model')
md = loadmodel_netcdf('../data/ANT_Mesh.nc')
racmo_melt, racmo_time = interpRACMO23p2_ERA5(md.mesh.lat, md.mesh.long, 'snowmelt',
                          'time',[datetime.datetime(1980,1,1), datetime.datetime(2001,1,1)])
print(f'Load ANT: extract melting region.')
md = md.extract(racmo_melt.mean(axis=1) > 0.)
print(f'   number of vertices: {md.mesh.numberofvertices}')


# In[ ]:


print(f'IMBIE2: interpolate imbie mask')
nc = netCDF4.Dataset('../data/IMBIE2/basinNumbers_8km.nc')
mask_imbie = interpXarrayGridToMesh('../data/IMBIE2/basinNumbers_8km.nc', md.mesh.x, md.mesh.y,
                             xname='x',yname='y',varname='basinNumber',method='nearest')
mask_imbie = np.array(mask_imbie.values, dtype=int) # extract only values!
if isnotebook():
    print(f'IMBIE2: {np.unique(mask_imbie)}')
    plotmodel(md,'data',mask_imbie)


# # Load RACMO23p2-ERA5

# In[ ]:


# m/yr -> water m/sec
rho_ice = 917
rho_freshwater = 1000;
yts = 365*24*3600
isunit = rho_ice/rho_freshwater*yts
racmo_smb, racmo_time = interpRACMO23p2_ERA5(md.mesh.lat, md.mesh.long, 'snowmelt',
                          'time',[datetime.datetime(1980,1,1), datetime.datetime(2001,1,1)])
racmo_smb = racmo_smb * isunit

racmo_melt, racmo_time = interpRACMO23p2_ERA5(md.mesh.lat, md.mesh.long, 'snowmelt',
                          'time',[datetime.datetime(1980,1,1), datetime.datetime(2001,1,1)])
racmo_melt = racmo_melt * isunit

# surface temperature (unit: K)
racmo_tsurf, racmo_time = interpRACMO23p2_ERA5(md.mesh.lat, md.mesh.long, 'tskin',
                          'time',[datetime.datetime(1980,1,1), datetime.datetime(2001,1,1)])

# net shortwave radiation (unit: K)
racmo_swsn, racmo_time = interpRACMO23p2_ERA5(md.mesh.lat, md.mesh.long, 'swsn',
                          'time',[datetime.datetime(1980,1,1), datetime.datetime(2001,1,1)],
                          'use_cftime',1)

print('Generate Xarray')
dsracmo = xarray.Dataset(data_vars={'smb':(['nx','time'], racmo_smb),
                                   'melt':(['nx','time'], racmo_melt),
                                   'swsn':(['nx','time'], racmo_swsn),
                                   'tsurf':(['nx','time'], racmo_tsurf)},
                        coords={'x':('nx',md.mesh.x),
                               'y':('nx',md.mesh.y),
                               'time':('time',racmo_time)})
del racmo_smb, racmo_melt, racmo_tsurf, racmo_swsn


# # Load ERA5 dataset

# In[ ]:


if 0: # interpoalte! dataset
    data_era5 = {}
    for varname in ['d2m','t2m','msr','mtpr','msdwlwrf','msdwswrf','sp','wind2m']:
        print(f'Processing... {varname}')
        ds = interpCDO(md.mesh.lat, md.mesh.long,
                       f'../data/ERA5/day/era5_{varname}_1980,2000.nc', varname=varname)
        data_era5[varname] = ds[varname].values
    
    print('Save dataset in matlab format.')
    scipy.io.savemat('../data/Prepare/ANT_InterpERA5.mat',data_era5)
else:
    if args.freq == 'day':
        print('Load saved dataset.')
        data_era5 = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980,2000.mat')
    elif args.freq == 'mon':
        data_era5 = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Month_1980,2000.mat')
    else:
        raise Exception(f'ERROR: Given freq(={args.freq}) is not available.')


# # Initialize forcing variables of ERA5 for SEMIC

# In[ ]:


rho_freshwater = 1000 # kg m-3

# extract ice shelf region
# pos_iceshelf = np.where(md.mask.ocean_levelset < 0)[0]
# pos_iceshelf = np.where(np.ones((md.mesh.numberofvertices,)))[0]
pos_iceshelf = md.mesh.extractedvertices-1

print('ERA5: Set forcing variable for SEMIC.')
sf = data_era5['msr'].T/rho_freshwater
rf = (data_era5['mtpr'].T - data_era5['msr'].T)/rho_freshwater

t2m = data_era5['t2m'].T
sp  = data_era5['sp'].T # surface pressure
lwd = data_era5['msdwlwrf'].T
swd = data_era5['msdwswrf'].T

wind = data_era5['wind2m'].T

# estimate the air density and specific humidity
print('ERA5: Prepare specific humidity and air density')
qs = Dewpoint2SpecificHumidity(data_era5['d2m'], data_era5['sp'])
rho_air = AirDensityFromSpecificHumidity(data_era5['sp'], data_era5['t2m'], q=qs)
qq = qs.T
rhoa = rho_air.T

print('ERA5: Extract ice shelf values.')
force = {}
force['sf']  = sf[pos_iceshelf,:]
force['rf']  = rf[pos_iceshelf,:]
force['t2m'] = t2m[pos_iceshelf,:]
force['sp']  = sp[pos_iceshelf,:]
force['lwd'] = lwd[pos_iceshelf,:]
force['swd'] = swd[pos_iceshelf,:]
force['wind']= wind[pos_iceshelf,:]
force['qq']  = qq[pos_iceshelf,:]
force['rhoa']= rhoa[pos_iceshelf,:]
nx    = np.shape(qq)[0]
ntime = np.shape(qq)[1]

print(f'ERA5: Available maximum size of time: {ntime}')


# # Prepare PSO optimization

# In[ ]:


import operator
import math

from deap import base
from deap import creator, tools
import random
random.seed(100)


# In[ ]:


# optimization with parameters
#              min, max, smin, smax
opts = {}
if args.freq == 'day':
    opts['amp']      = [1, 5]
    opts['hcrit']    = [np.log10(0.001), np.log10(1)]
    opts['rcrit']    = [0, 1]
    opts['mcrit']    = [np.log10(1e-9), np.log10(1e-7)]
    opts['alb_smax'] = [0.7, 0.95]
    opts['alb_smin'] = [0.2, 0.6]
    opts['albi']     = [0.2, 0.6]
elif args.freq == 'mon':
    if args.opt_method == 1:
        opts['amp']      = [1, 5]
        opts['hcrit']    = [np.log10(0.001), np.log10(1)]
        opts['rcrit']    = [0, 1]
        opts['mcrit']    = [np.log10(1e-9), np.log10(1e-7)]
        opts['alb_smax'] = [0.7, 0.95]
        opts['alb_smin'] = [0.2, 0.6]
        opts['albi']     = [0.2, 0.6]
    elif args.opt_method == 2:
        for nb in range(16):
            opts['amp%d'%(nb)] = [1, args.amp_max]
        opts['hcrit']    = [np.log10(0.001), np.log10(1)]
        opts['rcrit']    = [0, 1]
        opts['mcrit']    = [np.log10(1e-9), np.log10(1e-7)]
        opts['alb_smax'] = [0.7, 0.95]
        opts['alb_smin'] = [0.2, 0.6]
        opts['albi']     = [0.2, 0.6]
    else:
        raise Exception("ERROR: Given option is not available.")
else:
    raise Exception('ERROR: Not supported yet.')


# # Initialize particle information
# 
# ## SEMIC-day experiment
# 
# ## SEMIC-mon experiment
# 
# * Amplitude of maximum temperature should be larger than previous one (e.g., 5 K).

# In[ ]:


creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Particle", dict, fitness=creator.FitnessMin, speed=list, 
    smin=None, smax=None, pmin=None, pmax=None, best=None)


# In[ ]:


def generate(opts_bnd):
    _parts = {}
    smin   = []
    smax   = []
    # maximum and minimum position for each particle.
    pmax   = []
    pmin   = []
    for key, value in opts_bnd.items():
        _parts[key] = random.uniform(value[0], value[1])
        smin.append(-0.2*(value[1] - value[0]))
        smax.append( 0.2*(value[1] - value[0]))
        pmax.append(value[1])
        pmin.append(value[0])
    # _parts = [random.uniform(_pmin, _pmax) for _pmin, _pmax in zip(pmin, pmax)]
    part = creator.Particle(_parts) 
    part.speed = [random.uniform(_smin, _smax) for _smin, _smax in zip(smin, smax)]
    part.smin = smin
    part.smax = smax
    part.pmin = pmin
    part.pmax = pmax
    
    return part


# In[ ]:


def updateParticle(part, best, phi1=2, phi2=2, omega=args.omega):
    '''update particle location
    '''
    npts = len(part.speed)
    loc             = []
    loc_local_best  = []
    loc_global_best = []
    for key in part.keys():
        loc.append(part[key])
        loc_local_best.append(part.best[key])
        loc_global_best.append(best[key])
    
    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = map(operator.mul, u1, map(operator.sub, loc_local_best,  loc))
    v_u2 = map(operator.mul, u2, map(operator.sub, loc_global_best, loc))
    part.speed = list(map(operator.add, list(map(operator.mul, [omega]*npts, part.speed)),
                          map(operator.add, v_u1, v_u2)))

    # constrain its maximum and minimum speed.
    for i, speed in enumerate(part.speed):
        if abs(speed) < part.smin[i]:
            part.speed[i] = math.copysign(part.smin[i], speed)
        elif abs(speed) > part.smax[i]:
            part.speed[i] = math.copysign(part.smax[i], speed)
    loc = list(map(operator.add, loc, part.speed))

    # constrain its max. and min. position
    for idx, key, value in zip(range(len(part)), part.keys(), part.values()):
        if part[key] > part.pmax[idx]:
            part[key] = part.pmax[idx]
        if part[key] < part.pmin[idx]:
            part[key] = part.pmin[idx]

    # update dict
    for idx, key in enumerate(part.keys()):
        part[key] = loc[idx]


# In[ ]:


def evaluateSEMIC(part, md, mask_imbie, force, nx=12744, ntime=365, nloop=2):
    '''
    opt_amp      = [1, 5]
    opt_hcrit    = [np.log10(0.001), np.log10(1)]
    opt_rcrit    = [0, 1]
    opt_mcrit    = [np.log10(1e-9), np.log10(1e-7)]
    opt_alb_smax = [0.7, 0.95]
    opt_alb_smin = [0.2, 0.6]
    opt_albi     = [0.2, 0.6]
    '''
    # print('Initialize structure')
    semic = pyseb.SEMIC()
    
    semic.Initialize(nx)
    ONES = np.ones((nx,),dtype=float)
    semic.mask = 2*np.ones((nx,),dtype=int)
    semic.verbose = 0
    semic.alb_scheme = 3 # use isba albedo scheme
    
    semic.alb = 0.8*ONES[:]
    semic.alb_snow = 0.8*ONES[:]
    semic.tsurf = 260.*ONES[:]
    semic.hsnow = 5*ONES[:]
    semic.hice  = 0*ONES[:]
    
    # update parameter
    semic.Param.ceff = 2e+6
    semic.Param.csh  = 2e-3
    semic.Param.clh  = 5e-4
    if args.freq == 'day':
        semic.Param.amp  = part['amp']*ONES[:]
    elif args.freq == 'mon':
        if args.opt_method == 1:
            semic.Param.amp = part['amp']*ONES[:]
        elif args.opt_method == 2:
            amp = ONES[:]
            for nb in range(16):
                pos = (mask_imbie == nb)
                amp[pos] = part['amp%d'%(nb)]
            semic.Param.amp = amp[:]
    semic.Param.alb_smax = part['alb_smax']
    semic.Param.alb_smin = part['alb_smin']
    semic.Param.albi = part['albi']
    semic.Param.albl = 0.15
    semic.Param.hcrit = 10**part['hcrit']
    semic.Param.rcrit = part['rcrit']
    semic.Param.mcrit = 10**part['mcrit']
    
    # initialize results array
    # ntime = 365
    smb     = np.zeros((nx,ntime))
    # smb_snow= np.zeros((nx,ntime))
    melt    = np.zeros((nx,ntime))
    tsurf   = np.zeros((nx,ntime))
    alb     = np.zeros((nx,ntime))
    netswd  = np.zeros((nx,ntime))
    shf     = np.zeros((nx,ntime))
    lhf     = np.zeros((nx,ntime))

    for loop in range(nloop):
        for i in range(ntime):
            semic.sf   = force['sf'][:,i]
            semic.rf   = force['rf'][:,i]
            semic.swd  = force['swd'][:,i]
            semic.lwd  = force['lwd'][:,i]
            semic.wind = force['wind'][:,i]
            semic.sp   = force['sp'][:,i]
            semic.rhoa = force['rhoa'][:,i]
            semic.qq   = force['qq'][:,i]
            semic.t2m  = force['t2m'][:,i]
            semic.RunEnergyAndMassBalance()
            # break
            if loop == nloop-1:
                smb[:,i]   = semic.smb
                melt[:,i]  = semic.melt
                tsurf[:,i] = semic.tsurf
                alb[:,i]   = semic.alb
                netswd[:,i]= (1-alb[:,i])*force['swd'][:,i]

    # return ds file values
    semic_time = pandas.date_range(datetime.datetime(1980,1,1), periods=ntime)
    dssemic = xarray.Dataset(data_vars={'smb':(['nx','time'], smb),
                                       'melt':(['nx','time'], melt),
                                       'swsn':(['nx','time'], netswd),
                                       'tsurf':(['nx','time'], tsurf)},
                            coords={'x':('nx',md.mesh.x),
                                   'y':('nx',md.mesh.y),
                                   'time':('time',semic_time)})
    dssemic = dssemic.resample(time="1MS").mean(dim='time')

    del smb, melt, tsurf, alb, netswd

    # print('Evaluate matrix with cRMSE')
    matric1 = np.zeros((nx,))
    matric2 = np.zeros((nx,))
    matric3 = np.zeros((nx,))
    matric4 = np.zeros((nx,))
    
    for i in range(nx):
        if args.debug:
            matric1[i] = nCRMSE(dsracmo['tsurf'][i,:12].values, dssemic['tsurf'][:12,i].values, omitnan=1)
            matric2[i] = nCRMSE(dsracmo['melt'][i,:12].values,  dssemic['melt'][:12,i].values, omitnan=1)
            matric3[i] = nCRMSE(dsracmo['swsn'][i,:12].values,  dssemic['swsn'][:12,i].values, omitnan=1)
            matric4[i] = nCRMSE(dsracmo['smb'][i,:12].values,  dssemic['smb'][:12,i].values, omitnan=1)
        else:
            matric1[i] = nCRMSE(dsracmo['tsurf'][i,:].values, dssemic['tsurf'][:,i].values, omitnan=1)
            matric2[i] = nCRMSE(dsracmo['melt'][i,:].values,  dssemic['melt'][:,i].values, omitnan=1)
            matric3[i] = nCRMSE(dsracmo['swsn'][i,:].values,  dssemic['swsn'][:,i].values, omitnan=1)
            matric4[i] = nCRMSE(dsracmo['smb'][i,:].values,  dssemic['smb'][:,i].values, omitnan=1)

    matric_fin = np.sqrt(matric1**2 + matric2**2 + matric3**2 + matric4**2)
    areas   = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
    J = np.sum(areas*np.mean(matric_fin[md.mesh.elements-1],axis=1))

    return J

def evaluateSEMIC_wrap(part, md, mask_imbie, forc, queue, sema):
    if args.debug:
        nloop = 2
        ntime = 365
    else:
        nloop = 3
        ntime = 7671 # maximum available value
    # print(f'evalulate: ntime = {ntime}')
    # print(f'evalulate: nloop = {nloop}')
    J = evaluateSEMIC(part, md, mask_imbie, forc, ntime=ntime, nloop=nloop)
    queue.put(J)
    sema.release()


# In[ ]:


# Save each generation
def saveGeneration(fname, pop, best):
    data = {}
    for idx, part in enumerate(pop):
        data[idx] = {'part':part,
                  'fitness':part.fitness.values,
                  'speed':part.speed,
                  'local_best':part.best,
                  'local_best_fitness':part.best.fitness.values,
                  'global_best':best,
                  'global_best_fitness':best.fitness.values}
    with open(fname,'w') as fid:
        json.dump(data, fid)


# In[ ]:


toolbox = base.Toolbox()
toolbox.register("particle", generate, 
                 opts_bnd=opts)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=2.0, phi2=2.0)
toolbox.register("evaluate", evaluateSEMIC_wrap)

if isnotebook():
    pop = toolbox.population(n=10) # test
    print(pop[0])
    print(pop[0].smax)
    print(pop[0].smin)
    print(pop[0].speed)


# In[ ]:


random.seed(100)

npop = args.npop
GEN  = args.ngen
datadir = args.datadir #'./data/PSO-dummy/'
os.makedirs(datadir,exist_ok=1)

best = None
output_queue = multiprocessing.Queue()
semaphore = multiprocessing.Semaphore(args.ncpu)
pop = toolbox.population(n=npop)

# initialize log book
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean)
stats.register("std", numpy.std)
stats.register("min", numpy.min)
stats.register("max", numpy.max)

logbook = tools.Logbook()
logbook.header = ["gen", "evals"] + stats.fields

# initialize pandas dataframe
dflog = pandas.DataFrame(columns=['gen','evals','avg','std','min','max'])

tstart_glob = datetime.datetime.now()
for g in range(GEN):
    print('   Generation: %d/%d'%(g+1, GEN))
    procs  = []
    tstart = datetime.datetime.now()
    for part in tqdm.tqdm(pop, bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}', disable=(not args.tqdm)):
        semaphore.acquire()
        p = multiprocessing.Process(target=toolbox.evaluate,
                                    args=(part, md, mask_imbie, force, output_queue, semaphore),
                                   )
        procs.append(p)
        p.start()

    # print('   Join results....')
    for p in procs:
        p.join()

    # print('   get results')
    for part in pop:
        value = output_queue.get()
        part.fitness.values = [value]

    # print('   search global and local best.')
    for idx, part in enumerate(pop):
        if not part.best or part.best.fitness < part.fitness:
            part.best = creator.Particle(part)
            part.best.fitness.values = part.fitness.values
        if not best or best.fitness < part.fitness:
            print(f'-- Find global best in {idx}')
            best = creator.Particle(part)
            best.fitness.values = part.fitness.values

    # print('   update particle location.')
    # print('-- Previous')
    # print(pop[0])
    for part in pop:
        toolbox.update(part, best)
    # print('-- Updated')
    # print(pop[0])
    
    # Gather all the fitnesses in one list and print the stats
    logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
    print(logbook.stream)

    # okay, save each particle dataset
    saveGeneration(os.path.join(datadir,'PSO_gen%03d.json'%(g)), pop, best)

    dflog.loc[len(dflog)] = logbook[-1]
    dflog.to_csv(os.path.join(datadir,'PSO_summary.csv'))

    print(f'   Elapsed time: {datetime.datetime.now()-tstart}')

print(f'   Total elapsed time: {datetime.datetime.now()-tstart_glob}')

