#!/usr/bin/env python3
import pyseb
import numpy as np
import xarray, netCDF4
import os, sys
import scipy.io # load mat file
import argparse
import logging
import tqdm
import datetime
import pandas
import shutil
import json

# Go to specific directory
cwd = os.path.dirname(__file__)
os.chdir(cwd)
print(f'Current directory {cwd}')

logging.basicConfig(encoding='utf-8',
)
					# stream=sys.stdout)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
# logger = logger.StreamHandler(sys.stdout)

logger.info('Initialize argument parsing')
parser = argparse.ArgumentParser()
parser.add_argument('prefix',type=str,
                    help='Give which prefix do you want to run semic?')
parser.add_argument('-freq',default='mon',type=str,
                    help='Forcing variable frequency. "mon" and "day" are available. (default: mon)')
parser.add_argument('-alb_scheme',type=str,default='isba',
                    help='Select specific albedo scheme. (default: isba)')
parser.add_argument('-nloop',default=5,type=int,
                    help='number of literation for SEMIC. (default: 5)')
parser.add_argument('-debug',default=0,type=int,
                    help='Debugging current script. (default: False)')
args = parser.parse_args()

logger.info('Load Mat file')
if args.freq == 'day':
    if args.debug:
        forc = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980.mat')
    else:
        forc = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Day_1980,2000.mat')
elif args.freq == 'mon':
    if args.debug:
	    forc = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Month_1980.mat')
    else:
        forc = scipy.io.loadmat('../data/Prepare/ANT_InterpERA5_Month_1980,2000.mat')
else:
	logger.error("ERROR: Given freq vairlabe (={arg.freq}) is not available. Only day and mon are available.")

# check consistency
if not args.alb_scheme in ['slater','denby','isba']:
    raise Exception('ERROR: Given albedo scheme (=%s) is not available. slater, denby and isba are available.')

logger.info('-- calculate snow and rainfall flux')
rho_freshwater = 1000 # kg m-3
forc['sf'] = forc['msr']/rho_freshwater # unit: kg m-2 s-1 > water m/s
forc['rf'] = (forc['mtpr'] - forc['msr'])/rho_freshwater
forc['lwd'] = forc['msdwlwrf']
forc['swd'] = forc['msdwswrf']

logger.info('-- calculate specific humidity and air density')
forc['qq']   = pyseb.utils.Dewpoint2SpecificHumidity(forc['d2m'], forc['sp'])
forc['rhoa'] = pyseb.utils.AirDensityFromSpecificHumidity(forc['sp'], forc['t2m'], q=forc['qq'])

ntime = forc['t2m'].shape[0]
nx    = forc['t2m'].shape[1]
nloop = args.nloop

logger.info('Load ANT model domain')
nc = netCDF4.Dataset('../data/ANT_Mesh.nc')
# print(nc)
elements = np.array(nc['/mesh/elements'][:],dtype=int)
x = nc['/mesh/x'][:]
y = nc['/mesh/y'][:]
logger.info(f'-- mesh.x = {x}')

logger.info('Interpolate IMBIE2 basin.')
mask_imbie = pyseb.interp.interpXarrayGridToMesh('../data/IMBIE2/basinNumbers_8km.nc',x, y,xname='x',yname='y',varname='basinNumber',method='nearest')
mask_imbie = np.array(mask_imbie.values, dtype=int)

logger.info('Initialize SEMIC.')
semic = pyseb.prepare_semic_day(nx)

logger.info("Prepare SEMIC parameters depending on frequency")
ONES=np.ones((nx,)) # for grid information.
semic.alb_scheme = 3 # Use denby alb snow parameterziation
if args.freq == 'day':
    logger.info('-- Catch daily parameter!')
    '''
    Tamp:     3.7010
    hcrit:    0.0061
    rcrit:    0.8927
    mcrit:    6.22e-09
    alb_smax: 0.9125
    alb_smin: 0.4832
    albi:     0.4651
    '''
    # Okay, we load best PSO-SEMIC results forced with daily dataset.
    with open('./data/PSO_day_omega0.6_npop200.json','r') as fid:
        part = json.load(fid)['0']['global_best']
    semic.Param.amp   = part['amp']*ONES[:]
    semic.Param.hcrit = 10**part['hcrit']
    semic.Param.rcrit = part['rcrit']
    semic.Param.mcrit = 10**part['mcrit']
    semic.Param.alb_smax = part['alb_smax']
    semic.Param.alb_smin = part['alb_smin']
    semic.Param.albi 	= part['albi']
elif args.freq == 'mon':
    logger.info('-- load semic parameters for monthly forcing.')
    with open('./data/PSO_mon2_{args.prefix}.json','r') as fid:
          part = json.load(fid)['0']['global_best']
    amp = ONES.copy()
    for nb in range(16):
        pos = (mask_imbie == nb)
        amp[pos] = part['amp%d'%(nb)]
    semic.Param.amp = amp[:]

    # set albedo scheme
    if args.alb_scheme == 'slater':
        semic.alb_scheme = 1
    elif args.alb_scheme == 'denby':
        semic.alb_scheme = 2
    elif args.alb_scheme == 'isba':
        semic.alb_scheme = 3
    semic.Param.alb_smax = part['alb_smax']
    semic.Param.alb_smin = part['alb_smin']
    semic.Param.albi = part['albi']
    semic.Param.albl = 0.15
    semic.Param.hcrit = 10**part['hcrit']
    semic.Param.rcrit = part['rcrit']
    semic.Param.mcrit = 10**part['mcrit']
    
# initialize results array
smb     = np.zeros((nx,ntime))
melt    = np.zeros((nx,ntime))
tsurf   = np.zeros((nx,ntime))
alb     = np.zeros((nx,ntime))
netswd  = np.zeros((nx,ntime))
shf     = np.zeros((nx,ntime))
lhf     = np.zeros((nx,ntime))

for loop in range(nloop):
    print(f'   nloop: {loop}')
    for i in tqdm.tqdm(range(ntime)):
        semic.sf   = forc['sf'][i,:]
        semic.rf   = forc['rf'][i,:]
        semic.swd  = forc['swd'][i,:]
        semic.lwd  = forc['lwd'][i,:]
        semic.wind = forc['wind2m'][i,:]
        semic.sp   = forc['sp'][i,:]
        semic.rhoa = forc['rhoa'][i,:]
        semic.qq   = forc['qq'][i,:]
        semic.t2m  = forc['t2m'][i,:]
        semic.RunEnergyAndMassBalance()
        # break
        if loop == nloop-1:
            smb[:,i]   = semic.smb
            melt[:,i]  = semic.melt
            tsurf[:,i] = semic.tsurf
            alb[:,i]   = semic.alb
            netswd[:,i]= (1-alb[:,i])*forc['swd'][i,:]

logger.info('Export output in xarray.')
rho_freshwater = 1000
yts = 365*24*3600
# freshwater m3 sec-1 to Gt yr-1
TotalSmb  = pyseb.utils.TransientTotalTimeSeries(elements-1, x, y ,smb,) * rho_freshwater*yts/1e+12
TotalMelt = pyseb.utils.TransientTotalTimeSeries(elements-1, x, y, melt) * rho_freshwater*yts/1e+12

semic_time = pandas.date_range(datetime.datetime(1980,1,1), periods=ntime)
dssemic = xarray.Dataset(data_vars={'smb':(['nx','time'], smb),
                                    'melt':(['nx','time'], melt),
                                    'swsn':(['nx','time'], netswd),
                                    'tsurf':(['nx','time'], tsurf),
                                    'alb':(['nx','time'], alb),
                                    'TotalSmb':(['time'], TotalSmb),
                                    'TotalMelt':(['time'], TotalMelt),
									},
                            coords={'x':('nx',x),
                                    'y':('nx',y),
                                    'time':('time',semic_time)})
dssemic = dssemic.resample(time="1MS").mean(dim='time')

dssemic = dssemic.to_netcdf()
ofname = f'./ANT_SEMIC_ERA5_{args.freq}_{args.prefix}.nc'
logger.info('Export output: %s'%(ofname))
with open(ofname,'wb') as fid:
      fid.write(dssemic)
