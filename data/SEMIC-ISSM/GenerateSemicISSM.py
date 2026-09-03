#!/usr/bin/env python3
# Set modules
import os, sys
import numpy as np
import xarray
import matplotlib.pyplot as plt

# Set environments
CURRENT_DIR=os.path.dirname(__file__)
os.chdir(CURRENT_DIR)

def InitializeSEMIC_Transient(commit:str, albedo:int=2):
    """Initialize module depend on commit

    commit = 893a942 (SEMIC-fortran)
    commit = 4da5301 (SEMIC-C++)

    Example
    -------
    .. code-block:: python
        >> md = InitializeSEMIC_Transient(commit)

    Parameters
    ----------
    albedo: int 
        0: None, 1: Slater, 2: Denby
    """
    isverbose = 0
    
    md = model();
    yts = md.constants.yts;
    md.miscellaneous.name = f'testSEMIC_{commit}';
    # contruct mesh
    xy = np.array([[0,1,0,1,0,1,0],
            [0,0,1,1,2,2,3],
            ]);
    elements = np.array([[2,3,1], [2,4,3], [3,4,5], [4,6,5], [5,6,7]])
    segments = np.array([[1,2,1], [1,3,1], [2,4,2],
                    [3,5,3], [4,6,4], [5,7,5], [6,7,5],
                    ])
    md.mesh.x = np.transpose(xy[0,:]);
    md.mesh.y = np.transpose(xy[1,:]);
    md.mesh.elements = elements;
    md.mesh.segments = segments;
    md.mesh.numberofvertices = len(md.mesh.x);
    md.mesh.numberofelements = len(md.mesh.elements);
    numberofvertices = md.mesh.numberofvertices;

    #print(md.mesh)

    # set initial geometry.
    md.geometry.surface = 10*np.ones((md.mesh.numberofvertices,1));
    md.geometry.base     = 0*np.ones((md.mesh.numberofvertices,1));
    md.geometry.bed      = 0*np.ones((md.mesh.numberofvertices,1));
    md.geometry.thickness = md.geometry.surface - md.geometry.base;

    # set mask.
    md.mask.ice_levelset    = -1*np.ones((md.mesh.numberofvertices,1))
    md.mask.ocean_levelset = np.ones((md.mesh.numberofvertices,1))

    md.smb = SMBsemic();

    md.timestepping.time_step=1/365;
    md.timestepping.start_time=0;
    md.timestepping.final_time=md.timestepping.time_step*365;
    md.settings.output_frequency=1;

    # number of grid cell...
    nx = 7
    inputs = np.loadtxt(os.path.join(CURRENT_DIR,'../SEMIC-Krapp2017/transect_input.txt'))
    times = np.reshape((np.arange(1,365+1)/365),(1,365))
    md.smb.dailysnowfall     = np.transpose(inputs[:,0*nx:1*nx]);
    md.smb.dailyrainfall     = np.transpose(inputs[:,1*nx:2*nx]);
    md.smb.dailydsradiation = np.transpose(inputs[:,2*nx:3*nx]);
    md.smb.dailydlradiation = np.transpose(inputs[:,3*nx:4*nx]);
    md.smb.dailywindspeed    = np.transpose(inputs[:,4*nx:5*nx]);
    md.smb.dailypressure     = np.transpose(inputs[:,5*nx:6*nx]);
    md.smb.dailyairdensity  = np.transpose(inputs[:,6*nx:7*nx]);
    md.smb.dailyairhumidity = np.transpose(inputs[:,7*nx:8*nx]);
    md.smb.dailytemperature = np.transpose(inputs[:,8*nx:9*nx]);
    fields = ['dailysnowfall','dailyrainfall','dailydsradiation',
              'dailydlradiation','dailywindspeed','dailypressure',
                 'dailyairdensity','dailyairhumidity','dailytemperature'];

    print('    n fields = %d'%(len(fields)))
    # add time.
    for i in range(len(fields)):
        if isverbose:
            print(getattr(md.smb,fields[i]).shape)
        setattr(md.smb,fields[i], np.concatenate((getattr(md.smb,fields[i]), times),axis=0))

    md.smb = md.smb.initialize(md);
    md.smb.alb_smax = 0.81;
    md.smb.alb_smin = 0.60;
    md.smb.albi  = 0.41;
    md.smb.albl  = 0.15;
    md.smb.tmax  = 273.15;
    md.smb.tmin  = 273.15-10;
    md.smb.hcrit = 0.028;
    md.smb.rcrit = 0.79;
    md.smb.mcrit = 6e-8;
    md.smb.Tamp  = 3.5;

    md.smb.mask = np.array([1,1,1,2,2,2,2]).T

    md.smb.albedo      = 0.8 * np.ones((numberofvertices,1))
    md.smb.albedo_snow = 0.8 * np.ones((numberofvertices,1))
    md.smb.hsnow       = 5*np.ones((md.mesh.numberofvertices,1))
    md.smb.hice        = 0*np.ones((md.mesh.numberofvertices,1))

    md.smb.ismethod = 1;
    # md.smb.albedo_scheme=1; # slater
    # md.smb.albedo_scheme=2; # denby
    #md.smb.albedo_scheme=3; # isba
    md.smb.albedo_scheme=albedo
    if md.smb.albedo_scheme == 1:
         alb_scheme_name = 'Slater'
    elif md.smb.albedo_scheme == 2:
         alb_scheme_name = 'Denby'
    elif md.smb.albedo_scheme == 3:
         alb_scheme_name = 'Isba'
    else:
         raise Exception("ERROR: Given albedo scheme is not available.")
    md.smb.qmr = np.zeros((md.mesh.numberofvertices,1))
    md.smb.requested_outputs = md.smb.outputlists(md);
    md.smb.requested_outputs += ['SmbSemicQmr','SmbEvaporation']

    md.initialization.temperature = 260*np.ones((md.mesh.numberofvertices,1))

    # set...
    md.verbose.smb = 0;
    md.verbose.solution = 0;

    md.transient = md.transient.deactivateall();
    md.transient.issmb = 1;
    md.cluster.np = 1;

    return md

def RunSemic(md,commit:str,maxloop:int=1):
    """
    commit = 893a942 (SEMIC-fortran)
    commit = 4da5301 (SEMIC-C++)
    """
    
    for nloop in range(maxloop):
        print(f'Loop %i/%i'%(nloop+1, maxloop))
        md = solve(md,'tr','runtimename',1);

        # update initial values
        md.smb.albedo        = md.results.TransientSolution[-1].SmbAlbedo;
        md.smb.albedo_snow = md.results.TransientSolution[-1].SmbAlbedoSnow;
        md.smb.hsnow         = md.results.TransientSolution[-1].SmbHSnow;
        md.smb.hice          = md.results.TransientSolution[-1].SmbHIce;
        md.initialization.temperature = md.results.TransientSolution[-1].TemperatureSEMIC;
        md.smb.qmr            = md.results.TransientSolution[-1].SmbSemicQmr;

    # load results
    time    = np.array([step.time for step in md.results.TransientSolution.steps])
    alb     = np.concatenate([step.SmbAlbedo for step in md.results.TransientSolution.steps],axis=1)
    melt    = np.concatenate([step.SmbMelt for step in md.results.TransientSolution.steps],axis=1)
    smb     = np.concatenate([step.SmbMassBalance for step in md.results.TransientSolution.steps],axis=1)
    smb_snow= np.concatenate([step.SmbMassBalanceSnow for step in md.results.TransientSolution.steps],axis=1)
    subl    = np.concatenate([step.SmbEvaporation for step in md.results.TransientSolution.steps],axis=1)

    tsurf = np.concatenate([step.TemperatureSEMIC for step in md.results.TransientSolution.steps],axis=1)

    hsnow = np.concatenate([step.SmbHSnow for step in md.results.TransientSolution.steps],axis=1)
    hice  = np.concatenate([step.SmbHIce for step in md.results.TransientSolution.steps],axis=1)
    qmr   = np.concatenate([step.SmbSemicQmr for step in md.results.TransientSolution.steps],axis=1)

    # Only save result with main script. Not under `pytest`.
    print('Export Result of transect result in netcdf file format')
    nc = xarray.Dataset(coords={'nx':np.arange(7),'time':np.arange(365)})
    nc['smb']   = xarray.DataArray(smb.T, dims=('time','nx'))
    nc['melt']  = xarray.DataArray(melt.T, dims=('time','nx'))
    nc['hsnow'] = xarray.DataArray(hsnow.T, dims=('time','nx'))
    nc['tsurf'] = xarray.DataArray(tsurf.T, dims=('time','nx'))
    nc['alb']   = xarray.DataArray(alb.T, dims=('time','nx'))
    nc['subl']  = xarray.DataArray(subl.T, dims=('time','nx'))

    return nc

if __name__ == '__main__':
    """
    Availalbe commit

    commit = 893a942 (SEMIC-fortran)
    commit = 4da5301 (SEMIC-C++)
    """
    MAXLOOP=2

    # Run SEMIC-fortran
    commit = '893a942'
    ISSM_DIR=os.path.join(CURRENT_DIR,'./ISSM-893a942')
    os.environ['ISSM_DIR']=ISSM_DIR
    sys.path.append(os.path.join(ISSM_DIR,'src/m/dev'))
    from devpath import *
    from model import model
    from SMBsemic import SMBsemic
    from solve import solve

    md1 = InitializeSEMIC_Transient(commit)
    md1.cluster.codepath=os.path.join(ISSM_DIR,'bin')
    # print(md1.cluster)
    nc1 = RunSemic(md1,commit,maxloop=MAXLOOP)
    nc1.to_netcdf(f'test_semic_transect_AlbDenby_{commit}.nc')

    # Run SEMIC-c++
    commit = '4da5301'
    ISSM_DIR2=os.path.join(CURRENT_DIR,'./ISSM-4da5301')
    os.environ['ISSM_DIR']=ISSM_DIR2
    sys.path.append(os.path.join(ISSM_DIR2,'src/m/dev'))
    from devpath import *
    from model import model
    from SMBsemic import SMBsemic
    from solve import solve
    
    md2 = InitializeSEMIC_Transient(commit)
    md2.cluster.codepath=os.path.join(ISSM_DIR2,'bin')
    # print(md2.cluster)
    nc2 = RunSemic(md2,commit,maxloop=MAXLOOP)
    nc2.to_netcdf(f'test_semic_transect_AlbDenby_{commit}.nc')
    