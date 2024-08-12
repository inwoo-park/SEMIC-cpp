#!/usr/bin/env python3

def get_cmip6_modelname():
    '''return available CMIP6 models.

    Example
    -------
    code-block:: python
        import pyseb
        MODELNAMES = pyseb.io.get_cmip6_modelname()
    
    Returns
    -------
    output: list
        Available CMIP6 models.
    '''
    return [
    'ACCESS-CM2',
    'ACCESS-ESM1-5',
    'AWI-CM-1-1-MR',
    'BCC-CSM2-MR',
    # 'CAS-ESM2-0',
    'CMCC-CM2-SR5',
    'CMCC-ESM2',
    'CanESM5',
    'EC-Earth3',
    'EC-Earth3-Veg',
    'EC-Earth3-Veg-LR',
    'FGOALS-f3-L',
    'GFDL-ESM4',
    'INM-CM4-8',
    'INM-CM5-0',
    'IPSL-CM6A-LR',
    'KIOST-ESM',
    'MIROC6',
    'MPI-ESM1-2-HR',
    'MPI-ESM1-2-LR',
    'MRI-ESM2-0',
    ]