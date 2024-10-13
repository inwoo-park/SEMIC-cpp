#!/usr/bin/env python3

import warnings

__all__ = ['get_cmip6_modelname']
def get_cmip6_modelname(exclude:list=['FGOALS-f3-L']):
    '''return available CMIP6 models. CAS-ESM2-0 and FGOALS-f3-L provides weird values in monthly forcing value!

    Example
    -------
    code-block:: python
        import pyseb
        MODELNAMES = pyseb.io.get_cmip6_modelname()

    Parameters
    ----------
    exclude: list (default: None)
        Exclude specific model defined in return value
    
    Returns
    -------
    output: list
        Available CMIP6 models convering SEMIC with both SSP1-2.6 and SSP5-8.5 scenarios.
    '''

    modelname = [
    'ACCESS-CM2',
    'ACCESS-ESM1-5',
    'AWI-CM-1-1-MR',
    'BCC-CSM2-MR',
    #'CAS-ESM2-0',
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

    if exclude != None:
        if isinstance(exclude,str):
            exclude=[exclude]
        # First, check excluded model list
        for name in exclude:
            if not name in modelname:
                warnings.warn('WARNING: Given exclude model name (=%s) is not contained in pre-defined CMIP6 model list.'%(name))

        new_name = []
        for name in modelname:
            if not name in exclude:
                new_name.append(name)
            else:
                print(f'   exclude: {name}')
        modelname = new_name # update

    return modelname

if __name__ == '__main__':
    names = get_cmip6_modelname(exclude=['FGOALS-f3-L'])
    print(names)

    # trigger error!
    names = get_cmip6_modelname(exclude=['FGOALS-f3'])
