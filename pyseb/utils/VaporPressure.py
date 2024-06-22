#!/usr/bin/env python3
import numpy as np

__all__ = ['VaporPressure']
def VaporPressure(Ts,*args):
    '''
    Explain
     Estimate vapor pressure depending on surface temperature.
    
     pv = 611.2 * exp(a * (Ts-T0) / (Tb+Ts-T0) )
    
     a   : molar wiehgts of water vapour and dry air
     Tb  : back groundwater temperature (unit: K)
     T0  : freezing point of water (unit: K)
    
    Usage
     pv = VaporPressure(Ts);
     pv = VaporPressure(Ts,'water');
     pv = VaporPressure(Ts,'ice');
     
     # use Bolton 1980 relationship.
     pv = VaporPressure(Ts,'bolton');
    
    Inputs
     Ts:    surface temperature (unit: K)
    
    Outputs
     pv:    vapor pressure (Pa)
    
    References
    * Gill, A. E.: Atmosphere-Ocean Dynamics, International Geophysics Series, Academic Press, New York, Vol. 30, 1982
    * Krapp, M., Robinson, A., and Ganopolski, A.: SEMIC: an efficient surface energy and mass balance model applied to the Greenland ice sheet, The Cryosphere, 11, 1519–1535, https://doi.org/10.5194/tc-11-1519-2017, 2017.
    * Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    * Wallace, J. M. and Hobbs, P. V.: Atmospheric science: an introductory survey, Elsevier, 2006.
    '''

    ismethod = 'water';
    if len(args) == 2:
        ismethod = args[0]
    elif len(args) > 2:
        raise Exception('ERROR: current number of input argument is not available.')

    T0 = 273.15;
    if ismethod == 'water':
        a = 17.62;
        Tb = 243.12; # ~ 30.03
    elif ismethod == 'ice':
        a = 22.46;
        Tb = 272.62;
    elif ismethod == 'bolton': # bolton 1980
        a = 17.67;
        Tb = 273.15 - 29.65;
    else:
        raise Exception('ERROR: current method(=#s) is not available. available "ice", "water"'%(ismethod))

    # return output
    output = 611.2*np.exp(a*(Ts-T0)/(Ts+Tb-T0));

    return output
