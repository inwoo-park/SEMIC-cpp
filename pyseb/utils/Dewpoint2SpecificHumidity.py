#!/usr/bin/env python3
from VaporPressure import VaporPressure

__all__ = ['Dewpoint2SpecificHumidity']
def Dewpoint2SpecificHumidity(dewpoint, Ps):
    '''
    Explain
     calculate specific humidity based on dew point and temperature.
    
    Usage
     q = Dewpoint2SpecificHumidity(Td, Ps)
    
    Inputs
     Td         : dewpoint temperature (unit: K)
     Ps         : surface total pressure (unit: Pa)
    
    Outputs
     qs  - specific humidity (unit: kg kg-1)
    
    Reference
     Specific Humidity (Bolton 1980):
    * https://www.npl.co.uk/resources/q-a/dew-point-and-relative-humidity
    '''
    #epsilon = 0.622; # Mw/Md - pressure of water vapor
    epsilon = 0.62195691; # more float precision.

    # get water vapor pressure
    Pv = VaporPressure(dewpoint,'bolton');

    # relationship between total pressure and vapor pressure
    # Pv = w / (epsilon + w) * P

    w = epsilon / (Ps /Pv - 1);

    # relationship between mixing ratio, w, and specificy humidity, q
    # q = w / (1+w)

    qs = w/(1+w); 

    return qs

if __name__ == '__main__':
    d2m = 273.15-10
    ps  = 1e+5
    qs = Dewpoint2SpecificHumidity(d2m, ps)

    print(qs)
