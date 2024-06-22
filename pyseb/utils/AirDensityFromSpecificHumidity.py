#!/usr/bin/env python3
__all__ = ['AirDensityFromSpecificHumidity']
def AirDensityFromSpecificHumidity(Ps, T, q=None):
   '''
    Explain
     calculate air density with surface pressure, temperature, specific humidity.

     * for dry air condition.
       rho_air =  p_a / ( R_a * T )

       p_a : dry air pressure (unit: Pa)
       R_a : gas constant of dry air (unit: J/(kg K))
       T   : air temperature (unit: K)

     * for air vapor mixture
       rho_air = (p/R_a T) ( 1 + q ) / (1 + x R_w / R_a)
       R_(a,w) : gas constant of dry and vapor content (unit: J/(kg K))
       q       : specific humidity (kg/kg)
       T       : air temperature (unit: K)
       p       ; pressure in the moist air (Pa)

    Usage
     # calculate dry air density
     rho_air = AirDensitySpecificHumidity(Ps, T);

    Inputs
     Ps   - surface pressure (unit: Pa)
     T    - air temperature (unit: K)
     q    - specific humidity (unit: kg/kg)
             q = m_vapor / m_dry_air
             m_vapor   - water vapor mass
             m_dry_air - dry_air mass

    Outputs
     rho_air  - air density estimated with pressure and air temperature (unit: kg/m3)

    Reference
     https://www.engineeringtoolbox.com/density-air-d_680.html
   '''
   Ra = 287.058; # dry air constant J/kg/K
   Rw = 461.495; # water vapor constant J/kg/K
   R  = 8.31446; # J/K/mol
   if q is None: 
      # dry air density
      rho_air = Ps/Ra/T;
   else:
      # guess air density with specific humidity
      rho_air = Ps*(1+q)/T/(Ra+q*Rw);

   return rho_air

if __name__ == '__main__':
   import numpy as np
   import matplotlib.pyplot as plt
   Ps = 101325 # unit: Pa 
   T  = np.linspace(-50,0,50) + 273.15 # air tempeature (unit: K)
   # guess specific humidity in Antarctica
   q  = np.arange(5, 15,3)*1e-3 # specific humidity (unit: kg/kg)

   fig, ax = plt.subplots()
   rho_dry = AirDensityFromSpecificHumidity(Ps,T)
   ax.plot(T, rho_dry, label='dry air')
   for _q in q:
      rho_mix = AirDensitySpecificHumidity(Ps,T,_q)
      ax.plot(T,rho_mix, label='mixture {_q}')
   ax.legend()

   plt.show()

