#!/usr/bin/env python3

from .libpysemic import SEMIC
import numpy as np

def prepare_semic_day(nx=12744, amp=3.34, alb_smax=0.859, alb_smin=0.4753, albi=0.485,hcrit=0.0037, rcrit=0.7031, mcrit=1.43e-8,): #{{{
    # print('Initialize structure')
    semic = SEMIC()
    # nx    = 12744

    semic.Initialize(nx)
    ONES = np.ones((nx,),dtype=float)
    semic.mask = 2*np.ones((nx,),dtype=int)
    semic.verbose = 0
    semic.alb_scheme = 3

    # This values are initial condition for calculating SEMIC.
    # Especially, "albedo" (alb) value is required for 'ISBA' albedo scheme.
    semic.alb = 0.8*ONES[:]
    semic.alb_snow = 0.8*ONES[:]
    semic.tsurf = 260.*ONES[:]
    semic.hsnow = 5*ONES[:]
    semic.hice  = 0*ONES[:]

    # update parameter
    semic.Param.ceff = 2e+6
    semic.Param.csh  = 2e-3
    semic.Param.clh  = 5e-4
    semic.Param.amp  = amp*ONES[:]
    semic.Param.alb_smax = alb_smax
    semic.Param.alb_smin = alb_smin
    semic.Param.albi = albi
    semic.Param.albl = 0.15
    semic.Param.hcrit = hcrit
    semic.Param.rcrit = rcrit
    semic.Param.mcrit = mcrit

    return semic
