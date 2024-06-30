#!/usr/bin/env python3
import sys, os

experiment = []

experiment.append({'prefix':'omega0.6_npop100_Albdenby_ampMax10.00',
                   'alb_scheme':'denby',
                   'freq':'mon'})

experiment.append({'prefix':'omega0.6_npop100_ampMax10.00',
                   'alb_scheme':'isba',
                   'freq':'mon'})

for exp in experiment:
    print(exp)
    prefix = exp['prefix']
    freq   = exp['freq']
    alb_scheme = exp['alb_scheme']
    command = f'python run_semic.py {prefix} -freq {mon} -alb_scheme={alb_scheme}'
    os.system(command)