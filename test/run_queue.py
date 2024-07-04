#!/usr/bin/env python3
import numpy as np
import os, sys, socket

# get hostname
hostname = socket.gethostname().lower().replace('-','')

# Okay, generate qscript
pso_method = 'APSO' # or 'wPSO'
ngen  = 100 
#freq = 'day'
freq = 'mon'
if freq == 'day':
    amp_max = 5.
elif freq == 'mon':
    amp_max = 10.
#alb_scheme = 'denby'
alb_scheme = 'isba'

def script_default(prefix, ngen, npop, freq, omega, amp_max, alb_scheme='denby', pso_method='wPSO'): # {{{
    script  = ''
    script += '#!/bin/bash\n'
    script += '''\
#PBS -N {prefix}
#PBS -l nodes=1:ppn=30
#PBS -o {prefix}.outlog
#PBS -e {prefix}.errlog
'''.format(prefix=prefix)

    script += '''
source ~/.bashrc

source ~/.bashrc
if [ $(hostname) == 'simba00' ]; then
    module load python/3.9
    module load simba/intel/mpich/mpich-3.2
    module load simba/intel/mpich/petsc-3.18
fi
module load gcc/cdo-1.9.3

cd $HOME/SEMIC-cpp/test/

# okay, run script
python3 run_OptimizationERA5.py -ngen {ngen} -npop {npop} -debug 0 -omega {omega} -ncpu=30 \
        -tqdm=0 \
        -freq={freq} -alb_scheme {alb_scheme} -amp_max={amp_max}\
        -opt_method=2 -datadir "data/{prefix}"
'''.format(ngen=ngen, npop=npop, omega=omega, amp_max=amp_max, freq=freq, prefix=prefix, alb_scheme=alb_scheme)

    return script
    # }}}

# check consistency {{{
if not freq in ['day','mon']:
    raise Exception('ERROR! with freq')
# }}}

# prepare the list of optimization.
for npop in [100, 200]: # 100, 200
    for omega in [0.7, 0.8]: # 0.6, 0.8
        if freq == 'mon':
            prefix = f'PSO_{freq}2_omega{omega}_npop{npop}_Alb{alb_scheme}_ampMax{amp_max:.2f}' # initialize prefix of script
        elif freq == 'day':
            prefix = f'PSO_{freq}_omega{omega}_npop{npop}_Alb{alb_scheme}_ampMax{amp_max:.2f}' # initialize prefix of script

        # get script
        script = script_default(prefix, ngen, npop, freq, omega, amp_max)

        # okay, now write script
        with open(f'./run_{prefix}.sh','w') as fid:
            fid.write(script)

        # launch queue script.
        os.system(f'qsub ./run_{prefix}.sh')
