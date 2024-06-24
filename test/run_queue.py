#!/usr/bin/env python3
import numpy as np
import os, sys, socket

# get hostname
hostname = socket.gethostname().lower().replace('-','')

# Okay, generate qscript
ngen  = 50
#npop  = 100
#omega = 0.8

def script_default(prefix, ngen, npop, omega): # {{{
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
        -opt_method=2 -datadir "data/{prefix}"
'''.format(ngen=ngen, npop=npop, omega=omega, prefix=prefix)

    return script
    # }}}

for npop in [100, 200]:
    for omega in [0.6, 0.8]:
        prefix = f'PSO_mon2_omega{omega}_npop{npop}' # initialize prefix of script

        # get script
        script = script_default(prefix, ngen, npop, omega)

        # okay, now write script
        with open(f'./run_{prefix}.sh','w') as fid:
            fid.write(script)

        # launch queue script.
        os.system(f'qsub ./run_{prefix}.sh')

