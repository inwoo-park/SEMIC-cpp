#!/bin/bash
#PBS -N RunSEMIC_Optimization_day
#PBS -l nodes=1:ppn=36
#PBS -o RunSEMIC_Optimization_day.outlog 
#PBS -e RunSEMIC_Optimization_day.errlog 

module load python/3.9
module load simba/intel/mpich/mpich-3.2
module load simba/intel/mpich/petsc-3.18

# go to sepcific directory
cd $HOME/SEMIC-cpp/test

jupyter-nbconvert --to script ./run_OptimizationERA5.ipynb
#python3 run_OptimizationERA5.py -ngen 50 -npop 100 -opt_method=1 -datadir "data/PSO_mon1/" -freq "mon" -debug 0

#python run_OptimizationERA5.py -ngen 50 -npop 100 -debug 0 \
#	-opt_method=0 -datadir "data/PSO_day/" -freq "day"

#python3 run_OptimizationERA5.py -ngen 4 -npop 20 -debug 1 \
#	-opt_method=2 -datadir "data/PSO_mon2_debug/" -freq "mon"
#python3 run_OptimizationERA5.py -ngen 50 -npop 100 -debug 0 \
#	-opt_method=2 -datadir "data/PSO_mon2/" -freq "mon"

python3 run_OptimizationERA5.py -ngen 50 -npop 100 -debug 0 \
	-opt_method=2 -datadir "data/PSO_mon2/" -freq "mon"
