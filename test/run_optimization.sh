#!/bin/bash

jupyter-nbconvert --to script ./run_OptimizationERA5.ipynb
#python3 run_OptimizationERA5.py -ngen 50 -npop 100 -opt_method=1 -datadir "data/PSO_mon1/" -freq "mon" -debug 0

#python3 run_OptimizationERA5.py -ngen 4 -npop 20 -debug 1 \
#	-opt_method=2 -datadir "data/PSO_mon2_debug/" -freq "mon"
python3 run_OptimizationERA5.py -ngen 50 -npop 100 -debug 0 \
	-opt_method=2 -datadir "data/PSO_mon2/" -freq "mon"
