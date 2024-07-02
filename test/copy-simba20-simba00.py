#!/usr/bin/env python3
import sys, os
import pandas

expnames = []
expnames += ['PSO_day_omega0.6_npop100_Albdenby_ampMax10.00']
expnames += ['PSO_day_omega0.6_npop200_Albdenby_ampMax10.00']

expnames += ['PSO_mon2_omega0.6_npop100_ampMax10.00']
expnames += ['PSO_mon2_omega0.6_npop200_ampMax10.00']
expnames += ['PSO_mon2_omega0.6_npop100_Albdenby_ampMax10.00']
expnames += ['PSO_mon2_omega0.6_npop200_Albdenby_ampMax10.00']

expnames += ['PSO_mon2_omega0.6_npop100_Albisba_ampMax10.00']
expnames += ['PSO_mon2_omega0.6_npop200_Albisba_ampMax10.00']

if 1:
    SRC_DIR = os.path.abspath(f'./data')
    DEST_DIR= f'simba00:SEMIC-cpp/test/data/'
    command = f'rclone copy --progress --transfers=10 {SRC_DIR} {DEST_DIR}'
    for expname in expnames:
        print(expname)
        command += f' --include={expname}/*.json --include={expname}/*.csv'

    print(command)
    os.system(command)

# okay, now copy json file to .json file
for expname in expnames:
    dflog = pandas.read_csv(os.path.join('data',expname,'PSO_summary.csv'))
    maxGEN = dflog.gen.values[-1]
    print(maxGEN)

    SRC_FILE = os.path.join('data',expname,'PSO_gen%03d.json'%(maxGEN))
    DEST_FILE = os.path.join('data',expname + '.json')
    os.system(f'cp {SRC_FILE} {DEST_FILE}')
