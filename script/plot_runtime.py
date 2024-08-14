#!/usr/bin/env python3
import numpy as np
import pandas
import os, sys, platform
import matplotlib.pyplot as plt

# load runtime with openmp
df = pandas.read_csv('./data/simba00_openmp_ERA5_runtime.csv', index_col=0)

# string type to seconds.
df.time = pandas.to_timedelta(df.time)
df.time = df.time.dt.total_seconds()
#df.time = pandas.to_timedelta(df.time, unit='s')

#print(df.time)

fig, ax = plt.subplots()
df.plot(x='num_thread',y='time',ax=ax)
ax.grid(True)
ax.set_xlabel('numter of thread')
ax.set_ylabel('time (seconds)')
fig.savefig('./sibma00_openmp_runtime.png',dpi=300)
plt.show()
