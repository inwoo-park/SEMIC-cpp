#!/usr/bin/env python3

import os, sys

# change directory for window system.
cwd = os.path.dirname(__file__)
os.chdir(cwd)

# Do build script!
command = 'cd ../ && python setup.py build_ext --inplace'
os.system(command)
