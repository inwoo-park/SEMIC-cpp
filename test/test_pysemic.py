#!/usr/bin/env python3
import numpy as np
import os, sys, platform
sys.path.append('../build')
import pySEMIC

semic = pySEMIC.SEMIC()
semic.Initialize(10)
semic.Display()

