#!/usr/bin/env python3
import numpy as np
import os, sys, platform

def isnotebook():
   '''Check wheater notebook or not
   '''
   try:
      from IPython import get_ipython
      if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
         return False
   except ImportError:
      return False
   except AttributeError:
      return False
   return True

if __name__ == '__main__':
    print(isnotebook())
