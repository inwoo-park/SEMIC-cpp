#!/usr/bin/env python3

import pyseb
import numpy as np
import os
import pytest

# for checking memory usage
from memory_profiler import profile
import psutil

@profile
def test_SemicForcings():
    '''Test "SemicForcings" in "pyseb" module
    '''

    # Make random matrix
    nrow = 365
    ncol = 1000
    A0 = np.random.random((nrow, ncol))
    f = pyseb.SemicForcings()

    # Assign A value in t2m
    f.sf.set_value(A0)
    f.rf.set_value(A0)
    f.t2m.set_value(A0)

    # Get value!
    print(f"t2m size = ({f.t2m.nrow}, {f.t2m.ncol})")

    k = f.t2m.get_value_python()
    k = f.sf.get_value_python()
    assert(np.all(np.array(k) == A0))

    # clear memory
    del A0
    del k
    del f

if __name__ == '__main__':
    proc = psutil.Process(os.getpid())
    mem_before = proc.memory_info().rss / 1024**2
    test_SemicForcings()
    mem_after  = proc.memory_info().rss / 1024**2
    print('Memory usage')
    print(f'   before = {mem_before}')
    print(f'   after  = {mem_after}')
