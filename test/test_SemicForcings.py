#!/usr/bin/env python3

import pyseb
import numpy as np
import pytest

from memory_profiler import profile

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
    #print(f.t2m.get_value())

    k = f.t2m.get_value()
    assert(np.all(np.array(k) == A0))

    # clear memory
    del A0
    del k
    del f

if __name__ == '__main__':
    test_SemicForcings()
