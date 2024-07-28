#!/usr/bin/env python3

import pyseb
import numpy as np
import pytest

def test_SemicForcings():
    '''Test "SemicForcings" in "pyseb" module
    '''
    f = pyseb.SemicForcings()

    # Make random matrix
    nrow = 10
    ncol = 2
    A0 = np.random.random((nrow, ncol))

    # Assign A value in t2m
    f.t2m.set_value(A0)

    # Get value!
    print(f"t2m size = ({f.t2m.nrow}, {f.t2m.ncol})")
    print(f.t2m.get_value())

    assert(np.all(f.t2m.get_value() == A0))

if __name__ == '__main__':
    test_SemicForcings()
