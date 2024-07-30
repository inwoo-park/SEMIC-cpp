#!/usr/bin/env python3
import pyseb
import numpy as np

def test_SemicResult():
    nx = 10
    ntime = 365
    result = pyseb.SemicResult()

    print(result.output_request)

if __name__ == '__main__':
    test_SemicResult()