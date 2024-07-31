#!/usr/bin/env python3
import pyseb
import numpy as np

def test_SemicResult():
    nx = 10
    ntime = 365
    result = pyseb.SemicResult()

    output_request = result.output_request
    output_list    = result.output_list

    print(f"output list = {output_list}")

    print("Done!")

if __name__ == '__main__':
    test_SemicResult()
