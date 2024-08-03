#!/usr/bin/env python3
import pyseb
import numpy as np

def test_SemicResult(): # {{{
    nx = 10
    ntime = 20
    result = pyseb.SemicResult(nx, ntime)

    output_request = result.output_request
    output_list    = result.output_list

    evap = result.evap.get_value()
    subl = result.subl.get_value()

    hsnow = result.hsnow.get_value()
    alb_snow = result.alb_snow.get_value()
    alb = result.alb.get_value()

    assert(evap.shape == (nx, ntime))
    assert(subl.shape == (nx, ntime))

    assert(hsnow.shape == (nx, ntime))
    assert(alb_snow.shape == (nx, ntime))
    assert(alb.shape == (nx, ntime))
    # }}}

if __name__ == '__main__':
    test_SemicResult()
