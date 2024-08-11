#!/usr/bin/env python3
import platform

def iswindow():
    '''check wheater window or linux system. If window system, this function returns True value.

    Example
    -------
    code-block:: python
        
        imort pyseb
        iswindow = pyseb.utils.iswindow
    
    Returns
    -------
    flag: bool
        Return value true if window system.
    '''

    return platform.system() == 'Windows'

if __name__ == '__main__':
    print(iswindow())