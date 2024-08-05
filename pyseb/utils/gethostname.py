#!/usr/bin/env python3
import socket

def gethostname():
    ''' Get hostname of current machine

    Example
    -------
    code-block::python

        import pyseb
        hostname = pyseb.utils.gethostname() 

    Returns
    -------
    hostname: str
        current machine's name.
    '''
    return socket.gethostname().lower().replace('-','')

if __name__ == '__main__':
    gethostname()
