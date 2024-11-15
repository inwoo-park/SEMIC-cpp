#!/usr/bin/env python3
import pybind11
import sys
import sysconfig # check default compiler for python.
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages, Extension
from distutils.sysconfig import get_config_vars
 
'''

Reference
* https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used
'''

# Initialize version.
__version__ = '0.4.1'

copt = {'linux':[]}
lopt = {'linux':[]}

# Check python's default compiler set.
def CheckCompiler(): # {{{
    compiler = sysconfig.get_config_var('CC').split(' ')[0]
    if compiler ==  'icc':
        return 'intel'
    elif compiler == 'gcc':
        return 'gcc'
    else:
        raise Exception(f"ERROR: Given compiler ({compiler}) is not supported.")
    # }}}

# OpenMP flags depending on the compiler
extra_compile_args = []
#extra_compile_args = ['-std=c++11']
if sys.platform == "win32":
    extra_compile_args += ['/openmp']
else:
    compiler = CheckCompiler()
    if compiler == 'gcc':
        extra_compile_args += ['-fopenmp']
    elif compiler == 'intel':
        extra_compile_args += ['-qopenmp'] #,'-xHost','-O2']

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the `get_include()`
    method can be invoked."""

    def __str__(self):
        return pybind11.get_include()

ext_modules = [
        Pybind11Extension(
            'pyseb.libpysemic', # install libpysemic in "pyseb" directory.
            ['src/libpysemic.cpp','src/SurfaceEnergyBalance.cpp'],
            define_macros=[('VERSION_INFO',__version__), ('HAVE_PYBIND11',1)],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_compile_args,
        ),
        Pybind11Extension(
            'pyseb.libexample', # install libpysemic in "pyseb" directory.
            ['src/libexample.cpp'],
            define_macros=[('VERSION_INFO',__version__), ('HAVE_PYBIND11', 1)],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_compile_args,
        ),
        ]

setup(
        name='pyseb',
        version=__version__,
        author='In-Woo Park',
        author_email='inwoopark0415@gmail.com',
        description='SEMIC in python and cpp version',
        ext_modules=ext_modules,
        extras_require={'test':'pytest'},
        cmdclass={'build_ext':build_ext},
        zip_safe=False,
        python_requires='>=3.9',
        #setup_requires=['pybind11'],
        install_requires=[
            'pybind11',
            'numpy<2',
            'matplotlib',
            'tqdm',
            'xarray',
            'netCDF4',
            'memory-profiler', # for checking memory leakage.
            'scipy','pytest',
            'nbmake', # for testing ipynb file.
        ],
        # which required directory?
        packages=find_packages(include=['pyseb']),
        # package_dir={'./':["*.so","*.pyd"]},
    )
