#!/usr/bin/env python3
import pybind11
import sys
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages, Extension
# from setuptools.command.build_ext import build_ext
 
'''

Reference
* https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used
'''

# Initialize version.
__version__ = '0.1'

copt = {'linux':[]}
lopt = {'linux':[]}

# OpenMP flags depending on the compiler
if sys.platform == "win32":
    omp_flags = ['/openmp']
else:
    omp_flags = ['-fopenmp']

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
            define_macros=[('VERSION_INFO',__version__)],
            extra_compile_args=omp_flags,
            extra_link_args=omp_flags,
        ),
        Pybind11Extension(
            'pyseb.libexample', # install libpysemic in "pyseb" directory.
            ['src/libexample.cpp'],
            define_macros=[('VERSION_INFO',__version__)],
            extra_compile_args=omp_flags,
            extra_link_args=omp_flags,
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
        setup_requires=['pybind11'],
        install_requires=[
            'pybind11',
            'numpy',
            'matplotlib',
            'tqdm',
            'xarray',
            'netCDF4',
        ],
        # which required directory?
        packages=find_packages(include=['pyseb']),
        # package_dir={'./':["*.so","*.pyd"]},
    )
