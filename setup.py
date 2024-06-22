#!/usr/bin/env python3
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages
'''

Reference
* https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used
'''

# Initialize version.
__version__ = '0.1'

copt = {'linux':[]}
lopt = {'linux':[]}

ext_modules = [
        Pybind11Extension(
            '_pySEMIC',
            ['src/pySEMIC.cpp','src/SurfaceEnergyBalance.cpp'],
            define_macros=[('VERSION_INFO',__version__)],
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
        # which required directory?
        packages=find_packages(include=['pyseb']),
    )
