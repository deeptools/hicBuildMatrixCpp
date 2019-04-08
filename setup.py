from setuptools.command.build_ext import build_ext
from distutils.core import setup, Extension
import platform
import sys
from sysconfig import get_config_var, get_paths
import logging

__version__ = '0.0.1'

def get_include(): 
    info = get_paths()
    seqan_path = '/'.join(info['include'].split('/')[:-1])
    # seqan_path += '/seqan'
    return seqan_path

def __extra_compile_args():
    extra_compile_args = []

    if platform.system() == 'Darwin':
        extra_compile_args = ["-std=c++14"]
    else:
        extra_compile_args = ["-fopenmp", "-std=c++14", "-O3"]
    return extra_compile_args


def __extra_link_args():
    extra_link_args = []
    if platform.system() != 'Darwin':
        extra_link_args = ["-lgomp", "-lm", "-lrt", "-lpthread", "-lz",  "-lbz2"]
    return extra_link_args
# -O3 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt

sources_list = ['src/hicbuildmatrix_interface.cpp', 'src/hicbuildmatrix.cpp']
depends_list = ['src/hicbuildmatrix.hpp']

module_hicbuildmatrix = Extension('hicBuildMatrixCpp',
                      sources=sources_list,
                      depends = depends_list,
                      include_dirs=[
                          # Path to eigen3 headers
                          get_include()
                      ],
                      extra_link_args=__extra_link_args(),
                      extra_compile_args=__extra_compile_args()
                      )


setup(
    name='hicBuildMatrixCpp',
    version=__version__,
    author='Klesta Toti, Joachim Wolff',
    author_email='wolffj@informatik.uni-freiburg.de',
    description='A c++ extension for python to create a Hi-C interaction matrix',
    ext_modules=[module_hicbuildmatrix],
    install_requires=['pybind11>=2.2'],
   headers = ['src/krbalancing.hpp']
)