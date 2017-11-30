#!/usr/bin/env python

from distutils.core import setup, Extension
import sys

module32 = Extension( "Py4ti2int32",
                      define_macros = [('_4ti2_INT32_', '1')],
                      sources = [ "Py4ti2.cc" ],
                      extra_link_args=['-l4ti2int32', '-lglpk', '-l4ti2common', '-lzsolve' ] )

module64 = Extension( "Py4ti2int64",
                      define_macros = [('_4ti2_INT64_', '1')],
                      sources = [ "Py4ti2.cc" ],
                      extra_link_args=['-l4ti2int64', '-lglpk', '-l4ti2common', '-lzsolve' ] )

setup(
    name = 'Py4ti2',
    version = '0.1',
    description = 'An interface to some 4ti2 functions',
    author = 'Alfredo Sanchez-R.-Navarro',
    author_email = 'alfredo.sanchez@uca.es',
#    url = 'https://github.com/Normaliz/PyNormaliz',
    ext_modules = [ module32, module64 ],
#    package_data = {'': [ "COPYING", "GPLv2" ] },
)
