#!/usr/bin/env python

from distutils.core import setup, Extension
import os

module32 = Extension( "Py4ti2int32_cc",
                      define_macros = [('_4ti2_INT32_', '1')],
                      sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                          "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                          "qsolstuf.cc" ],
                      extra_link_args=['-l4ti2int32', '-lglpk', '-l4ti2common', '-lzsolve' ] )

module64 = Extension( "Py4ti2int64_cc",
                       define_macros = [('_4ti2_INT64_', '1')],
                      sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                          "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                          "qsolstuf.cc" ],
                      extra_link_args=['-l4ti2int64', '-lglpk', '-l4ti2common', '-lzsolve' ] )

modulegmp = Extension( "Py4ti2gmp_cc",
                       define_macros = [('_4ti2_GMP_', '1')],
                      sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                          "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                          "qsolstuf.cc" ],
                      extra_link_args=['-l4ti2gmp', '-lglpk', '-lgmp', '-l4ti2common', '-lzsolve' ] )

modules_cc=[ module32, module64, modulegmp ]
modules_py=[ "Py4ti2int32", "Py4ti2int64", "Py4ti2gmp" ]
if os.getenv('_4ti2_INT32', None) == None:
    modules_py.remove("Py4ti2int32")
    modules_cc.remove(module32)
    print("Py4ti2int32 will not be built")

if os.getenv('_4ti2_INT64', None) == None:
    modules_py.remove("Py4ti2int64")
    modules_cc.remove(module64)
    print("Py4ti2int64 will not be built")

if os.getenv('_4ti2_HAVE_GMP', None)==None:
    modules_py.remove("Py4ti2gmp")
    modules_cc.remove(modulegmp)
    print("Py4ti2gmp will not be built")

if len(modules_cc)!=0:
    setup(
        name = 'Py4ti2',
        version = '0.5',
        description = 'An interface to some 4ti2 functions',
        author = 'Alfredo Sanchez-R.-Navarro',
        author_email = 'alfredo.sanchez@uca.es',
        url = 'https://github.com/alfsan/Py4ti2',
        py_modules = modules_py, 
        ext_modules = modules_cc,
        package_data = {'': [ "COPYING", "GPLv2", "README.md" ] }
    )
