#!/usr/bin/env python

from setuptools import setup
from setuptools.command.build_ext import build_ext
from distutils.core import Extension
import os

FTI2_DIR = os.environ.get('FTI2_DIR', '/usr')

module32 = Extension( "Py4ti2int32_cc",
                      define_macros = [('_4ti2_INT32_', '1')],
                      sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                          "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                          "qsolstuf.cc" ],
#                      include_dirs = [FTI2_DIR + "/include"],
#                      library_dirs = [FTI2_DIR + "/lib"], 
                      include_dirs = [os.path.join(FTI2_DIR, 'include')],
                      library_dirs = [os.path.join(FTI2_DIR, 'lib')], 
                      extra_link_args = ['-l4ti2int32', '-lglpk', '-l4ti2common', '-lzsolve' ] )

module64 = Extension( "Py4ti2int64_cc",
                       define_macros = [('_4ti2_INT64_', '1')],
                      sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                          "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                          "qsolstuf.cc" ],
                      include_dirs = [os.path.join(FTI2_DIR, 'include')],
                      library_dirs = [os.path.join(FTI2_DIR, 'lib')], 
                      extra_link_args=['-l4ti2int64', '-lglpk', '-l4ti2common', '-lzsolve' ] )

modulegmp = Extension( "Py4ti2gmp_cc",
                        define_macros = [('_4ti2_GMP_', '1')],
                        sources = [ "Py4ti2.cc", "datatran.cc", "4ti2mcnv.cc", 
                                    "zsolstuf.cc", "vecarcnv.cc", "groestuf.cc", 
                                    "qsolstuf.cc" ],
                        include_dirs = [os.path.join(FTI2_DIR, 'include')],
                        library_dirs = [os.path.join(FTI2_DIR, 'lib')], 
                        extra_link_args=['-l4ti2gmp', '-lglpk', '-lgmp', '-l4ti2common', '-lzsolve' ] )

modules_cc=[ module32, module64, modulegmp ]
modules_py=[ "Py4ti2int32", "Py4ti2int64", "Py4ti2gmp" ]

print(os.getenv('CC'))

CONFIG_HEADER = os.path.join(FTI2_DIR, 'include', '4ti2', '4ti2_config.h')
def check_4ti2_prec(p) :
    _4ti2_prec = "no"
    try:
        with open(CONFIG_HEADER, 'r', encoding='utf-8') as f:
            content = f.read()
            count = content.count(p)

        if count == 1:
            _4ti2_prec = "yes"
    except FileNotFoundError:
        print(f"Error: config file '{CONFIG_HEADER}' not found.")
    except Exception as e:
        print(f"Error: can not read config file.")
    return _4ti2_prec

if check_4ti2_prec("int32_t") == "no":
    modules_py.remove("Py4ti2int32")
    modules_cc.remove(module32)
    print("Py4ti2int32 will not be built")

if check_4ti2_prec("int64_t") == "no":
    modules_py.remove("Py4ti2int64")
    modules_cc.remove(module64)
    print("Py4ti2int64 will not be built")

if check_4ti2_prec("HAVE_GMP") == "no":
    modules_py.remove("Py4ti2gmp")
    modules_cc.remove(modulegmp)
    print("Py4ti2gmp will not be built")

setup(
    py_modules = modules_py, 
    ext_modules = modules_cc,
    package_data = {'': [ "COPYING", "GPLv2", "README.md" ] }
)
