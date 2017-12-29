#ifndef _PY4TI2_DATATRAN_

#define _PY4TI2_DATATRAN_

#include <string>

#include <Python.h>

#include "groebner/DataType.h"
#include "4ti2/4ti2.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

extern std::string whathappened;

std::string PyUnicodeToString( PyObject* in );

bool PyLongToIntegerType( PyObject * in, IntegerType& out );

#ifdef _4ti2_GMP_
PyObject *IntegerTypeToPyLong( mpz_class input );
#else
PyObject *IntegerTypeToPyLong( IntegerType input );
#endif

#endif
