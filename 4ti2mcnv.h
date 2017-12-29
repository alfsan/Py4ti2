#ifndef _PY4TI2_4TI2MCNV_
#define _PY4TI2_4TI2MCNV_

#include <string>

#include <Python.h>

#include "4ti2/4ti2.h"

extern std::string whathappened;

bool PyIntListListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm );

bool PyIntListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm );

PyObject *_4ti2matrixToPyIntListList( _4ti2_matrix *m );

#endif
