#ifndef _PY4TI2_VECARCNV_

#define _PY4TI2_VECARCNV_

#include <string>

#include <Python.h>

#include "groebner/Vector.h"
#include "groebner/VectorArray.h"

extern std::string whathappened;

PyObject *VectorToPyIntList( _4ti2_::Vector& v );

PyObject *VectorArrayToPyIntListList( _4ti2_::VectorArray& va );

bool PyIntListToVector( PyObject *input, _4ti2_::Vector& outv );

bool PyIntListListToVectorArray( PyObject *input, _4ti2_::VectorArray& outva );
#endif
