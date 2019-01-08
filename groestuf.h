#ifndef _PY4TI2_GROEBNERSTUF_

#define _PY4TI2_GROEBNERSTUF_

#include <Python.h>

#ifndef _4ti2_HAVE_GMP
PyObject *_4ti2ParticularSolution( PyObject *self, PyObject *args );
#endif

PyObject *_4ti2Minimize( PyObject *self, PyObject *args );
PyObject *_4ti2GroebnerBasis( PyObject *self, PyObject *args );
PyObject *_4ti2NormalForm( PyObject *self, PyObject *args );
PyObject *_4ti2MarkovBasis( PyObject *self, PyObject *args );
PyObject *_4ti2ZBasis( PyObject *self, PyObject *args );
PyObject *_4ti2Walk( PyObject *self, PyObject *args );

#endif

