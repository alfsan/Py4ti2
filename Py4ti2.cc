#include "zsolstuf.h"
#include "groestuf.h"
#include "qsolstuf.h"

#include <Python.h>

PyObject *Py4ti2Error;

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyMethodDef Py4ti2Methods[] = {
//#ifndef _4ti2_HAVE_GMP
//    {"solve", (PyCFunction)_4ti2ParticularSolution, METH_VARARGS, 
//        "Computes a particular solution of a linear diophantine equations system" },
//#endif
    {"_minimize", (PyCFunction)_4ti2Minimize, METH_VARARGS, 
        "Computes de minimal solutions of an integer linear program or, more general, a lattice program, using Groebner basis." },
    {"_groebner", (PyCFunction)_4ti2GroebnerBasis, METH_VARARGS, 
        "Computes a Groebner basis of the toric ideal of a matrix, or, more general, of the lattice ideal of a lattice." },
    {"_normalform", (PyCFunction)_4ti2NormalForm, METH_VARARGS, 
        "Computes the normal form of a list of feasible points." },
    {"_markov", (PyCFunction)_4ti2MarkovBasis, METH_VARARGS, 
        "Computes a Markov basis (generating set) of the toric ideal of a matrix or, more general, of the lattice ideal of a lattice." },
    {"_zbasis", (PyCFunction)_4ti2ZBasis, METH_VARARGS, 
        "Computes an integer lattice basis." },
    {"_walk", (PyCFunction)_4ti2Walk, METH_VARARGS, 
        "Computes the minimal solution of an integer linear program or, more general, a lattice program using a Groebner basis." },
    {"_hilbert", (PyCFunction)_4ti2Hilbert, METH_VARARGS,
        "Computes the Hilbert basis of a matrix or a given lattice."},
    {"_graver", (PyCFunction)_4ti2Graver, METH_VARARGS,
        "Computes the Graver basis of a matrix or a given lattice."},
    {"_zsolve", (PyCFunction)_4ti2Zsolve, METH_VARARGS,
        "Solves linear inequality and equation systems over the integers."},
    {"_circuits", (PyCFunction)_4ti2Circuits, METH_VARARGS,
        "Computes the circuits of a cone."},
    {"_qsolve", (PyCFunction)_4ti2Qsolve, METH_VARARGS,
        "Computes a generator description of a cone."},
    {"_rays", (PyCFunction)_4ti2Rays, METH_VARARGS,
        "Computes the extreme rays of a cone."},
    {NULL, NULL, 0, NULL }        /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static int Py4ti2_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int Py4ti2_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "Py4ti2_cc",
        NULL,
        sizeof(struct module_state),
        Py4ti2Methods,
        NULL,
        Py4ti2_traverse,
        Py4ti2_clear,
        NULL
};

#define INITERROR return NULL

#else

#define INITERROR return

#endif


#if PY_MAJOR_VERSION >= 3

PyMODINIT_FUNC PyInit_Py4ti2int32_cc(void)

#else

extern "C" void initPy4ti2int32_cc(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("Py4ti2int32", Py4ti2Methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(const_cast<char*>("Py4ti2.INITError"), NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    
    Py4ti2Error = PyErr_NewException(const_cast<char*>("4ti2.interface_error"), NULL, NULL );
    Py_INCREF( Py4ti2Error );
    
    PyModule_AddObject( module, "error", Py4ti2Error );
    
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

#if PY_MAJOR_VERSION >= 3

PyMODINIT_FUNC PyInit_Py4ti2int64_cc(void)

#else

extern "C" void initPy4ti2int64_cc(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("Py4ti2int64", Py4ti2Methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(const_cast<char*>("Py4ti2.INITError"), NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    
    Py4ti2Error = PyErr_NewException(const_cast<char*>("4ti2.interface_error"), NULL, NULL );
    Py_INCREF( Py4ti2Error );
    
    PyModule_AddObject( module, "error", Py4ti2Error );
    
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

#if PY_MAJOR_VERSION >= 3

PyMODINIT_FUNC PyInit_Py4ti2gmp_cc(void)

#else

extern "C" void initPy4ti2gmp_cc(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("Py4ti2gmp", Py4ti2Methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(const_cast<char*>("Py4ti2.INITError"), NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    
    Py4ti2Error = PyErr_NewException(const_cast<char*>("4ti2.interface_error"), NULL, NULL );
    Py_INCREF( Py4ti2Error );
    
    PyModule_AddObject( module, "error", Py4ti2Error );
    
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
