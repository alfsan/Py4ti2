#include "datatran.h"
#include "4ti2mcnv.h"

#include <Python.h>

#include <iostream>
#include <string>
#include <vector>

#include "4ti2/4ti2.h"

// #include "glpk.h"

// #ifdef _4ti2_HAVE_GMP
// #include <gmp.h>
// #endif

#if PY_MAJOR_VERSION >= 3
#define string_check PyUnicode_Check
#else
#define string_check PyString_Check
#endif

extern PyObject * Py4ti2Error;

PyObject *_4ti2Circuits( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }

    _4ti2_state *circuits_api;
#ifdef _4ti2_INT32_
    circuits_api = _4ti2_circuits_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    circuits_api = _4ti2_circuits_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    circuits_api = _4ti2_circuits_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[5] = {"mat", "sign", "rel"};
    _4ti2matrix_InputData circuits_input( 3, input_types, circuits_api );
    
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( circuits_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !circuits_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( circuits_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    char *argv[2] = { (char*)"circuits", (char*)"-q" };

    if ( _4ti2_state_set_options( circuits_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( circuits_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(circuits_api) != _4ti2_OK ) {
        _4ti2_state_delete( circuits_api );
        PyErr_SetString(Py4ti2Error, "circuits computation error");
        return NULL;
    }

    PyObject *result = PyTuple_New(4);
    int elresul = 0;

    _4ti2_matrix* cir_matrix;
    _4ti2_state_get_matrix(circuits_api, "cir", &cir_matrix);
    PyObject *cir_list = NULL;
    if ( cir_matrix != 0 ) {
        cir_list = _4ti2matrixToPyIntListList( cir_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("cir") );
        PyTuple_SetItem( result, elresul++, cir_list );
    }

    _4ti2_matrix* qfree_matrix;
    _4ti2_state_get_matrix(circuits_api, "qfree", &qfree_matrix);
    PyObject *qfree_list = NULL;
    if ( qfree_matrix != 0 ) {
        qfree_list = _4ti2matrixToPyIntListList( qfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("qfree") );
        PyTuple_SetItem( result, elresul++, qfree_list );
    }

    _4ti2_state_delete( circuits_api );

    return result;
}

PyObject *_4ti2Qsolve( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }

    _4ti2_state *qsolve_api;
#ifdef _4ti2_INT32_
    qsolve_api = _4ti2_qsolve_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    qsolve_api = _4ti2_qsolve_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    qsolve_api = _4ti2_qsolve_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[5] = {"mat", "sign", "rel"};
    _4ti2matrix_InputData qsolve_input( 3, input_types, qsolve_api );
    
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( qsolve_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !qsolve_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( qsolve_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    char *argv[2] = { (char*)"qsolve", (char*)"-q" };

    if ( _4ti2_state_set_options( qsolve_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( qsolve_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(qsolve_api) != _4ti2_OK ) {
        _4ti2_state_delete( qsolve_api );
        PyErr_SetString(Py4ti2Error, "qsolve computation error");
        return NULL;
    }
    
    PyObject *result = PyTuple_New(4);
    int elresul = 0;

    _4ti2_matrix* qhom_matrix;
    _4ti2_state_get_matrix(qsolve_api, "qhom", &qhom_matrix);
    PyObject *qhom_list = NULL;
    if ( qhom_matrix != 0 ) {
        qhom_list = _4ti2matrixToPyIntListList( qhom_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("qhom") );
        PyTuple_SetItem( result, elresul++, qhom_list );
    }

    _4ti2_matrix* qfree_matrix;
    _4ti2_state_get_matrix(qsolve_api, "qfree", &qfree_matrix);
    PyObject *qfree_list = NULL;
    if ( qfree_matrix != 0 ) {
        qfree_list = _4ti2matrixToPyIntListList( qfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("qfree") );
        PyTuple_SetItem( result, elresul++, qfree_list );
    }

    _4ti2_state_delete( qsolve_api );

    return result;
}

PyObject *_4ti2Rays( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }

    _4ti2_state *rays_api;
#ifdef _4ti2_INT32_
    rays_api = _4ti2_rays_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    rays_api = _4ti2_rays_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    rays_api = _4ti2_rays_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[5] = {"mat", "sign", "rel"};
    _4ti2matrix_InputData rays_input( 3, input_types, rays_api );
    
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( rays_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !rays_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( rays_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    char *argv[2] = { (char*)"rays", (char*)"-q" };

    if ( _4ti2_state_set_options( rays_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( rays_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(rays_api) != _4ti2_OK ) {
        _4ti2_state_delete( rays_api );
        PyErr_SetString(Py4ti2Error, "rays computation error");
        return NULL;
    }
    
    PyObject *result = PyTuple_New(4);
    int elresul = 0;

    _4ti2_matrix* ray_matrix;
    _4ti2_state_get_matrix(rays_api, "ray", &ray_matrix);
    PyObject *ray_list = NULL;
    if ( ray_matrix != 0 ) {
        ray_list = _4ti2matrixToPyIntListList( ray_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("ray") );
        PyTuple_SetItem( result, elresul++, ray_list );
    }

    _4ti2_matrix* qfree_matrix;
    _4ti2_state_get_matrix(rays_api, "qfree", &qfree_matrix);
    PyObject *qfree_list = NULL;
    if ( qfree_matrix != 0 ) {
        qfree_list = _4ti2matrixToPyIntListList( qfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("qfree") );
        PyTuple_SetItem( result, elresul++, qfree_list );
    }

    _4ti2_state_delete( rays_api );

    return result;
}

