#include "datatran.h"
#include "4ti2mcnv.h"

#include <Python.h>

#include <iostream>
#include <string>
// #include <vector>

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

PyObject *_4ti2Graver( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }

    _4ti2_state *graver_api;
#ifdef _4ti2_INT32_
    graver_api = _4ti2_graver_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    graver_api = _4ti2_graver_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    graver_api = _4ti2_graver_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[5] = {"mat", "lat", "sign", "ub", "lb"};
    _4ti2matrix_InputData graver_input( 5, input_types, graver_api );
    
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( graver_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !graver_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( graver_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    char *argv[2] = { (char*)"graver", (char*)"-q" };

    if ( _4ti2_state_set_options( graver_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( graver_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(graver_api) != _4ti2_OK ) {
        _4ti2_state_delete( graver_api );
        PyErr_SetString(Py4ti2Error, "graver computation error");
        return NULL;
    }
    
    _4ti2_matrix* zhom_matrix;
    _4ti2_state_get_matrix(graver_api, "zhom", &zhom_matrix);
    PyObject *zhom_list = NULL;
    if ( zhom_matrix != 0 ) {
       zhom_list = _4ti2matrixToPyIntListList( zhom_matrix );
    }

    _4ti2_state_delete( graver_api );

    return zhom_list;
}

PyObject *_4ti2Hilbert( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }
    
    _4ti2_state *hilbert_api;
#ifdef _4ti2_INT32_
    hilbert_api = _4ti2_hilbert_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    hilbert_api = _4ti2_hilbert_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    hilbert_api = _4ti2_hilbert_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[5] = {"mat", "lat", "sign", "rel", "ub"}; 
    _4ti2matrix_InputData hilbert_input( 5, input_types, hilbert_api );

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( hilbert_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !hilbert_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( hilbert_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    char *argv[2] = { (char*)"hilbert", (char*)"-q" };

    if ( _4ti2_state_set_options( hilbert_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( hilbert_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(hilbert_api) != _4ti2_OK ) {
        _4ti2_state_delete( hilbert_api );
        PyErr_SetString(Py4ti2Error, "hilbert computation error");
        return NULL;
    }

    int elresul = 0;
    PyObject *result = PyTuple_New(4);

    _4ti2_matrix* zhom_matrix;
    _4ti2_state_get_matrix(hilbert_api, "zhom", &zhom_matrix);
    PyObject *zhom_list;
    if ( zhom_matrix != 0 ) {
        zhom_list = _4ti2matrixToPyIntListList( zhom_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zhom") );
        PyTuple_SetItem( result, elresul++, zhom_list );
    }

    _4ti2_matrix* zfree_matrix;
    _4ti2_state_get_matrix(hilbert_api, "zfree", &zfree_matrix);
    PyObject *zfree_list;
    if ( zfree_matrix != 0 ) {
        zfree_list = _4ti2matrixToPyIntListList( zfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zfree") );
        PyTuple_SetItem( result, elresul, zfree_list );
    }

    _4ti2_state_delete( hilbert_api );

    return result;
}

PyObject *_4ti2Zsolve( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }
    
    _4ti2_state *zsolve_api;
#ifdef _4ti2_INT32_
    zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_64);
#elif defined(_4ti2_GMP_)
    zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_ARB);
#endif

    const char *input_types[7] = {"mat", "lat", "rhs", "sign", "rel", "ub", "lb"};
    _4ti2matrix_InputData zsolve_input( 7, input_types, zsolve_api );

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( zsolve_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);

        if ( !zsolve_input.read( typeofinp, eol2 ) ) {
            _4ti2_state_delete( zsolve_api );
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }
    
    char *argv[2] = { (char*)"zsolve", (char*)"-q" };

    if ( _4ti2_state_set_options( zsolve_api, 2, argv) != _4ti2_OK ) {
        _4ti2_state_delete( zsolve_api );
        PyErr_SetString(Py4ti2Error, "Unexpected error");
        return NULL;
    }

    if ( _4ti2_state_compute(zsolve_api) != _4ti2_OK ) {
        _4ti2_state_delete( zsolve_api );
        PyErr_SetString(Py4ti2Error, "zsolve computation error");
        return NULL;
    }

    PyObject *result = PyTuple_New(6);
    int elresul = 0;

    _4ti2_matrix* zinhom_matrix;
    _4ti2_state_get_matrix(zsolve_api, "zinhom", &zinhom_matrix);
    PyObject *zinhom_list;
    if ( zinhom_matrix != 0 ) {
        zinhom_list = _4ti2matrixToPyIntListList( zinhom_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zinhom") );
        PyTuple_SetItem( result, elresul++, zinhom_list );
    } 

    _4ti2_matrix* zhom_matrix;
    _4ti2_state_get_matrix(zsolve_api, "zhom", &zhom_matrix);
    PyObject *zhom_list;
    if ( zhom_matrix != 0 ) {
        zhom_list = _4ti2matrixToPyIntListList( zhom_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zhom") );
        PyTuple_SetItem( result, elresul++, zhom_list );
    }

    _4ti2_matrix* zfree_matrix;
    _4ti2_state_get_matrix(zsolve_api, "zfree", &zfree_matrix);
    PyObject *zfree_list;
    if ( zfree_matrix != 0 ) {
        zfree_list = _4ti2matrixToPyIntListList( zfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zfree") );
        PyTuple_SetItem( result, elresul, zfree_list );
    }

    _4ti2_state_delete( zsolve_api );

    return result;
}

