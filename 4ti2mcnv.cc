#include "Py4ti2.h"
#include "datatran.h"

#include <Python.h>

#include <string>

#include "groebner/DataType.h"

#include "4ti2/4ti2.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

bool PyIntListListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm )
{
    if ( !PyList_Check(input) ) {
        whathappened = "a list is expected for conversion to 4ti2 object";     
        return false;
    }

    const int nrows = PyList_Size(input);
    if ( nrows <= 0 ) {
        whathappened = "an non empty list was expected";
        return false;
    }

    PyObject *lstasrow = PyList_GetItem(input, 0);
    if ( !PyList_Check(lstasrow) ) {
        whathappened = "a list of lists was expected";
        return false;
    }
    const int length = PyList_Size(lstasrow);

    if ( _4ti2_state_create_matrix(state, nrows, 
                length, name, outm) != _4ti2_OK ) {
        whathappened = "this is serious, unable to create a 4ti2 matrix object";
        return false;
    }
#ifdef _4ti2_INT32_
//    _4ti2_int32_t value;
    int32_t value;
#elif defined(_4ti2_INT64_)
//    _4ti2_int64_t value;
    int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
    mpz_ptr ptrvalue = value.get_mpz_t();
#endif

    for (int i = 0; i < nrows; ++i ) {
        lstasrow = PyList_GetItem(input, i);
        if ( length != PyList_Size(lstasrow) ) {
            whathappened = "length of sublist differs while converting a list of lists to a 4ti2 object";
            return false;
        }
        for (int j = 0; j < length; ++j) {
            PyObject *lstval = PyList_GetItem(lstasrow, j);
            if (!PyLongToIntegerType(lstval, value) ) {
                whathappened += ", while list of lists to 4ti2 object conversion";
                return false;
            }
#ifdef _4ti2_INT32_
            _4ti2_matrix_set_entry_int32_t(*outm, i, j, value);
#elif defined(_4ti2_INT64_)
            _4ti2_matrix_set_entry_int64_t(*outm, i, j, value);
#elif defined(_4ti2_HAVE_GMP)
            _4ti2_matrix_set_entry_mpz_ptr(*outm, i, j, ptrvalue);
#endif
        }
    }
    return true;
}

bool PyIntListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm ) {
    if ( !PyList_Check(input) ) {
        whathappened = "a list is expected for conversion to 4ti2 matrix object"; 
        return false;
    }

    const int length = PyList_Size(input);
    if ( length < 0 ) {
        whathappened = "an non empty list was expected";
        return false;
    }
    PyObject *lstasrow = PyList_GetItem(input, 0);
    if ( PyList_Check(lstasrow) ) {
        whathappened = "a list of lists is not valid as input for conversion to 4ti2 matrix object"; 
        return false;
    }

    if ( _4ti2_state_create_matrix(state, 1, 
                length, name, outm) != _4ti2_OK ) {
        whathappened = "this is serious, unable to create a 4ti2 matrix object";
        return false;
    }
#ifdef _4ti2_INT32_
//    _4ti2_int32_t value;
    int32_t value;
#elif defined(_4ti2_INT64_)
//    _4ti2_int64_t value;
    int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
    mpz_ptr ptrvalue = value.get_mpz_t();
#endif

    for (int j = 0; j < length; ++j) {
        PyObject *lstval = PyList_GetItem( input, j );
        if (!PyLongToIntegerType(lstval, value) ) {
            whathappened += ", while list to 4ti2 object conversion";
            return false;
        }
#ifdef _4ti2_INT32_
        _4ti2_matrix_set_entry_int32_t(*outm, 0, j, value);
#elif defined(_4ti2_INT64_)
        _4ti2_matrix_set_entry_int64_t(*outm, 0, j, value);
#elif defined(_4ti2_HAVE_GMP)
        _4ti2_matrix_set_entry_mpz_ptr(*outm, 0, j, ptrvalue);
#endif
    }
    return true;
}

PyObject *_4ti2matrixToPyIntListList( _4ti2_matrix *m )
{
    const int nrows = _4ti2_matrix_get_num_rows( m );
    
    PyObject *intlstlst = PyList_New( nrows );

    const int ncols = _4ti2_matrix_get_num_cols( m );

#ifdef _4ti2_INT32_
//    _4ti2_int32_t value;
    int32_t value;
#elif defined(_4ti2_INT64_)
//    _4ti2_int64_t value;
    int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
    mpz_ptr ptrvalue = value.get_mpz_t();
#endif
    for ( int i = 0; i < nrows; ++i ) {
        PyObject *intlst = PyList_New( ncols );
        for ( int j = 0; j < ncols; ++j ) {
#ifdef _4ti2_INT32_
            _4ti2_matrix_get_entry_int32_t(m, i, j, &value);
#elif defined(_4ti2_INT64_)
            _4ti2_matrix_get_entry_int64_t(m, i, j, &value);
#elif defined(_4ti2_HAVE_GMP)
            _4ti2_matrix_get_entry_mpz_ptr(m, i, j, ptrvalue);
#endif
            PyList_SetItem( intlst, j, IntegerTypeToPyLong( value ) );
        }
        PyList_SetItem( intlstlst, i, intlst );
    }
    return intlstlst;
}


