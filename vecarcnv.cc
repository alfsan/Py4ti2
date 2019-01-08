#include "Py4ti2.h"
#include "datatran.h"

#include <Python.h>

#include <string>

#include "groebner/Vector.h"
#include "groebner/VectorArray.h"
#include "groebner/DataType.h"

PyObject *VectorToPyIntList( _4ti2_::Vector& v )
{
    const int size = v.get_size();
    PyObject *intlst = PyList_New( size );
    for ( int i = 0; i< size; ++i )
        PyList_SetItem( intlst, i, IntegerTypeToPyLong( v[i] ) );

    return intlst;
}

PyObject *VectorArrayToPyIntListList( _4ti2_::VectorArray& va )
{
    const int numb = va.get_number();
    const int size = va.get_size();
    IntegerType value;
    PyObject *intlstlst = PyList_New( numb );
    for ( int i = 0; i< numb; ++i ) {
        PyObject *intlst = PyList_New( size );

        for ( int j = 0; j < size; ++j ) {
            value = va[i][j];
            PyList_SetItem( intlst, j, IntegerTypeToPyLong( value ) );
        }
        
        PyList_SetItem( intlstlst, i, intlst );
    }
    return intlstlst;
}

bool PyIntListToVector( PyObject *input, _4ti2_::Vector& outv )
{
    if ( !PyList_Check(input) ) {
        whathappened = "a list is expected for conversion to 4ti2 object";
        return false;
    }
    int length = PyList_Size(input);
    if ( length != outv.get_size() ) {
        whathappened = "conversion from list to 4ti2 object failed, \
                        dimension of input data mismatch";
        return false;
    }
    for (int i = 0; i < length; ++i) {
        PyObject* tmp = PyList_GetItem(input, i);
        if ( !PyLongToIntegerType(tmp, outv[i]) ) {
            whathappened += ", while integer list to 4ti2 object conversion";
            return false;
        }
    }
    return true;
}

bool PyIntListListToVectorArray( PyObject *input, _4ti2_::VectorArray& outva )
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
    PyObject *tmp = PyList_GetItem(input, 0);
    if( !PyList_Check(tmp) ) {
        whathappened = "a list of lists was expected";
        return false;
    }

    int length = PyList_Size(tmp);
    _4ti2_::Vector rowv(length);

    _4ti2_::VectorArray _va(0, length);

    if ( !PyIntListToVector(tmp, rowv) ) {
        whathappened += ", while list of lists 4ti2 object conversion";
        return false;
    }
    _va.insert( rowv );


    for (int i = 1; i < nrows; ++i ) {
        tmp = PyList_GetItem(input, i);
        if ( length != PyList_Size(tmp) ) {
            whathappened = "length of sublist differs while converting a list of lists to a 4ti2 object";
            return false;
        }
        if ( !PyIntListToVector(tmp, rowv) ) {
            whathappened += ", while list of lists 4ti2 object conversion";
            return false;
        }
        _va.insert( rowv );
    }

    outva = _va;
    
    return true;
}

