#ifndef _PY4TI2_4TI2MCNV_
#define _PY4TI2_4TI2MCNV_

#include <vector>
#include <string>

#include <Python.h>

#include "4ti2/4ti2.h"

extern std::string whathappened;

bool PyIntListListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm );

bool PyIntListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
        const char* name, _4ti2_matrix **outm );

PyObject *_4ti2matrixToPyIntListList( _4ti2_matrix *m );

extern PyObject * Py4ti2Error;

struct _4ti2matrix_InputData {
    std::vector<std::string> input_type_str;
    _4ti2_matrix **data; // *mat, *lat, *rhs, *rel, *sign, *ub, *lb;
    _4ti2_state *state;

    // n is the length of x
    _4ti2matrix_InputData( int n, const char *x[], _4ti2_state *s ) {
        assert( n > 0);
        for (int i = 0; i < n ; ++i) input_type_str.push_back( std::string(x[i]) );
        state = s;
        data = new _4ti2_matrix *[n];
        for (int i = 0; i < n; ++i) data[i] = 0;
    }

    ~_4ti2matrix_InputData() {
        delete[] data;
    }


    int input_type2index( std::string& i_t ) {
        unsigned int i = 0;
        do { 
            if (i_t == input_type_str[i++] ) 
                return i-1;
        } while (i < input_type_str.size());
        return -1;
    }

    bool read( std::string& input_type, PyObject *pydata ) {
        int idx = input_type2index( input_type );
        if ( idx < 0 ) {
//            PyErr_SetString(Py4ti2Error, "Unexpected argument");
            whathappened = "unexpected argument "+input_type;
            return false;
        }
        PyObject *fiel = PyList_GetItem(pydata, 0);
        if ( !PyList_Check( fiel ) ) {
            if ( !PyIntListTo4ti2matrix( pydata, state, input_type.c_str(), 
                        &(data[idx]) ) ) {
                std::string pref = "\'" + input_type + "\' argument: ";
                whathappened = pref + whathappened;
                return false;
            }
        }
        else { 
            if ( !PyIntListListTo4ti2matrix( pydata, state, input_type.c_str(), 
                        &(data[idx]) ) ) {
                std::string pref = "\'" + input_type + "\' argument: ";
                whathappened = pref + whathappened;
                return false;
            }
        }
        return true;
    }
};

#endif
