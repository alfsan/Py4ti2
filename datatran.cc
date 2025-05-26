#include <Python.h>

// #include <iostream>
#include <string>

#include "groebner/Vector.h"
#include "groebner/VectorArray.h"
#include "groebner/DataType.h"

#include "4ti2/4ti2.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

std::string whathappened;

// Conversion routine PyUnicodeToString from 
// PyNormaliz module by Sebastian Gutsche.
std::string PyUnicodeToString( PyObject* in ){
#if PY_MAJOR_VERSION >= 3
    std::string out = "";
    int length = PyUnicode_GET_LENGTH( in );
    for( int i = 0; i < length; ++i )
        out += PyUnicode_READ_CHAR( in, i );
    return out;
#else
    char* out = PyString_AsString( in );
    return std::string(out);
#endif
}

// Conversion routine PyLongToIntegerType from PyLongToNmz in 
// PyNormaliz module by Sebastian Gutsche.

#ifdef _4ti2_GMP_
// IntegerType is of type mpz_class
bool PyLongToIntegerType( PyObject * in, IntegerType& out ){
  PyObject * in_as_string = PyObject_Str( in );
  const char* in_as_c_string = PyUnicodeToString( in_as_string ).c_str();
  out.set_str( in_as_c_string, 10 );
  return true;
}

// Conversion routine IntegerTypeToPyLong from NmzToPyLong in 
// PyNormaliz module by Sebastian Gutsche.

PyObject *IntegerTypeToPyLong( mpz_class input )
{
    std::string mpz_as_string = input.get_str();
    char *mpz_as_c_string = const_cast<char*>( mpz_as_string.c_str() );
    char *pend;
    PyObject *ret_val = PyLong_FromString( mpz_as_c_string, &pend, 10 );
    return ret_val;
}
#else

bool PyLongToIntegerType( PyObject *input, IntegerType& outi )
{
    long long val;
    int overflow;
    val = PyLong_AsLongLongAndOverflow( input, &overflow );
    if ( overflow == -1 ) {
        whathappened = "overflow detected converting to 4ti2 value";
        return false;
    }
#ifdef _4ti2_INT32_
    int nbits = 32;
#elif defined(_4ti2_INT64_)
    int nbits = 64;
#endif
    long long chk = val;
    int lbits = 0;
    while (chk > 0) {
        chk = chk >> 1;
        ++lbits;
    }
    if ( lbits > nbits ) {
#ifdef _4ti2_INT32_
        whathappened = "value size in bits exceeds 32 architecture";
#elif defined(_4ti2_INT64_)
        whathappened = "value size in bits exceeds 64 architecture";
#endif
        return false;
    }
    outi = (IntegerType) val;
    return true;
}

PyObject *IntegerTypeToPyLong( IntegerType input )
{
#ifdef _4ti2_INT32_
    return PyLong_FromLong( input );
#else
    return PyLong_FromLongLong( input );
#endif
}
#endif


