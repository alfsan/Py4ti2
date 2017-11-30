#include <Python.h>

#include <iostream>
#include <sstream>
#include <streambuf>
#include <string>

#include "groebner/VectorStream.h"
#include "groebner/Vector.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/VectorArray.h"
#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/GeneratingSet.h"
#include "groebner/Globals.h"
#include "groebner/GroebnerBasis.h"
#include "groebner/DataType.h"
#include "groebner/Minimize.h"
#include "groebner/LatticeBasis.h"


#include "4ti2/4ti2.h"

#include "glpk.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

#if PY_MAJOR_VERSION >= 3
#define string_check PyUnicode_Check
#else
#define string_check PyString_Check
#endif

static PyObject * Py4ti2Error;

std::string whathappened;

IntegerType mod( IntegerType a, IntegerType b )
{
    return ( a % b + b) % b;
}

void diolinsys_instance_col( _4ti2_::VectorArray& A, _4ti2_::Vector& b, 
        int torsion, _4ti2_::Size j, IntegerType val, _4ti2_::Vector& bp, 
        bool& no_hom )
{
    assert( 0 <= j && j < A.get_size() );

    _4ti2_::Index colt = A.get_size() - torsion;
    no_hom = false;
    for( _4ti2_::Size k = 0; k < A.get_number(); ++k ) {
        if ( k < torsion ) 
            bp[k] = mod( b[k] - ( val * A[k][j] ), A[k][colt + k] );
        else bp[k] = b[k] - val * A[k][j];
        if ( bp[k] ) no_hom = true;
    }
}

void insert_1st_column( const _4ti2_::Vector& vs1, 
        const _4ti2_::VectorArray& vs2, _4ti2_::VectorArray& vs )
{
    assert(vs1.get_size() == vs2.get_number());
    assert(vs.get_size() == vs2.get_size()+1);
    for( _4ti2_::Index i = 0; i < vs2.get_number(); ++i )
        vs[i][0] = vs1[i];
    _4ti2_::VectorArray::lift( vs2, 1, vs2.get_size()+1, vs );
}

// Orden producto: lexicografico en la 1a variable
// y reverso lexicográfico graduado en el resto
int cost_definition_gmdp( _4ti2_::VectorArray& cost )
{
    _4ti2_::Size ncol, nfil;
    
    ncol = cost.get_size();
    nfil = cost.get_number();
    assert( ncol == nfil );
    
    _4ti2_::Size i, j;
    for (i=0; i<ncol; i++) cost[0][i] = 0;
    cost[1][0] = 0;
    cost[0][0] = 1;
    for (i=1; i<ncol; i++) cost[1][i] = 1;
    
    for (j=2; j<nfil; j++) {
        for (i=0; i<ncol; i++) cost[j][i] = 0;
        cost[j][ncol-j+1]=-1;
    }
    return 1;
}

// Orden reverso lexicografico graduado 
int cost_definition_dp( _4ti2_::VectorArray& cost )
{
    _4ti2_::Size ncol, nfil;
    
    ncol = cost.get_size();
    nfil = cost.get_number();
    assert( ncol == nfil );
    
    _4ti2_::Size i, j;
    for (i=0; i<ncol; i++) cost[0][i] = 1;
    
    for (j=1; j<nfil; j++) {
        for (i=0; i<ncol; i++) cost[j][i] = 0;
        cost[j][ncol-j]=-1;
    }
    return 1;
}

bool solve_by_ideal_computations( _4ti2_::VectorArray& A, int torsion, 
        _4ti2_::Vector& b, bool no_hom, _4ti2_::Vector& _particular )
{
    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::VectorArray *bA, *cost;
    _4ti2_::Feasible *feasible;

    if ( no_hom ) {
        bA = new _4ti2_::VectorArray( A.get_number(), A.get_size()+1 );
        insert_1st_column( b, A, *bA );
    }
    else bA = &A;
    _4ti2_::Size numbelem, sizeelem;
    _4ti2_::Index i, j;
    _4ti2_::VectorArray fbas;
    if ( torsion > 0 ) {
        _4ti2_::VectorArray lat( 0, bA->get_size() );
        _4ti2_::lattice_basis( *bA, lat );

        numbelem = lat.get_number();

        if ( numbelem == 0) {
            // Restore normal standard output
            std::cout.rdbuf( old );

            if ( no_hom ) delete bA;
            return false;
        }

        sizeelem = lat.get_size() - torsion;
        _4ti2_::VectorArray latwotorsion( numbelem, sizeelem );
        for( i = 0; i < numbelem; ++i )
            for( j = 0; j < sizeelem; ++j )
                latwotorsion[i][j] = lat[i][j];
    
        feasible = new _4ti2_::Feasible( &latwotorsion );
        fbas = feasible->get_basis(); 
    }
    else {
        sizeelem = bA->get_size();
        feasible = new _4ti2_::Feasible (0, bA);
        fbas = feasible->get_basis(); 
    }

    if ( fbas.get_number() == 0 ) {
        // Restore normal standard output
        std::cout.rdbuf( old );
        if ( no_hom ) delete bA;
        delete feasible;
        return false;
    }

    _4ti2_::Globals::minimal = false;
    _4ti2_::GeneratingSet gs( *feasible, 0 );

    cost = new _4ti2_::VectorArray( sizeelem, sizeelem );
    if ( no_hom ) cost_definition_gmdp( *cost );
    else cost_definition_dp( *cost );

    _4ti2_::GroebnerBasis gb( gs, cost );

    // Restore normal standard output
    std::cout.rdbuf( old );

    _4ti2_::VectorArray binomios = gb.get_groebner_basis();

    if (no_hom) delete bA;
    delete cost;
    delete feasible;

    IntegerType binx;

    numbelem = binomios.get_number();
    if ( numbelem == 0 ) return false;

    sizeelem = binomios[0].get_size();
    bool vale = false;
    if ( no_hom ) {
        // Hay solucion si solo uno de los generadores tiene todas las
        // coordenadas del mismo signo (salvo las que sean nulas)
        i = 0;
        while ( i < numbelem && !vale ) {
            binx = binomios[i][0];
            if ( binx == 1 || binx == -1 ) {
                for( j = 1; j < sizeelem && binx * binomios[i][j] <= 0; ++j ) {}
                vale = ( j == sizeelem );
            }
            ++i;
        }
        if ( vale )
            for( --i, j = 1; j < sizeelem; ++j )
                _particular[j-1] = -binx * binomios[i][j];
    }
    else {
        i = 0;
        while ( i < numbelem && !vale ) {
            for ( j = 0; j < sizeelem && binomios[i][j] == 0; ++j) {}
            if ( j < sizeelem ) {
                binx = binomios[i][j];
                for( ++j; j < sizeelem && binx * binomios[i][j] >= 0; ++j ) {}
                vale = ( j == sizeelem );
            }
            ++i;
        }
        if ( vale ) {
            binx = binx > 0 ? 1 : -1;
            for( --i, j = 0; j< sizeelem; ++j )
                _particular[j] = binx * binomios[i][j];
        }
    }
    return vale;
}

// Conversion routine PyLongToIntegerType from PyLongToNmz in 
// PyNormaliz module by Sebastian Gutsche.

#ifdef _4ti2_WITH_GMP_
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
    string mpz_as_string = input.get_str();
    char *mpz_as_c_string = const_cast<char*>( mpz_as_string.c_str() );
    char *pend;
    PyObject *ret_val = PyLong_FromString( mpz_as_c_string, &pend, 10 );
    return ret_val;
}
#else

static bool PyLongToIntegerType( PyObject *input, IntegerType& outi )
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

static bool PyIntListToVector( PyObject *input, _4ti2_::Vector& outv )
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

PyObject *VectorToPyIntList( _4ti2_::Vector& v )
{
    const int size = v.get_size();
    PyObject *intlst = PyList_New( size );
    for ( int i = 0; i< size; ++i )
        PyList_SetItem( intlst, i, IntegerTypeToPyLong( v[i] ) );

    return intlst;
}

PyObject *_4ti2matrixToPyIntListList( _4ti2_matrix *m )
{
    const int nrows = _4ti2_matrix_get_num_rows( m );
    PyObject *intlstlst = PyList_New( nrows );

    const int ncols = _4ti2_matrix_get_num_cols( m );

#ifdef _4ti2_INT32_
    _4ti2_int32_t value;
#elif defined(_4ti2_INT64_)
    _4ti2_int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
#endif
    for ( int i = 0; i < nrows; ++i ) {
        PyObject *intlst = PyList_New( ncols );
        for ( int j = 0; j < ncols; ++j ) {
#ifdef _4ti2_INT32_
            _4ti2_matrix_get_entry_int32_t(m, i, j, &value);
#elif defined(_4ti2_INT64_)
            _4ti2_matrix_get_entry_int64_t(m, i, j, &value);
#elif defined(_4ti2_HAVE_GMP)
            _4ti2_matrix_get_entry_mpz_ptr(m, i, j, &value);
#endif
            PyList_SetItem( intlst, j, IntegerTypeToPyLong( value ) );
        }
        PyList_SetItem( intlstlst, i, intlst );
    }
    return intlstlst;
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

static bool PyIntListListToVectorArray( PyObject *input, _4ti2_::VectorArray& outva )
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

static bool PyIntListListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
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
    _4ti2_int32_t value;
#elif defined(_4ti2_INT64_)
    _4ti2_int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
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
            _4ti2_matrix_set_entry_mpz_ptr(*outm, i, j, value);
#endif
        }
    }
    return true;
}

static bool PyIntListTo4ti2matrix( PyObject *input, _4ti2_state* state, 
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
    _4ti2_int32_t value;
#elif defined(_4ti2_INT64_)
    _4ti2_int64_t value;
#elif defined(_4ti2_HAVE_GMP)
    mpz_class value;
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
        _4ti2_matrix_set_entry_mpz_ptr(*outm, 0, j, value);
#endif
    }
    return true;
}

// Conversion routine PyUnicodeToString from 
// PyNormaliz module by Sebastian Gutsche.
std::string PyUnicodeToString( PyObject* in ){
#if PY_MAJOR_VERSION >= 3
    std::string out = "";
    int length = PyUnicode_GET_SIZE( in );
    for( int i = 0; i < length; ++i )
        out += PyUnicode_READ_CHAR( in, i );
    return out;
#else
    char* out = PyString_AsString( in );
    return std::string(out);
#endif
}

/*
 * Python module stuff
 */

static PyObject *_4ti2ParticularSolution( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    _4ti2_::VectorArray A;
    PyObject *eol = PyTuple_GetItem(args, 0);

    if ( !PyIntListListToVectorArray(eol, A) ) {
        std::string pref = "1st argument: ";
        whathappened = pref + whathappened;
        PyErr_SetString(Py4ti2Error, whathappened.c_str());
//                "1st argument must be a list of lists of integers");
        return NULL;
    }
    _4ti2_::Vector b( A.get_number(), 0 );
    bool no_hom = false;
    if ( nargs > 1 ) {
        eol = PyTuple_GetItem(args, 1);
        if ( !PyIntListToVector(eol, b) ) {
            std::string pref = "2nd argument: ";
            whathappened = pref + whathappened;
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
//                    "2nd argument must be a list of integers");
            return NULL;
        }
        no_hom = true;
    }
    _4ti2_::Vector sol( A.get_size() );
    if ( solve_by_ideal_computations( A, 0, b, no_hom, sol) ) {
        return VectorToPyIntList(sol);
    }
    else 
        return PyList_New(0);
}

struct _4ti2GroebnerBasisInput {
    _4ti2_::VectorArray *mat, *lat, *weights, *mar, *cost; 
    _4ti2_::Vector *sign, *weightsmax, *zsol;

    _4ti2GroebnerBasisInput() {
        mat = lat = weights = mar = cost = 0; 
        sign = weightsmax = zsol = 0;
    }

    void clearallbutmar() {
        if ( mat != 0 ) delete mat;
        if ( lat != 0 ) delete lat;
        if ( weights != 0 ) delete weights;
        if ( weightsmax != 0 ) delete weightsmax;
        if ( cost != 0 ) delete cost;
        if ( sign != 0 ) delete sign;
        if ( zsol != 0 ) delete zsol;
    }

    void clearall() {
        clearallbutmar();
        if ( mar!=0 ) delete mar;
    }
};

static PyObject *_4ti2GroebnerBasis( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    struct _4ti2GroebnerBasisInput data;

    int dim = -1;

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if (typeofinp.compare("mat") == 0) {
            if ( data.mat != 0 ) {
                data.clearall();
                PyErr_SetString(Py4ti2Error, "just one \'mat\' argument is possible");
                return NULL;
            }
            data.mat = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(eol2, *data.mat) ) {
                data.clearall();
                std::string pref = "\'mat\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
            if (dim >= 0 && dim != data.lat->get_size() ) {
                data.clearall();
                PyErr_SetString(Py4ti2Error, "size mismatch in matrix and lattice");
                return NULL;
            }
            else
                dim = data.mat->get_size();
        }
        else if (typeofinp.compare("lat") == 0) {
            data.lat = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(eol2, *data.lat) ) {
                data.clearall();
                std::string pref = "\'lat\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
            if (dim >= 0 && dim != data.mat->get_size() ) {
                data.clearall();
                PyErr_SetString(Py4ti2Error, "size mismatch in matrix and lattice");
                return NULL;
            }
            else
                dim = data.lat->get_size();
        }
        else if (typeofinp.compare("sign") == 0) {
            data.sign = new _4ti2_::Vector(dim);
            if ( !PyIntListToVector(eol2, *data.sign) ) {
                data.clearall();
                std::string pref = "\'sign\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("weights") == 0) {
            data.weights = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(eol2, *data.weights) ) {
                data.clearall();
                std::string pref = "\'weights\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("weightsmax") == 0) {
            data.weightsmax = new _4ti2_::Vector(dim);
            if ( !PyIntListToVector(eol2, *data.weightsmax) ) {
                data.clearall();
                std::string pref = "\'weightsmax\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("zsol") == 0) {
            data.zsol = new _4ti2_::Vector(dim);
            if ( !PyIntListToVector(eol2, *data.zsol) ) {
                data.clearall();
                std::string pref = "\'zsol\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("mar") == 0) {
            data.mar = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(eol2, *data.mar) ) {
                data.clearall();
                std::string pref = "\'mar\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("cost") == 0) {
            data.cost = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(eol2, *data.cost) ) {
                data.clearall();
                std::string pref = "\'cost\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else {
            data.clearall();
            PyErr_SetString(Py4ti2Error, "Unexpected argument");
            return NULL;
        }
    }
    // Tests input data according to input_Feasible.cpp 
    if ( data.mat == 0 && data.lat == 0 ) {
        data.clearall();
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }

    _4ti2_::BitSet urs(dim);
    if ( data.sign != 0 ) {
        for (int i = 0; i < dim; ++i) {
            IntegerType value = (*data.sign)[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                data.clearall();
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                data.clearall();
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }

    if ( (data.weights != 0 && data.weightsmax == 0) || 
            (data.weights == 0 && data.weightsmax != 0) ) {
        data.clearall();
        PyErr_SetString(Py4ti2Error, 
                "weightsmax and weights are both required or none should be given");
        return NULL;
    }
    else if (data.weights != 0 && data.weightsmax == 0) {
        if ( dim != data.weights->get_size() || dim != data.weightsmax->get_size() ) {
            data.clearall();
            PyErr_SetString(Py4ti2Error, 
                    "weights input dimension and problem dimension are different");
            return NULL;
        }
    }

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.lat, 
            data.mat, &urs, data.zsol, data.weights, data.weightsmax );
   
    _4ti2_::Globals::minimal = false;
    _4ti2_::GeneratingSet gs( *feasible, data.mar );

    _4ti2_::GroebnerBasis gb(gs, data.cost);

    _4ti2_::VectorArray binomios = gb.get_groebner_basis();

    delete feasible;
    data.clearallbutmar();
    
    if ( binomios.get_number() > 0 )
        return VectorArrayToPyIntListList( binomios );
    else
        return PyList_New(0); 
}

static PyObject *_4ti2Graver( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }
    
    _4ti2_state *graver_api;
    _4ti2_matrix *matlat, *sign, *lb, *ub;
#ifdef _4ti2_INT32_
    graver_api = _4ti2_graver_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    graver_api = _4ti2_graver_create_state(_4ti2_PREC_INT_64);
#endif
     PySys_WriteStdout( "number of arguments %i ", nargs / 2);
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( graver_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if (typeofinp.compare("mat") == 0 || typeofinp.compare("lat") == 0) {
            if ( !PyIntListListTo4ti2matrix(eol2, graver_api, typeofinp.c_str(), &matlat) ) {
                PySys_WriteStdout( "matlat rows %d\n", _4ti2_matrix_get_num_rows( matlat ));
                _4ti2_state_delete( graver_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given mat or lat argument, next, a list of lists of integers is expected");
//                return NULL;
                std::string pref = "\'mat\' or \'lat\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("sign") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, graver_api, typeofinp.c_str(), &sign) ) {
                _4ti2_state_delete( graver_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given sign argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'sign\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("lb") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, graver_api, typeofinp.c_str(), &lb) ) {
                _4ti2_state_delete( graver_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given lb argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'lb\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("ub") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, graver_api, typeofinp.c_str(), &ub) ) {
                _4ti2_state_delete( graver_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given ub argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'ub\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else {
            _4ti2_state_delete( graver_api );
            PyErr_SetString(Py4ti2Error, "Unexpected argument");
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
    
    PyObject *result = PyTuple_New( 4 );

    int elresul = 0;

    _4ti2_matrix* zhom_matrix;
    _4ti2_state_get_matrix(graver_api, "zhom", &zhom_matrix);
    PyObject *zhom_list;
    if ( zhom_matrix != 0 ) {
       zhom_list = _4ti2matrixToPyIntListList( zhom_matrix );
       PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zhom") );
       PyTuple_SetItem( result, elresul++, zhom_list );
    }

    _4ti2_matrix* zfree_matrix;
    _4ti2_state_get_matrix(graver_api, "zfree", &zfree_matrix);
    PyObject *zfree_list = 0;
    if (zfree_matrix != 0) {
        zfree_list = _4ti2matrixToPyIntListList( zfree_matrix );
        PyTuple_SetItem( result, elresul++, PyUnicode_FromString("zfree") );
        PyTuple_SetItem( result, elresul, zfree_list );
    }
    
    _4ti2_state_delete( graver_api );

    return result;
}

static PyObject *_4ti2Hilbert( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }
    
    _4ti2_state *hilbert_api;
    _4ti2_matrix *matlat, *rel, *sign, *ub;
    char quiet[] = "-q";
#ifdef _4ti2_INT32_
    char prec[] = "--precision=32";
    hilbert_api = _4ti2_hilbert_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    char prec[] = "--precision=64";
    hilbert_api = _4ti2_hilbert_create_state(_4ti2_PREC_INT_64);
#endif
//     PySys_WriteStdout( "number of arguments %i ", nargs / 2);
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( hilbert_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if (typeofinp.compare("mat") == 0 || typeofinp.compare("lat") == 0) {
            if ( !PyIntListListTo4ti2matrix(eol2, hilbert_api, typeofinp.c_str(), &matlat) ) {
//                PySys_WriteStdout( "matlat rows %d\n", _4ti2_matrix_get_num_rows( matlat ));
                _4ti2_state_delete( hilbert_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given mat or lat argument, next, a list of lists of integers is expected");
//                return NULL;
                std::string pref = "\'mat\' or \'lat\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("rel") == 0) { 
            if ( !PyIntListTo4ti2matrix(eol2, hilbert_api, typeofinp.c_str(), &rel) ) {
                _4ti2_state_delete( hilbert_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given rel argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'rel\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("sign") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, hilbert_api, typeofinp.c_str(), &sign) ) {
                _4ti2_state_delete( hilbert_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given sign argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'sign\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("ub") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, hilbert_api, typeofinp.c_str(), &ub) ) {
                _4ti2_state_delete( hilbert_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given ub argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'ub\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else {
            _4ti2_state_delete( hilbert_api );
            PyErr_SetString(Py4ti2Error, "Unexpected argument");
            return NULL;
        }
    }
    

    char **argv = new char* [3];
    argv[1] = quiet;
    argv[2] = prec;

    if ( _4ti2_state_set_options( hilbert_api, 3, argv) != _4ti2_OK ) {
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
    delete[] argv;

    return result;
}

static PyObject *_4ti2Zsolve( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args);

    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "Incorrect arguments");
        return NULL;
    }
    
    _4ti2_state *zsolve_api;
    _4ti2_matrix *matlat, *rhs, *rel, *sign, *ub, *lb;
    char quiet[] = "-q";
#ifdef _4ti2_INT32_
    char prec[] = "--precision=32";
    zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_32);
#elif defined(_4ti2_INT64_)
    char prec[] = "--precision=64";
    zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_64);
#endif
//     PySys_WriteStdout( "number of arguments %i ", nargs / 2);
    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            _4ti2_state_delete( zsolve_api );
            PyErr_SetString(Py4ti2Error, "Incorrect arguments");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if (typeofinp.compare("mat") == 0 || typeofinp.compare("lat") == 0) {
            if ( !PyIntListListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &matlat) ) {
//                PySys_WriteStdout( "matlat rows %d\n", _4ti2_matrix_get_num_rows( matlat ));
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given mat or lat argument, next, a list of lists of integers is expected");
//                return NULL;
                std::string pref = "\'mat\' or \'lat\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("rhs") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &rhs) ) {
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given rhs argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'rhs\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("rel") == 0) { 
            if ( !PyIntListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &rel) ) {
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given rel argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'rel\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("sign") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &sign) ) {
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given sign argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'sign\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("lb") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &lb) ) {
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given lb argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'lb\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else if (typeofinp.compare("ub") == 0) {
            if ( !PyIntListTo4ti2matrix(eol2, zsolve_api, typeofinp.c_str(), &ub) ) {
                _4ti2_state_delete( zsolve_api );
//                PyErr_SetString(Py4ti2Error, 
//                    "given ub argument, next, a lists of integers is expected");
//                return NULL;
                std::string pref = "\'ub\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return NULL;
            }
        }
        else {
            _4ti2_state_delete( zsolve_api );
            PyErr_SetString(Py4ti2Error, "Unexpected argument");
            return NULL;
        }
    }
    

    char **argv = new char* [3];
    argv[1] = quiet;
    argv[2] = prec;

    if ( _4ti2_state_set_options( zsolve_api, 3, argv) != _4ti2_OK ) {
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
    delete[] argv;

    return result;
}

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
    {"solve", (PyCFunction)_4ti2ParticularSolution, METH_VARARGS, 
        "Computes a particular solution of a linear diophantine equations system" },
    {"groebner", (PyCFunction)_4ti2GroebnerBasis, METH_VARARGS, 
        "Computes a Groebner basis of the toric ideal of a matrix, or, more general, of the lattice ideal of a lattice." },
    {"hilbert", (PyCFunction)_4ti2Hilbert, METH_VARARGS,
        "Computes the Hilbert basis of a matrix or a given lattice"},
    {"graver", (PyCFunction)_4ti2Graver, METH_VARARGS,
        "Computes the Graver basis of a matrix or a given lattice"},
    {"zsolve", (PyCFunction)_4ti2Zsolve, METH_VARARGS,
        "Solves linear inequality and equation systems over the integers"},
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
        "Py4ti2sol",
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

PyMODINIT_FUNC PyInit_Py4ti2int32(void)

#else

extern "C" void initPy4ti2int32(void)
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

PyMODINIT_FUNC PyInit_Py4ti2int64(void)

#else

extern "C" void initPy4ti2int64(void)
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
