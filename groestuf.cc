#include "datatran.h"
#include "vecarcnv.h"

#include <Python.h>

#include <iostream>
#include <sstream>
#include <streambuf>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <stdexcept>

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
#include "groebner/Optimise.h"
#include "groebner/WalkAlgorithm.h"

#include "glpk.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

#if PY_MAJOR_VERSION >= 3
#define string_check PyUnicode_Check
#else
#define string_check PyString_Check
#endif

extern PyObject * Py4ti2Error;

#ifndef _4ti2_GMP_

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
#endif

#ifndef _4ti2_HAVE_GMP

PyObject *_4ti2ParticularSolution( PyObject *self, PyObject *args )
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

#endif

struct VecAr_InputData {
    int dim;
    std::map<std::string, int> str_va_map, str_v_map;
    _4ti2_::VectorArray *va_mar; // *mar for the special case of Markov input 
                                 // for GeneratingSet
    _4ti2_::VectorArray **va_data; // *mat, *lat, *sign, *gro, *cost, 
                                   // *feas, *weigths, *gro_start, *cost_start;
    _4ti2_::Vector **v_data; // *sign, *weightsmax, *zsol;

    // n is the length of x
    VecAr_InputData( int n, const char *va[], int m, const char *v[] ) {
        assert( n >= 0 && m >= 0 );
        dim = -1;
        va_mar = 0;
        va_data = new _4ti2_::VectorArray *[n];
        for (int i = 0; i < n ; ++i ) {
            str_va_map[std::string(va[i])] = i;
            va_data[i] = 0;
        }
        v_data = new _4ti2_::Vector *[m];
        for (int i = 0; i < m ; ++i ) {
            str_v_map[std::string(v[i])] = i;
            v_data[i] = 0;
        }
    }

    void delete_mar() {
        if ( va_mar != 0 ) delete va_mar;
    }

    ~VecAr_InputData() {
        int lim = str_va_map.size();
        for ( int i = 0; i < lim; ++i )
            if (va_data[i] != 0) delete va_data[i];
        
        lim = str_v_map.size();
        for ( int i = 0; i < lim; ++i )
            if (v_data[i] != 0)  delete v_data[i];
    }

    int v_index( const char *i_t ) {
        int ret;
        try {
            std::string strit(i_t);
            ret = str_v_map.at(strit);
        }
        catch (std::out_of_range& e) {
            return -1;
        }
        return ret;
    }

    int va_index( const char *i_t ) {
        int ret;
        try {
            std::string strit(i_t);
            ret = str_va_map.at(strit);
        }
        catch (std::out_of_range& e) {
            return -1;
        }
        return ret;
    }

    bool read( std::string& input_type, PyObject *pydata ) {
        // mar input matrix: weird processing due to _4ti2_::GeneratingSet
        if (input_type.compare("mar") == 0) {
            va_mar = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(pydata, *va_mar) ) {
                std::string pref = "\'mar\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return false;
            }
            if (dim >= 0 && dim != va_mar->get_size() ) {
                whathappened = "\'mar\' argument: data input size mismatch";
                return false;
            }
            else if (dim <= 0) dim = va_mar->get_size();
            return true;
        }
        std::map<std::string, int>::iterator itva = str_va_map.find(input_type);
        int idx;
        if ( itva != str_va_map.end()) {
            idx = itva->second;
            
            // By default, we destroy repeated inputs readed yet
            if (va_data[idx] != 0) delete va_data[idx];
            va_data[idx] = new _4ti2_::VectorArray();
            if ( !PyIntListListToVectorArray(pydata, *(va_data[idx])) ) {
                std::string pref = "\'"+ input_type + "\' argument: ";
                whathappened = pref + whathappened;
                PyErr_SetString(Py4ti2Error, whathappened.c_str());
                return false;
            }
            if (dim >= 0 && dim != va_data[idx]->get_size() ) {
                whathappened = "\'" + input_type 
                                + "\' argument: data input size mismatch";
                return false;
            }
            else if (dim <= 0) 
                dim = va_data[idx]->get_size();
        }
        else {
            itva = str_v_map.find(input_type);
            if ( itva != str_v_map.end() ) {
                if (dim <= 0) {
                    if (!PyList_Check(pydata)) return false;
                    else dim = PyList_Size(pydata);
                }
                else if (dim != PyList_Size(pydata)) {
                    whathappened = "\'" + input_type 
                                + "\' argument: data input size mismatch";
                    return false; 
                }
                idx = itva->second;
                if (v_data[idx] != 0) delete v_data[idx];
                v_data[idx] = new _4ti2_::Vector(dim);
                if ( !PyIntListToVector(pydata, *(v_data[idx])) ) {
                    std::string pref = "\'" + input_type + "\' argument: ";
                    whathappened = pref + whathappened;
                    return false;
                }
            }
            else {
                whathappened = "\'" + input_type + "\' argument: unknown data input";
                return false;
            }
        }
        return true;
    }
};

PyObject *_4ti2Minimize( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    const char *va_types[] = {"mat", "lat", "cost"};
    const char *v_types[] = {"sign", "zsol"};
    struct VecAr_InputData data(3, va_types, 2, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, 
                "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }
    // Tests input data according to input_Feasible.cpp 
    int idx_mat = data.va_index("mat");
    int idx_lat = data.va_index("lat");
    if ( data.va_data[idx_mat] == 0 && data.va_data[idx_lat] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }

    _4ti2_::BitSet urs(data.dim);
    int idx_sign = data.v_index("sign");
    if ( data.v_data[idx_sign] != 0 ) {
        for (int i = 0; i < data.dim; ++i) {
            IntegerType value = (*data.v_data[idx_sign])[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }

    // Test input according to minimize_main.cpp
    int idx_cost = data.va_index("cost");
    if ( data.va_data[idx_cost] == 0 || 
         data.va_data[idx_cost]->get_number() != 1 ) {
        PyErr_SetString(Py4ti2Error, 
            "there should be a single cost function");
        return NULL;
    }

    int idx_zsol = data.v_index("zsol");
    if ( data.va_data[idx_zsol] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a fiber (zsol) is needed as input data");
        return NULL;
    }

    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.va_data[idx_lat], 
            data.va_data[idx_mat], &urs, data.v_data[idx_zsol] );

    _4ti2_::Vector sol(*feasible->get_rhs());
    _4ti2_::Optimise opt;
    opt.compute( *feasible, (*(data.va_data[idx_cost]))[0], sol );

    // Restore normal standard output
    std::cout.rdbuf( old );

    delete feasible;
    
    return VectorToPyIntList( sol );
}

PyObject *_4ti2GroebnerBasis( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    const char *va_types[] = {"mat", "lat", "cost", "weights", "mar"};
    const char *v_types[] = {"sign", "zsol", "weightsmax"};
    struct VecAr_InputData data(5, va_types, 3, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            data.delete_mar();
            PyErr_SetString(Py4ti2Error, 
                    "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            data.delete_mar();
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    // Tests input data according to input_Feasible.cpp 
    int idx_mat = data.va_index("mat");
    int idx_lat = data.va_index("lat");
    if ( data.va_data[idx_mat] == 0 && data.va_data[idx_lat] == 0 ) {
        data.delete_mar();
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }
    
    _4ti2_::BitSet urs(data.dim);
    int idx_sign = data.v_index("sign");

    if ( data.v_data[idx_sign] != 0 ) {
        for (int i = 0; i < data.dim; ++i) {
            IntegerType value = (*data.v_data[idx_sign])[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                data.delete_mar();
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                data.delete_mar();
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }
    
    // Test input data according to function specs
    int idx_weights = data.va_index("weights");
    int idx_weightsmax = data.v_index("weightsmax");

    if ( (data.va_data[idx_weights] != 0 && 
          data.v_data[idx_weightsmax] == 0) || 
            (data.va_data[idx_weights] == 0 && 
             data.v_data[idx_weightsmax] != 0) ) {
        data.delete_mar();
        PyErr_SetString(Py4ti2Error, 
          "weightsmax and weights are both required or none should be given");
        return NULL;
    }

    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.va_data[idx_lat], 
            data.va_data[idx_mat], &urs, data.v_data[data.v_index("zsol")],
            data.va_data[idx_weights], data.v_data[idx_weightsmax] );

    _4ti2_::Globals::minimal = false;
    _4ti2_::GeneratingSet gs( *feasible, data.va_mar );

    _4ti2_::GroebnerBasis gb(gs, data.va_data[data.va_index("cost")]);

    _4ti2_::VectorArray gbe = gb.get_groebner_basis();

    // Restore normal standard output
    std::cout.rdbuf( old );

    delete feasible;
    
    if ( gbe.get_number() > 0 )
        return VectorArrayToPyIntListList( gbe );
    else
        return PyList_New(0); 
}

PyObject *_4ti2NormalForm( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }

    const char *va_types[] = {"mat", "lat", "gro", "cost", "feas"};
    const char *v_types[] = {"sign"};
    struct VecAr_InputData data(5, va_types, 1, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, 
                    "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    // Tests input data according to input_Feasible.cpp 
    int idx_mat = data.va_index("mat");
    int idx_lat = data.va_index("lat");
    if ( data.va_data[idx_mat] == 0 && data.va_data[idx_lat] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }

    _4ti2_::BitSet urs(data.dim);
    int idx_sign = data.v_index("sign");
    if ( data.v_data[idx_sign] != 0 ) {
        for (int i = 0; i < data.dim; ++i) {
            IntegerType value = (*data.v_data[idx_sign])[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }

    // Test input as in normalform_main.cpp
    int idx_gro = data.va_index("gro");
    if ( data.va_data[idx_gro] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "there should be a Groebner basis input");
        return NULL;
    }

    int idx_feas = data.va_index("feas");
    if ( data.va_data[idx_feas] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a list of feasible solutions is needed");
        return NULL;
    }

    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.va_data[idx_lat], 
            data.va_data[idx_mat], &urs );

    int idx_cost = data.va_index("cost");
    if ( data.va_data[idx_cost] == 0 ) 
        data.va_data[idx_cost] = 
            new _4ti2_::VectorArray(0, feasible->get_dimension());

    _4ti2_::Minimize algorithm;
    algorithm.minimize( *feasible, *(data.va_data[idx_cost]), 
                        *(data.va_data[idx_gro]), *(data.va_data[idx_feas]) );

    // Restore normal standard output
    std::cout.rdbuf( old );

    delete feasible;

    PyObject *nf_list = VectorArrayToPyIntListList(*(data.va_data[idx_feas]));
    
    return nf_list;
}

PyObject *_4ti2MarkovBasis( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    const char *va_types[] = {"mat", "lat", "weights"};
    const char *v_types[] = {"sign", "zsol", "weightsmax"};
    struct VecAr_InputData data(3, va_types, 3, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, 
                    "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    // Tests input data according to input_Feasible.cpp 
    int idx_mat = data.va_index("mat");
    int idx_lat = data.va_index("lat");
    if ( data.va_data[idx_mat] == 0 && data.va_data[idx_lat] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }
    
    _4ti2_::BitSet urs(data.dim);
    int idx_sign = data.v_index("sign");

    if ( data.v_data[idx_sign] != 0 ) {
        for (int i = 0; i < data.dim; ++i) {
            IntegerType value = (*data.v_data[idx_sign])[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }

    // Test input data according to function specs
    int idx_weights = data.va_index("weights");
    int idx_weightsmax = data.v_index("weightsmax");

    if ( (data.va_data[idx_weights] != 0 && 
          data.v_data[idx_weightsmax] == 0) || 
            (data.va_data[idx_weights] == 0 && 
             data.v_data[idx_weightsmax] != 0) ) {
        PyErr_SetString(Py4ti2Error, 
          "weightsmax and weights are both required or none should be given");
        return NULL;
    }

    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.va_data[idx_lat], 
            data.va_data[idx_mat], &urs, data.v_data[data.v_index("zsol")],
            data.va_data[idx_weights], data.v_data[idx_weightsmax] );

    _4ti2_::GeneratingSet gs( *feasible, 0 );
    gs.standardise();

    _4ti2_::VectorArray mar = gs.get_generating_set();

    // Restore normal standard output
    std::cout.rdbuf( old );

    delete feasible;
    
    if ( mar.get_number() > 0 )
        return VectorArrayToPyIntListList( mar );
    else
        return PyList_New(0); 
}

PyObject *_4ti2ZBasis( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    const char *va_types[] = {"mat"};
    const char *v_types[] = {};
    struct VecAr_InputData data(1, va_types, 0, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    int idx_mat = data.va_index("mat");
    if ( data.va_data[idx_mat] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a matrix is needed as input data");
        return NULL;
    }
    
    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::VectorArray zbasis(0, data.va_data[idx_mat]->get_size());
    lattice_basis( *(data.va_data[idx_mat]), zbasis );
    
    // Restore normal standard output
    std::cout.rdbuf( old );

    if ( zbasis.get_number() > 0 )
        return VectorArrayToPyIntListList( zbasis );
    else
        return PyList_New(0); 
}


PyObject *_4ti2Walk( PyObject *self, PyObject *args )
{
    const int nargs = PyTuple_Size(args); 
    
    if ( nargs % 2 != 0 ) {
        PyErr_SetString(Py4ti2Error, "incorrect arguments: an even number of arguments is expected");
        return NULL;
    }
    const char *va_types[] = {"mat", "lat", "cost", "gro.start", "cost.start"};
    const char *v_types[] = {"sign", "zsol"};
    struct VecAr_InputData data(5, va_types, 2, v_types);

    for (int i = 0; i < nargs / 2; ++i) {
        PyObject *eol1 = PyTuple_GetItem(args, 2*i);
        if ( !string_check(eol1) ) {
            PyErr_SetString(Py4ti2Error, "incorrect arguments: odd arguments should be strings");
            return NULL;
        }
        std::string typeofinp = PyUnicodeToString(eol1);
        
        PyObject *eol2 = PyTuple_GetItem(args, 2*i+1);
        if ( !data.read(typeofinp, eol2) ) {
            PyErr_SetString(Py4ti2Error, whathappened.c_str());
            return NULL;
        }
    }

    // Tests input data according to input_Feasible.cpp 
    int idx_mat = data.va_index("mat");
    int idx_lat = data.va_index("lat");
    if ( data.va_data[idx_mat] == 0 && data.va_data[idx_lat] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a matrix and/or lattice is needed as input data");
        return NULL;
    }
    
    _4ti2_::BitSet urs(data.dim);
    int idx_sign = data.v_index("sign");

    if ( data.v_data[idx_sign] != 0 ) {
        for (int i = 0; i < data.dim; ++i) {
            IntegerType value = (*data.v_data[idx_sign])[i];
            if (value == 0) urs.set(i); 
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1) {
                data.delete_mar();
                PyErr_SetString(Py4ti2Error, 
                    "some value in sign is not yet supported");
                return NULL;
            }
            else {
                data.delete_mar();
                PyErr_SetString(Py4ti2Error, 
                    "unsupported number in sign vector");
                return NULL;
            }
        }
    }
    
    // Tests input data according to walk_main.cpp 
    int idx_grostart = data.va_index("gro.start");
    if ( data.va_data[idx_grostart] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "the starting Groebner basis is needed as input data");
        return NULL;
    }

    int idx_cost = data.va_index("cost");
    if ( data.va_data[idx_cost] == 0 ) {
        PyErr_SetString(Py4ti2Error, 
            "a target cost is needed as input data");
        return NULL;
    }

    int idx_coststart = data.va_index("cost.start");
    if ( data.va_data[idx_coststart] == 0 ) {
       data.va_data[idx_coststart] = new _4ti2_::VectorArray(0, data.dim);
    }

    // Suppress GLPK/4ti2 output
    glp_term_out( GLP_OFF );  
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    std::cout.rdbuf( ss.rdbuf() );

    _4ti2_::Feasible *feasible = new _4ti2_::Feasible( data.va_data[idx_lat], 
            data.va_data[idx_mat], &urs, data.v_data[data.v_index("zsol")] );

    _4ti2_::WalkAlgorithm algorithm;
    algorithm.compute( *feasible, *(data.va_data[idx_coststart]), 
                        *(data.va_data[idx_grostart]), 
                        *(data.va_data[idx_cost]) );

    // Restore normal standard output
    std::cout.rdbuf( old );

    delete feasible;
    
    if ( data.va_data[idx_grostart]->get_number() > 0 )
        return VectorArrayToPyIntListList( *(data.va_data[idx_grostart]) );
    else
        return PyList_New(0); 
}
