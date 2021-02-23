/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**@file
 * These is a templatized, C++ callable interface to LAPACK and BLAS.
 * Each method must be explicitly specialized for the supported precisions.
 */


#include "SimTKcommon.h"
#include "SimTKlapack.h"
#include "SimTKmath.h"
#include "LapackInterface.h"
#include "WorkSpace.h"
#include <cstring>

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strcat, sprintf, etc.
#endif

static const double EPS = .000001;
namespace SimTK {

int LapackInterface::getLWork( float* work) { return( (int)work[0] ); }
int LapackInterface::getLWork( double* work) { return( (int)work[0] ); }
int LapackInterface::getLWork( std::complex<float>* work) { return( (int)work[0].real() ); }
int LapackInterface::getLWork( std::complex<double>* work) { return( (int)work[0].real() ); }

// eigenvlaues for nonsymmetric matrices
template <class T> 
void LapackInterface::geev (char jobvl, char jobvr,
    int n, T* a, int lda, std::complex<typename CNT<T>::TReal>* values, 
    T* vl, int ldvl, std::complex<typename CNT<T>::TReal>* vr, 
    int ldvr, T* work, int lwork, int& info ) {
    assert(false);
}

template <> void LapackInterface::geev<double>
   (char jobvl, char jobvr,
    int n, double* a, int lda, std::complex<double>* values, 
    double* vl, int ldvl, std::complex<double>* rightVectors, int ldvr, double* work,
    int lwork, int& info )
{
    TypedWorkSpace<double> wr(n);
    TypedWorkSpace<double> wi(n);
    TypedWorkSpace<double> vr(n*n);

    // avoid valgrind uninitialized warnings
    for(int i=0;i<n;i++) wi.data[i] = 0;  
    dgeev_( jobvl, jobvr, 
//             n, a, lda, wr.data, wi.data, vl, ldvl, vr.data, ldvr, 
             n, a, lda, wr.data, wi.data, vl, ldvl, vr.data, ldvr, 
             work, lwork, info, 1, 1);

    if( info < 0 ) {
        SimTK_THROW2( SimTK::Exception::IllegalLapackArg, "dgeev", info );
    }

    for(int i=0;i<n;i++) {
        values[i] = std::complex<double>(wr.data[i], wi.data[i] );
    }

    /*
    ** LAPACK returns the eigen vectors as complex conjuate pairs 
    ** if the eigen value is real  ( imaginary part == 0 ) then the eigen vector is real
    ** else the vectors are returned with the real part in the jth column and the
    ** imaginary part in the j+1 column
    */
    for(int j=0;j<n;j++) {
        if( fabs(wi.data[j]) < EPS ) {
            for(int i=0;i<n;i++) {
                rightVectors[j*n+i] = std::complex<double>(vr.data[j*n+i], 0.0 );
             }
        } else {
            for(int i=0;i<n;i++) {
                rightVectors[j*n+i] = std::complex<double>(vr.data[j*n+i], vr.data[(j+1)*n+i]);
                rightVectors[(j+1)*n+i] = std::complex<double>(vr.data[j*n+i], -vr.data[(j+1)*n+i]);
            }
            j++;
        }
    } 
/*
    for(int j=0;j<n;j++) { 
        for(int i=0;i<n;i++) printf("%f %f    ", rightVectors(i,j).real(), rightVectors(i,j).imag() ); 
        printf("\n");
    }
*/
}

template <> void LapackInterface::geev<float>
   (char jobvl, char jobvr,
    int n, float* a, int lda, std::complex<float>* values,
    float* vl, int ldvl, std::complex<float>* rightVectors, int ldvr, float* work,
    int lwork, int& info )
{
    TypedWorkSpace<float> wr(n);
    TypedWorkSpace<float> wi(n);
    TypedWorkSpace<float> vr(n*n);

    // avoid valgrind uninitialized warnings
    for(int i=0;i<n;i++) wi.data[i] = 0;  

    sgeev_( jobvl, jobvr, 
             n, a, lda, wr.data, wi.data, vl, ldvl, vr.data, ldvr, 
             work, lwork, info, 
             1, 1);

    if( info < 0 ) {
        SimTK_THROW2( SimTK::Exception::IllegalLapackArg, "sgeev", info );
    }

    for(int i=0;i<n;i++) {
        values[i] = std::complex<float>(wr.data[i], wi.data[i] );
    }
    /*
    ** LAPACK returns the eigen vectors as complex conjuate pairs 
    ** if the eigen value is real  ( imaginary part == 0 ) then the eigen vector is real
    ** else the vectors are returned with the real part in the jth column and the
    ** imaginary part in the j+1 column
    */
    for(int j=0;j<n;j++) {
        if( fabs(wi.data[j]) < (float)EPS ) {
            for(int i=0;i<n;i++) {
                rightVectors[j*n+i] = std::complex<float>(vr.data[j*n+i], 0.0 );
//printf(" %f ",vr.data[j*n+i] );
             }
//printf("\n");
        } else {
            for(int i=0;i<n;i++) {
                rightVectors[j*n+i] = std::complex<float>(vr.data[j*n+i], vr.data[(j+1)*n+i]);
                rightVectors[(j+1)*n+i] = std::complex<float>(vr.data[j*n+i], -vr.data[(j+1)*n+i]);
            }
            j++;
    }
    } 
}
template <> void LapackInterface::geev<std::complex<float> >
   (char jobvl, char jobvr,
    int n, std::complex<float>* a, int lda, std::complex<float>* values, 
    std::complex<float>* vl, int ldvl, std::complex<float>* rightVectors, int ldvr, std::complex<float>* work,
    int lwork, int& info )
{

    TypedWorkSpace<float> Rwork(2*n);
    cgeev_( jobvl, jobvr, 
             n, a, lda, values,  vl, ldvl, rightVectors, ldvr, 
             work, lwork, Rwork.data, info, 
             1, 1);

    if( info < 0 ) {
        SimTK_THROW2( SimTK::Exception::IllegalLapackArg, "cgeev", info );
    }

}

template <> void LapackInterface::geev<std::complex<double> >
   (char jobvl, char jobvr,
    int n, std::complex<double>* a, int lda, std::complex<double>* values, 
    std::complex<double>* vl, int ldvl, std::complex<double>* rightVectors, int ldvr, std::complex<double>* work,
    int lwork, int& info )
{

    TypedWorkSpace<double> Rwork(2*n);
    zgeev_( jobvl, jobvr, 
             n, a, lda, values,  vl, ldvl, rightVectors, ldvr, 
             work, lwork, Rwork.data, info, 
             1, 1);

    if( info < 0 ) {
        SimTK_THROW2( SimTK::Exception::IllegalLapackArg, "zgeev", info );
    }

}


}   // namespace SimTK

