
#ifndef SimTK_SIMMATH_LAPACK_INTERFACE_H_
#define SimTK_SIMMATH_LAPACK_INTERFACE_H_

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

namespace SimTK {

class LapackInterface { 
   
public:

static int getLWork( float* work);
static int getLWork( double* work);
static int getLWork( std::complex<float>* work);
static int getLWork( std::complex<double>* work);

template <class T> static
void geev(char jobvl, char jobvr, int n, T* a, int lda, 
    std::complex<typename CNT<T>::TReal>* values, 
    T* vl, int ldvl, std::complex<typename CNT<T>::TReal>* vr, 
    int ldvr, T* work, int lwork, int& info );

}; // class LapackInterface

}   // namespace SimTK

#endif // SimTK_SIMMATH_LAPACK_INTERFACE_H_
