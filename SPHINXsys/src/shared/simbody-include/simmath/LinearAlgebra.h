#ifndef SimTK_LINEAR_ALGEBRA_H_
#define SimTK_LINEAR_ALGEBRA_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath linear algebra tools.
 */


#include "SimTKcommon.h"
#include "simmath/internal/common.h"


namespace SimTK {

//  default for reciprocal of the condition number
// TODO: sherm 080128 I changed this from 0.01 to a more reasonable
// value but it is still wrong because the default should depend
// on the matrix size, something like max(m,n)*eps^(7/8) where
// eps is machine precision for float or double as appropriate.
static const double DefaultRecpCondition = 1e-12;

/**
 * Class to compute Eigen values and Eigen vectors of a matrix
 */
class SimTK_SIMMATH_EXPORT Eigen {
    public:

    ~Eigen();

    Eigen();
    Eigen( const Eigen& c );
    Eigen& operator=(const Eigen& rhs);

    /// create a default eigen class
    template <class ELT> Eigen( const Matrix_<ELT>& m );
    /// supply matrix which eigen values will be computed for  
    template <class ELT> void factor( const Matrix_<ELT>& m );
    /// get all the eigen values and eigen vectors of a matrix 
    template <class VAL, class VEC> void getAllEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors);
    /// get all the eigen values of a matrix 
    template <class T> void getAllEigenValues( Vector_<T>& values);

    /// get a few eigen values  and eigen vectors of a symmetric matrix which are within a range of indices
    template <class VAL, class VEC> void getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, int ilow, int ihi);
    /// get a few eigen vectors of a symmetric matrix which are within a range of indices
    template <class T> void getFewEigenVectors( Matrix_<T>& vectors, int ilow, int ihi );
    /// get a few eigen values of a symmetric matrix which are within a range of indices
    template <class T> void getFewEigenValues( Vector_<T>& values, int ilow, int ihi );

    /// get a few eigen values  and eigen vectors of a symmetric matrix which are within a range of eigen values
    template <class VAL, class VEC> void getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, typename CNT<VAL>::TReal rlow, typename CNT<VAL>::TReal rhi);
    /// get a few eigen vectors of a symmetric matrix which are within a range of eigen values
    template <class T> void getFewEigenVectors( Matrix_<T>& vectors, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi );
    /// get a few eigen values of a symmetric matrix which are within a range of eigen values
    template <class T> void getFewEigenValues( Vector_<T>& values, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi );

     
    protected:
    class EigenRepBase *rep;

}; // class Eigen

} // namespace SimTK 

#endif //SimTK_LINEAR_ALGEBRA_H_
