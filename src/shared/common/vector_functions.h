/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	vector_functions.h
 * @brief 	Basic functions for vector type data.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef VECTOR_FUNCTIONS_H
#define VECTOR_FUNCTIONS_H

#include "data_type.h"

namespace SPH
{
Vec2d FirstAxisVector(const Vec2d &zero_vector);
Vec3d FirstAxisVector(const Vec3d &zero_vector);

Vec3d upgradeToVec3d(const Real &input);
Vec3d upgradeToVec3d(const Vec2d &input);
Vec3d upgradeToVec3d(const Vec3d &input);
Mat3d upgradeToMat3d(const Mat2d &input);
Mat3d upgradeToMat3d(const Mat3d &input);

Vecd degradeToVecd(const Vec3d &input);
Matd degradeToMatd(const Mat3d &input);

Mat2d getInverse(const Mat2d &A);
Mat3d getInverse(const Mat3d &A);
Mat2d getAverageValue(const Mat2d &A, const Mat2d &B);
Mat3d getAverageValue(const Mat3d &A, const Mat3d &B);
Mat2d inverseCholeskyDecomposition(const Mat2d &A);
Mat3d inverseCholeskyDecomposition(const Mat3d &A);
Mat2d getDiagonal(const Mat2d &A);
Mat3d getDiagonal(const Mat3d &A);

/** Real dot product between two matrices, resulting in a scalar value (sum of products of element-wise) */
Real CalculateBiDotProduct(Mat2d Matrix1, Mat2d Matrix2); // calculate Real dot
Real CalculateBiDotProduct(Mat3d Matrix1, Mat3d Matrix2); // calculate Real dot

/** get transformation matrix. */
Mat2d getTransformationMatrix(const Vec2d &direction_of_y);
Mat3d getTransformationMatrix(const Vec3d &direction_of_z);
Mat3d getTransformationMatrix(const Vec3d &direction_of_z, const Vec3d &direction_of_y);

template <typename VecType>
Real getCosineOfAngleBetweenTwoVectors(const VecType &vector_1, const VecType &vector_2)
{
    return vector_1.dot(vector_2) / (vector_1.norm() * vector_2.norm() + TinyReal);
};

/** get the projection of the vector_1 on vector 2,
 *  which is parallel to the vector_2, meaning it is the vector_2 * scalar */
template <typename VecType>
VecType getVectorProjectionOfVector(const VecType &vector_1, const VecType &vector_2)
{
    return vector_1.dot(vector_2) * vector_2 / (vector_2.squaredNorm() + TinyReal);
};

template <typename Datatype>
Real getSquaredNorm(const Datatype &variable) { return variable.squaredNorm(); };

inline Real getSquaredNorm(const Real &variable) { return variable * variable; };

/** von Mises stress from stress matrix */
Real getVonMisesStressFromMatrix(const Mat2d &sigma);
Real getVonMisesStressFromMatrix(const Mat3d &sigma);

/** principal strain or stress from strain or stress matrix */
Vec2d getPrincipalValuesFromMatrix(const Mat2d &A);
Vec3d getPrincipalValuesFromMatrix(const Mat3d &A);

/** get transformation matrix. */
Real getCrossProduct(const Vec2d &vector_1, const Vec2d &vector_2);
Vec3d getCrossProduct(const Vec3d &vector_1, const Vec3d &vector_2);
/** Modulo operation for Arrayi */
Array2i mod(const Array2i &input, int modulus);
Array3i mod(const Array3i &input, int modulus);
} // namespace SPH
#endif // VECTOR_FUNCTIONS_H
