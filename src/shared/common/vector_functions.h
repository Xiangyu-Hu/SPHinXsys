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

inline Vec3d upgradeToVec3d(const Real &input)
{
    return Vec3d(input, 0.0, 0.0);
};
inline Vec3d upgradeToVec3d(const Vec2d &input)
{
    return Vec3d(input[0], input[1], 0.0);
};
inline Vec3d upgradeToVec3d(const Vec3d &input) { return input; };

inline Mat3d upgradeToMat3d(const Mat2d &input)
{
    Mat3d output = Mat3d::Zero();
    output.block<2, 2>(0, 0) = input;
    return output;
};
inline Mat3d upgradeToMat3d(const Mat3d &input) { return input; };

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

inline Real first_component(const Real &input) { return input; };
template <int Dim1, int Dim2>
Real first_component(const Eigen::Matrix<Real, Dim1, Dim2> &input) { return input(0, 0); };

inline Real component_square(const Real &input) { return input * input; };

template <int Dim1, int Dim2>
Eigen::Matrix<Real, Dim1, Dim2> component_square(const Eigen::Matrix<Real, Dim1, Dim2> &input) { return input.cwiseAbs2(); };

template <typename ComponentFunction, typename... Args>
Real transform_component(const Real &input, const ComponentFunction &function, Args &&...args)
{
    return function(input, std::forward<Args>(args)...);
};

template <int Dim1, int Dim2, typename ComponentFunction, typename... Args>
Eigen::Matrix<Real, Dim1, Dim2> transform_component(
    const Eigen::Matrix<Real, Dim1, Dim2> &input, const ComponentFunction &function, Args &&...args)
{
    Eigen::Matrix<Real, Dim1, Dim2> output;
    for (int i = 0; i < Dim1; ++i)
        for (int j = 0; j < Dim2; ++j)
            output(i, j) = function(input(i, j), std::forward<Args>(args)(i, j)...);
    return output;
};

template <typename ComponentFunction, typename... Args>
void for_each_component(const Real &input, const ComponentFunction &function, Args &&...args)
{
    function(input, std::forward<Args>(args)...);
};

template <int Dim1, int Dim2, typename ComponentFunction, typename... Args>
void for_each_component(const Eigen::Matrix<Real, Dim1, Dim2> &input,
                        const ComponentFunction &function, Args &&...args)
{
    Eigen::Matrix<Real, Dim1, Dim2> output;
    for (int i = 0; i < Dim1; ++i)
        for (int j = 0; j < Dim2; ++j)
            function(input(i, j), std::forward<Args>(args)(i, j)...);
};

inline Vec3d vectorizeTensorSquare(const Vec2d &input)
{
    return Vec3d(input[0] * input[0], input[1] * input[1], input[0] * input[1]);
};

inline Vec6d vectorizeTensorSquare(const Vec3d &input)
{
    return Vec6d(input[0] * input[0], input[1] * input[1], input[2] * input[2],
                 input[0] * input[1], input[1] * input[2], input[2] * input[0]);
};

inline Vec3d vectorizeSymMatrix(const Mat2d &input)
{
    return Vec3d(input(0, 0), input(1, 1), input(0, 1));
};

inline Vec6d vectorizeSymMatrix(const Mat3d &input)
{
    return Vec6d(input(0, 0), input(1, 1), input(2, 2), input(0, 1), input(1, 2), input(2, 0));
};

template <typename DataType>
DataType tensorProduct(const DataType &value1, const Real &value2)
{
    return value1 * value2;
};

template <int N, int M, int O>
Eigen::Matrix<Real, N, M> tensorProduct(const Eigen::Matrix<Real, N, O> &value1,
                                        const Eigen::Matrix<Real, M, O> &value2)
{
    return value1 * value2.transpose();
};
} // namespace SPH
#endif // VECTOR_FUNCTIONS_H
