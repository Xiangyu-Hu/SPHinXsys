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
#ifndef SMALL_VECTORS_H
#define SMALL_VECTORS_H

#include "base_data_type.h"

namespace SPH
{
Vec2d FirstAxisVector(const Vec2d &zero_vector);
Vec3d FirstAxisVector(const Vec3d &zero_vector);

Vec3d upgradeToVec3d(const Real &input);
Vec3d upgradeToVec3d(const Vec2d &input);
Vec3d upgradeToVec3d(const Vec3d &input);
Mat3d upgradeToMat3d(const Mat2d &input);
Mat3d upgradeToMat3d(const Mat3d &input);

void degradeToVecd(const Vec3d &input, Vec2d &output);
inline void degradeToVecd(const Vec3d &input, Vec3d &output) { output = input; };
void degradeToMatd(const Mat3d &input, Mat2d &output);
inline void degradeToMatd(const Mat3d &input, Mat3d &output) { output = input; };

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
/** get angle between two vectors. */
Real getCosineOfAngleBetweenTwoVectors(const Vec2d &vector_1, const Vec2d &vector_2);
Real getCosineOfAngleBetweenTwoVectors(const Vec3d &vector_1, const Vec3d &vector_2);

/** get orthogonal projection of a vector. */
Vec2d getVectorProjectionOfVector(const Vec2d &vector_1, const Vec2d &vector_2);
Vec3d getVectorProjectionOfVector(const Vec3d &vector_1, const Vec3d &vector_2);

/** von Mises stress from stress matrix */
Real getVonMisesStressFromMatrix(const Mat2d &sigma);
Real getVonMisesStressFromMatrix(const Mat3d &sigma);

/** principal strain or stress from strain or stress matrix */
Vec2d getPrincipalValuesFromMatrix(const Mat2d &A);
Vec3d getPrincipalValuesFromMatrix(const Mat3d &A);

/** get transformation matrix. */
Real getCrossProduct(const Vec2d &vector_1, const Vec2d &vector_2);
Vec3d getCrossProduct(const Vec3d &vector_1, const Vec3d &vector_2);

/** Bounding box for system, body, body part and shape, first: lower bound, second: upper bound. */
template <typename VecType>
class BaseBoundingBox
{
  public:
    VecType first_, second_;
    int dimension_;

    BaseBoundingBox() : first_(VecType::Zero()), second_(VecType::Zero()), dimension_(VecType::Zero().size()){};
    BaseBoundingBox(const VecType &lower_bound, const VecType &upper_bound)
        : first_(lower_bound), second_(upper_bound), dimension_(lower_bound.size()){};
    /** Check the bounding box contain. */
    bool checkContain(const VecType &point)
    {
        bool is_contain = true;
        for (int i = 0; i < dimension_; ++i)
        {
            if (point[i] < first_[i] || point[i] > second_[i])
            {
                is_contain = false;
                break;
            }
        }
        return is_contain;
    };
};
/** Operator define. */
template <class T>
bool operator==(const BaseBoundingBox<T> &bb1, const BaseBoundingBox<T> &bb2)
{
    return bb1.first_ == bb2.first_ && bb1.second_ == bb2.second_ ? true : false;
};
/** Intersection fo bounding box.*/
template <class BoundingBoxType>
BoundingBoxType getIntersectionOfBoundingBoxes(const BoundingBoxType &bb1, const BoundingBoxType &bb2)
{
    /** Check that the inputs are correct. */
    int dimension = bb1.dimension_;
    /** Get the Bounding Box of the intersection of the two meshes. */
    BoundingBoxType bb(bb1);
    /** #1 check that there is overlap, if not, exception. */
    for (int i = 0; i < dimension; ++i)
        if (bb2.first_[i] > bb1.second_[i] || bb2.second_[i] < bb1.first_[i])
            std::runtime_error("getIntersectionOfBoundingBoxes: no overlap!");
    /** #2 otherwise modify the first one to get the intersection. */
    for (int i = 0; i < dimension; ++i)
    {
        /** If the lower limit is inside change the lower limit. */
        if (bb1.first_[i] < bb2.first_[i] && bb2.first_[i] < bb1.second_[i])
            bb.first_[i] = bb2.first_[i];
        /**  If the upper limit is inside, change the upper limit. */
        if (bb1.second_[i] > bb2.second_[i] && bb2.second_[i] > bb1.first_[i])
            bb.second_[i] = bb2.second_[i];
    }
    return bb;
}

/** obtain minimum dimension of a bounding box */
template <class BoundingBoxType>
Real MinimumDimension(const BoundingBoxType &bbox)
{
    return (bbox.second_ - bbox.first_).cwiseAbs().minCoeff();
};
} // namespace SPH
#endif // SMALL_VECTORS_H
