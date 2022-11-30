/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	vector_functions.h
 * @brief 	Basic functions for vector type data.
 * @author	Chi ZHang and Xiangyu Hu
 */
#ifndef SMALL_VECTORS_H
#define SMALL_VECTORS_H

#include "data_type.h"

namespace SPH
{

	SimTK::Vec2 EigenToSimTK(const Vec2d &eigen_vector);
	SimTK::Vec3 EigenToSimTK(const Vec3d &eigen_vector);
	Vec2d SimTKToEigen(const SimTK::Vec2 &simTK_vector);
	Vec3d SimTKToEigen(const SimTK::Vec3 &simTK_vector);

	Vec2d FirstAxisVector(const Vec2d &zero_vector);
	Vec3d FirstAxisVector(const Vec3d &zero_vector);

	Vec3d upgradeToVector3D(const Real &input);
	Vec3d upgradeToVector3D(const Vec2d &input);
	Vec3d upgradeToVector3D(const Vec3d &input);
	Mat3d upgradeToMatrix3D(const Mat2d &input);
	Mat3d upgradeToMatrix3D(const Mat3d &input);

	Mat2d getInverse(const Mat2d &A);
	Mat3d getInverse(const Mat3d &A);
	Mat2d getAverageValue(const Mat2d &A, const Mat2d &B);
	Mat3d getAverageValue(const Mat3d &A, const Mat3d &B);
	Mat2d inverseCholeskyDecomposition(const Mat2d &A);
	Mat3d inverseCholeskyDecomposition(const Mat3d &A);
	Mat2d getDiagonal(const Mat2d &A);
	Mat3d getDiagonal(const Mat3d &A);

	/** double dot product between two matrices, resulting in a scalar value (sum of products of element-wise) */
	Real CalculateDoubleDotProduct(Mat2d Matrix1, Mat2d Matrix2); // calculate double dot
	Real CalculateDoubleDotProduct(Mat3d Matrix1, Mat3d Matrix2); // calculate double dot

	/** get transformation matrix. */
	Mat2d getTransformationMatrix(const Vec2d &direction_of_y);
	Mat3d getTransformationMatrix(const Vec3d &direction_of_z);

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

	/** obtain minimum dimension of a bounding box */
	Real MinimumDimension(const BoundingBox &bbox);
}
#endif // SMALL_VECTORS_H
