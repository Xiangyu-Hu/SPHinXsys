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

#include "base_data_type.h"

namespace SPH
{

	SimTK::Vec2 EigenToSimTK(const Vec2d &eigen_vector);
	SimTK::Vec3 EigenToSimTK(const Vec3d &eigen_vector);
	Vec2d SimTKToEigen(const SimTK::Vec2 &simTK_vector);
	Vec3d SimTKToEigen(const SimTK::Vec3 &simTK_vector);

	SimTK::Mat22 EigenToSimTK(const Mat2d &eigen_matrix);
	SimTK::Mat33 EigenToSimTK(const Mat3d &eigen_matrix);
	Mat2d SimTKToEigen(const SimTK::Mat22 &simTK_matrix);
	Mat3d SimTKToEigen(const SimTK::Mat33 &simTK_matrix);

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

	/**
	 * @class Rotation2d
	 * @brief Rotation Coordinate transform (around the origin)
	 * in 2D with an angle.
	 */
	class Rotation2d
	{
		Real cosine_angle_, sine_angle_;

	public:
		explicit Rotation2d(Real angle)
			: cosine_angle_(std::cos(angle)), sine_angle_(std::sin(angle)){};
		virtual ~Rotation2d(){};

		/** Forward transformation. */
		Vec2d xformFrameVecToBase(const Vec2d &origin)
		{
			return Vec2d(origin[0] * cosine_angle_ - origin[1] * sine_angle_,
						 origin[1] * cosine_angle_ + origin[0] * sine_angle_);
		};
		/** Inverse transformation. */
		Vec2d xformBaseVecToFrame(const Vec2d &target)
		{
			return Vec2d(target[0] * cosine_angle_ + target[1] * sine_angle_,
						 target[1] * cosine_angle_ - target[0] * sine_angle_);
		};
	};

	/**
	 * @class Transform2d
	 * @brief Coordinate transform in 2D
	 * Note that the rotation is around the frame (or local) origin.
	 */
	class Transform2d
	{
	private:
		Rotation2d rotation_;
		Vec2d translation_;

	public:
		Transform2d() : rotation_(Rotation2d(0)), translation_(Vec2d::Zero()){};
		explicit Transform2d(const Vec2d &translation) : rotation_(Rotation2d(0)), translation_(translation){};
		explicit Transform2d(const Rotation2d &rotation, const Vec2d &translation = Vec2d::Zero())
			: rotation_(rotation), translation_(translation){};

		/** Forward rotation. */
		Vec2d xformFrameVecToBase(const Vec2d &origin)
		{
			return rotation_.xformFrameVecToBase(origin);
		};

		/** Forward transformation. Note that the rotation operation is carried out first. */
		Vec2d shiftFrameStationToBase(const Vec2d &origin)
		{
			return translation_ + xformFrameVecToBase(origin);
		};

		/** Inverse rotation. */
		Vec2d xformBaseVecToFrame(const Vec2d &target)
		{
			return rotation_.xformBaseVecToFrame(target);
		};

		/** Inverse transformation. Note that the inverse translation operation is carried out first. */
		Vec2d shiftBaseStationToFrame(const Vec2d &target)
		{
			return xformBaseVecToFrame(target - translation_);
		};
	};

	/**
	 * @class Transform3d
	 * @brief Wrapper for SimTK::Transform
	 * Note that the rotation is around the frame (or local) origin.
	 */
	class Transform3d
	{
	private:
		Mat3d rotation_;
		Vec3d translation_;
		SimTK::Transform transform_;

	public:
		Transform3d() : rotation_(Mat3d::Zero()), translation_(Vec3d::Zero())
		{
			// Default constructor gives an identity transform.
			transform_ = SimTK::Transform();
		};
		explicit Transform3d(const Vec3d &translation) : rotation_(Mat3d::Zero()), translation_(translation)
		{
			// Construct or default-convert a translation (expressed as a Vec3) into a transform with that translation and a zero rotation.
			transform_ = SimTK::Transform(EigenToSimTK(translation));
		};
		explicit Transform3d(const Mat3d &rotation, const Vec3d &translation = Vec3d::Zero())
			: rotation_(rotation), translation_(translation)
		{
			// Combine a rotation and a translation into a transform.
			transform_ = SimTK::Transform(SimTK::Rotation_<Real>(EigenToSimTK(rotation)), EigenToSimTK(translation));
		};

		/** Forward rotation. */
		Vec3d xformFrameVecToBase(const Vec3d &origin)
		{
			SimTK::Vec3 x_to_base = transform_.xformFrameVecToBase(EigenToSimTK(origin));
			return SimTKToEigen(x_to_base);
		};

		/** Forward transformation. Note that the rotation operation is carried out first. */
		Vec3d shiftFrameStationToBase(const Vec3d &origin)
		{
			SimTK::Vec3 frame_to_base = transform_.shiftFrameStationToBase(EigenToSimTK(origin));
			return SimTKToEigen(frame_to_base);
		};

		/** Inverse rotation. */
		Vec3d xformBaseVecToFrame(const Vec3d &target)
		{
			SimTK::Vec3 base_to_frame = transform_.xformBaseVecToFrame(EigenToSimTK(target));
			return SimTKToEigen(base_to_frame);
		};

		/** Inverse transformation. Note that the inverse translation operation is carried out first. */
		Vec3d shiftBaseStationToFrame(const Vec3d &target)
		{
			SimTK::Vec3 base_to_base = transform_.shiftBaseStationToFrame(EigenToSimTK(target));
			return SimTKToEigen(base_to_base);
		};
	};
}
#endif // SMALL_VECTORS_H
