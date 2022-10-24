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
 *  HU1527/12-1 and Hu1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	base_data_type.h
 * @brief 	This is the date type definition for SPHinXsys. 
 * @author	Chi ZHang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
 */

#ifndef BASE_DATA_TYPE_H
#define BASE_DATA_TYPE_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <climits>
#include <algorithm>
#include <vector>
#include <map>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Simbody.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"

#include "scalar_functions.h"

namespace SPH
{
	/**
	 * Matrix<T, 2, 1>::Identity  return {1,0}
	 * Matrix<T, 2, 2>::Identity  return {{1,0},
	 *								      {0,1},}
	 * Matrix<T, n, n>::Ones  Set all element to One. 						
	 */
	using Real = double;
	/** Vector with integers. */
	using Vec2i = Eigen::Matrix<int, 2, 1>;
	using Vec3i = Eigen::Matrix<int, 3, 1>;
	/** Vector with unsigned int. */
	using Vec2u = Eigen::Matrix<size_t, 2, 1>;
	using Vec3u = Eigen::Matrix<size_t, 3, 1>;
	/** Useful float point constants. */
	constexpr size_t MaxSize_t = std::numeric_limits<size_t>::max();
	constexpr double MinRealNumber = std::numeric_limits<double>::min();
	constexpr double MaxRealNumber = std::numeric_limits<double>::max();
	/** Vector with float point number.*/
	using Vec2d = Eigen::Matrix<Real, 2, 1>;
	using Vec3d = Eigen::Matrix<Real, 3, 1>;
	/** Small, 2*2 and 3*3, matrix with float point number. */
	using Mat2d = Eigen::Matrix<Real, 2, 2>;
	using Mat3d = Eigen::Matrix<Real, 3, 3>;
	/** Unified initialize to zero for all data type. */
	template <typename DataType> 
	struct DataTypeInitializer
	{
		static inline DataType zero;
	};
	template<> struct DataTypeInitializer<int>   {static inline int   zero = 0;};
	template<> struct DataTypeInitializer<Real>  {static inline Real  zero = 0.0;};
	template<> struct DataTypeInitializer<Vec2d> {static inline Vec2d zero = Vec2d::Zero();};
	template<> struct DataTypeInitializer<Vec3d> {static inline Vec3d zero = Vec3d::Zero();};
	template<> struct DataTypeInitializer<Mat2d> {static inline Mat2d zero = Mat2d::Zero();};
	template<> struct DataTypeInitializer<Mat3d> {static inline Mat3d zero = Mat3d::Zero();};
	/** Type trait for data type index. */
	template <typename T>
	struct DataTypeIndex
	{
		static constexpr int value = std::numeric_limits<int>::max();
	};
	template <>
	struct DataTypeIndex<Real>
	{
		static constexpr int value = 0;
	};
	template <>
	struct DataTypeIndex<Vec2d>
	{
		static constexpr int value = 1;
	};
	template <>
	struct DataTypeIndex<Vec3d>
	{
		static constexpr int value = 1;
	};
	template <>
	struct DataTypeIndex<Mat2d>
	{
		static constexpr int value = 2;
	};
	template <>
	struct DataTypeIndex<Mat3d>
	{
		static constexpr int value = 2;
	};
	template <>
	struct DataTypeIndex<int>
	{
		static constexpr int value = 3;
	};
	/** Verbal boolean for positive and negative axis directions. */
	const int xAxis = 0;
	const int yAxis = 1;
	const int zAxis = 2;
	const bool positiveDirection = true;
	const bool negativeDirection = false;
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
		/**
	 * @class Rotation2d
	 * @brief Rotation Coordinate transform (around the origin)
	 * in 2D with an angle.
	 */
	class Rotation2d
	{
		Real angle_, cosine_angle_, sine_angle_;

	public:
		explicit Rotation2d(Real angle)
			: angle_(angle), cosine_angle_(std::cos(angle)), sine_angle_(std::sin(angle)){};
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
	
	using Transform3d = SimTK::Transform;

	/**
	 * @class Transform3d
	 * @brief Using the SimTK::Transform from simbody library. 
	 */
	// class Transform3d
	// {
	// private:
	// 	Vec3d translation_;
	// 	SimTK::Transform simtk_transform_;

	// public:
	// 	Transform3d() : rotation_(Rotation3d(0)), translation_(Vec3d(0)){};
	// 	explicit Transform3d(const Vec3d &translation) : rotation_(Rotation3d(0)), translation_(translation){};
	// 	explicit Transform3d(const Rotation3d &rotation, const Vec3d &translation = Vec2d(0))
	// 		: rotation_(rotation), translation_(translation){};

	// 	/** Forward rotation. */
	// 	Vec3d xformFrameVecToBase(const Vec3d &origin)
	// 	{
	// 		return rotation_.xformFrameVecToBase(origin);
	// 	};

	// 	/** Forward transformation. Note that the rotation operation is carried out first. */
	// 	Vec3d shiftFrameStationToBase(const Vec3d &origin)
	// 	{
	// 		return translation_ + xformFrameVecToBase(origin);
	// 	};

	// 	/** Inverse rotation. */
	// 	Vec3d xformBaseVecToFrame(const Vec3d &target)
	// 	{
	// 		return rotation_.xformBaseVecToFrame(target);
	// 	};

	// 	/** Inverse transformation. Note that the inverse translation operation is carried out first. */
	// 	Vec3d shiftBaseStationToFrame(const Vec3d &target)
	// 	{
	// 		return xformBaseVecToFrame(target - translation_);
	// 	};

	/**
	 * @brief Convert any input to string and pad the output with zeros
	 * @todo Use external library for general string formatting, e.g. abseil, fmt library, or std::format
	 */
	template <typename T>
	std::string padValueWithZeros(T &&value, size_t max_string_width = 10)
	{
		std::ostringstream s_time;
		s_time << std::setw(max_string_width) << std::setfill('0') << value;
		return s_time.str();
	}
	
	/** Constant parameters. */
	const Real Pi = Real(M_PI);
	const Real Eps = 2.22045e-16;
	const Real TinyReal = 2.71051e-20;
	const Real Infinity = std::numeric_limits<double>::max();
}

#endif // BASE_DATA_TYPE_H
