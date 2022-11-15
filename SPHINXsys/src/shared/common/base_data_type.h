/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
#ifndef BASE_DATA_TYPE_H
#define BASE_DATA_TYPE_H

#include "Simbody.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "scalar_functions.h"

#include <cmath>
#include <iostream>
#include <cassert>
#include <climits>
#include <algorithm>
#include <vector>
#include <map>

namespace SPH
{

	template <int N, class T>
	class SVec
	{
	private:
		T v[N];

	public:
		SVec<N, T>()
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)0;
		}

		explicit SVec<N, T>(T value_for_all)
		{
			for (int i = 0; i < N; ++i)
				v[i] = value_for_all;
		}

		template <class S>
		explicit SVec<N, T>(const S *source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		template <class S>
		explicit SVec<N, T>(const SVec<N, S> &source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		template <class S>
		explicit SVec<N, T>(const SimTK::Vec<N, S> &source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		SVec<N, T>(T v0, T v1)
		{
			v[0] = v0;
			v[1] = v1;
		}

		SVec<N, T>(T v0, T v1, T v2)
		{
			v[0] = v0;
			v[1] = v1;
			v[2] = v2;
		}

		T &operator[](int index)
		{
			assert(index >= 0 && index < N);
			return v[index];
		}

		const T &operator[](int index) const
		{
			assert(index >= 0 && index < N);
			return v[index];
		}

		bool nonzero(void) const
		{
			for (int i = 0; i < N; ++i)
				if (v[i])
					return true;
			return false;
		}

		bool operator==(const SVec<N, T> &b) const
		{
			bool res = true;
			for (int i = 0; i < N; ++i)
				res = res && (v[i] == b[i]);
			return res;
		}

		bool operator!=(const SVec<N, T> &b) const
		{
			bool res = false;
			for (int i = 0; i < N; ++i)
				res = res || (v[i] != b[i]);
			return res;
		}

		// Arithmetic operators
		SVec<N, T> operator=(const SVec<N, T> &b)
		{
			for (int i = 0; i < N; ++i)
				v[i] = b.v[i];
			return *this;
		}

		SVec<N, T> operator+(void) const
		{
			return SVec(*this);
		}

		SVec<N, T> operator+=(T a)
		{
			for (int i = 0; i < N; ++i)
				v[i] += a;
			return *this;
		}

		SVec<N, T> operator+(T a) const
		{
			SVec<N, T> w(*this);
			w += a;
			return w;
		}

		SVec<N, T> operator+=(const SVec<N, T> &w)
		{
			for (int i = 0; i < N; ++i)
				v[i] += w[i];
			return *this;
		}

		SVec<N, T> operator+(const SVec<N, T> &w) const
		{
			SVec<N, T> sum(*this);
			sum += w;
			return sum;
		}

		SVec<N, T> operator-=(T a)
		{
			for (int i = 0; i < N; ++i)
				v[i] -= a;
			return *this;
		}

		SVec<N, T> operator-(T a) const
		{
			SVec<N, T> w(*this);
			w -= a;
			return w;
		}

		SVec<N, T> operator-=(const SVec<N, T> &w)
		{
			for (int i = 0; i < N; ++i)
				v[i] -= w[i];
			return *this;
		}

		// unary minus
		SVec<N, T> operator-(void) const
		{
			SVec<N, T> negative;
			for (int i = 0; i < N; ++i)
				negative.v[i] = -v[i];
			return negative;
		}

		// minus
		SVec<N, T> operator-(const SVec<N, T> &w) const
		{
			SVec<N, T> diff(*this);
			diff -= w;
			return diff;
		}

		// scalar product
		SVec<N, T> operator*=(T a)
		{
			for (int i = 0; i < N; ++i)
				v[i] *= a;
			return *this;
		}

		SVec<N, T> operator*(T a) const
		{
			SVec<N, T> w(*this);
			w *= a;
			return w;
		}

		SVec<N, T> operator*=(const SVec<N, T> &w)
		{
			for (int i = 0; i < N; ++i)
				v[i] *= w.v[i];
			return *this;
		}

		SVec<N, T> operator*(const SVec<N, T> &w) const
		{
			SVec<N, T> componentwise_product;
			for (int i = 0; i < N; ++i)
				componentwise_product[i] = v[i] * w.v[i];
			return componentwise_product;
		}

		SVec<N, T> operator/=(T a)
		{
			for (int i = 0; i < N; ++i)
				v[i] /= a;
			return *this;
		}

		SVec<N, T> operator/(T a) const
		{
			SVec<N, T> w(*this);
			w /= a;
			return w;
		}

		SVec<N, T> operator/=(const SVec<N, T> &w)
		{
			for (int i = 0; i < N; ++i)
				v[i] /= w.v[i];
			return *this;
		}

		SVec<N, T> operator/(const SVec<N, T> &w) const
		{
			SVec<N, T> componentwise_divide;
			for (int i = 0; i < N; ++i)
				componentwise_divide[i] = v[i] / w.v[i];
			return componentwise_divide;
		}
	};

	template <int N, class T>
	std::ostream &operator<<(std::ostream &out, const SVec<N, T> &v)
	{
		out << '[' << v[0];
		for (int i = 1; i < N; ++i)
			out << ' ' << v[i];
		out << ']';
		return out;
	}

	template <int N, class T>
	std::istream &operator>>(std::istream &in, SVec<N, T> &v)
	{
		in >> v[0];
		for (int i = 1; i < N; ++i)
			in >> v[i];
		return in;
	}

	// vector with integers
	using Vec2i = SVec<2, int>;
	using Vec3i = SVec<3, int>;

	// vector with unsigned int
	using Vec2u = SVec<2, size_t>;
	using Vec3u = SVec<3, size_t>;

	// float point number
	using Real = SimTK::Real;

	// useful float point constants s
	const Real Pi = Real(M_PI);
	using SimTK::Eps;
	using SimTK::SqrtEps;
	using SimTK::Infinity;
	using SimTK::TinyReal;
	constexpr size_t MaxSize_t = std::numeric_limits<size_t>::max();
	constexpr double MinRealNumber = std::numeric_limits<double>::min();
	constexpr double MaxRealNumber = std::numeric_limits<double>::max();


	// vector with float point number
	using Vec2d = SimTK::Vec2;
	using Vec3d = SimTK::Vec3;

	// small matrix with float point number
	using Mat2d = SimTK::Mat22;
	using Mat3d = SimTK::Mat33;
	// small symmetric matrix with float point number
	using SymMat2d = SimTK::SymMat22;
	using SymMat3d = SimTK::SymMat33;

	// type trait for data type index
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

	// verbal boolean for positive and negative axis directions
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
		VecType first, second;
		int dimension_;

		BaseBoundingBox() : first(VecType(0)), second(VecType(0)), dimension_(VecType(0).size()){};
		BaseBoundingBox(const VecType &lower_bound, const VecType &upper_bound)
			: first(lower_bound), second(upper_bound),
			  dimension_(lower_bound.size()){};

		bool checkContain(const VecType &point)
		{
			bool is_contain = true;
			for (int i = 0; i < dimension_; ++i)
			{
				if (point[i] < first[i] || point[i] > second[i])
				{
					is_contain = false;
					break;
				}
			}
			return is_contain;
		};
	};

	template <class T>
	bool operator==(const BaseBoundingBox<T> &bb1, const BaseBoundingBox<T> &bb2)
	{
		return bb1.first == bb2.first && bb1.second == bb2.second ? true : false;
	}

	template <class BoundingBoxType>
	BoundingBoxType getIntersectionOfBoundingBoxes(const BoundingBoxType &bb1, const BoundingBoxType &bb2)
	{
		// check that the inputs are correct
		int dimension = bb1.dimension_;
		// Get the Bounding Box of the intersection of the two meshes
		BoundingBoxType bb(bb1);
		// #1 check that there is overlap, if not, exception
		for (int i = 0; i < dimension; ++i)
			if (bb2.first[i] > bb1.second[i] || bb2.second[i] < bb1.first[i])
				std::runtime_error("getIntersectionOfBoundingBoxes: no overlap!");
		// #2 otherwise modify the first one to get the intersection
		for (int i = 0; i < dimension; ++i)
		{
			// if the lower limit is inside change the lower limit
			if (bb1.first[i] < bb2.first[i] && bb2.first[i] < bb1.second[i])
				bb.first[i] = bb2.first[i];
			// if the upper limit is inside, change the upper limit
			if (bb1.second[i] > bb2.second[i] && bb2.second[i] > bb1.first[i])
				bb.second[i] = bb2.second[i];
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
		explicit Rotation2d(SimTK::Real angle)
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
		Transform2d() : rotation_(Rotation2d(0)), translation_(Vec2d(0)){};
		explicit Transform2d(const Vec2d &translation) : rotation_(Rotation2d(0)), translation_(translation){};
		explicit Transform2d(const Rotation2d &rotation, const Vec2d &translation = Vec2d(0))
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
	 * @brief Coordinate transform in 3D from SimTK
	 */
	using Transform3d = SimTK::Transform;
}

#endif // BASE_DATA_TYPE_H
