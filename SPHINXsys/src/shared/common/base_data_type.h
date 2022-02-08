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

	//vector with integers
	using Vec2i = SVec<2, int>;
	using Vec3i = SVec<3, int>;

	//vector with unsigned int
	using Vec2u = SVec<2, size_t>;
	using Vec3u = SVec<3, size_t>;

	//float point number
	using Real = SimTK::Real;

	//useful float point constants s
	const Real Pi = Real(M_PI);
	using SimTK::Eps;
	using SimTK::Infinity;
	using SimTK::TinyReal;
	constexpr size_t MaxSize_t = std::numeric_limits<size_t>::max();

	//vector with float point number
	using Vec2d = SimTK::Vec2;
	using Vec3d = SimTK::Vec3;

	//small matrix with float point number
	using Mat2d = SimTK::Mat22;
	using Mat3d = SimTK::Mat33;
	//small symmetric matrix with float point number
	using SymMat2d = SimTK::SymMat22;
	using SymMat3d = SimTK::SymMat33;

	//type trait for particle data type index
	template <typename T>
	struct ParticleDataTypeIndex
	{
		static constexpr int value = std::numeric_limits<int>::max();
	};
	template <>
	struct ParticleDataTypeIndex<Real>
	{
		static constexpr int value = 0;
	};
	template <>
	struct ParticleDataTypeIndex<Vec2d>
	{
		static constexpr int value = 1;
	};
	template <>
	struct ParticleDataTypeIndex<Vec3d>
	{
		static constexpr int value = 1;
	};
	template <>
	struct ParticleDataTypeIndex<Mat2d>
	{
		static constexpr int value = 2;
	};
	template <>
	struct ParticleDataTypeIndex<Mat3d>
	{
		static constexpr int value = 2;
	};
	template <>
	struct ParticleDataTypeIndex<int>
	{
		static constexpr int value = 3;
	};

	//verbal boolean for positive and negative axis directions
	const int xAxis = 0;
	const int yAxis = 1;
	const int zAxis = 2;
	const bool positiveDirection = true;
	const bool negativeDirection = false;

	/**
	 * @class Transform2d
	 * @brief Coordinate transfrom in 2D
	 */
	class Transform2d
	{
		Real rotation_angle_;
		Vec2d translation_;

	public:
		Transform2d(SimTK::Real rotation_angle)
			: rotation_angle_(rotation_angle), translation_(0){};
		Transform2d(SimTK::Real rotation_angle, Vec2d translation)
			: rotation_angle_(rotation_angle), translation_(translation){};
		/** Forward tranformation. */
		Vec2d imposeTransform(Vec2d &origin)
		{
			Vec2d target(origin[0] * cos(rotation_angle_) - origin[1] * sin(rotation_angle_),
						 origin[1] * cos(rotation_angle_) + origin[0] * sin(rotation_angle_));
			return target + translation_;
		};
		/** Inverse tranformation. */
		Vec2d imposeInverseTransform(Vec2d &target)
		{
			Vec2d origin(target[0] * cos(-rotation_angle_) - target[1] * sin(-rotation_angle_),
						 target[1] * cos(-rotation_angle_) + target[0] * sin(-rotation_angle_));
			return origin - translation_;
		};
	};

	/**
	 * @class Transform3d
	 * @brief Coordinate transfrom in 3D from SimTK
	 */
	using Transform3d = SimTK::Transform;
}

#endif //BASE_DATA_TYPE_H
