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
#ifndef SPHINXSYS_BASE_SMALLVEC_H
#define SPHINXSYS_BASE_SMALLVEC_H



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

using namespace std;
namespace SPH {

	template<int N, class T>
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

		template<class S>
		explicit SVec<N, T>(const S *source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		template <class S>
		explicit SVec<N, T>(const SVec<N, S>& source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		template <class S>
		explicit SVec<N, T>(const SimTK::Vec<N, S>& source)
		{
			for (int i = 0; i < N; ++i)
				v[i] = (T)source[i];
		}

		SVec<N, T>(T v0, T v1)
		{
			v[0] = v0; v[1] = v1;
		}

		SVec<N, T>(T v0, T v1, T v2)
		{
			v[0] = v0; v[1] = v1; v[2] = v2;
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
				if (v[i]) return true;
			return false;
		}

		bool operator==(const SVec<N, T>& b) const
		{
			bool res = true;
			for (int i = 0; i < N; ++i)
				res = res && (v[i] == b[i]);
			return res;
		}

		bool operator!=(const SVec<N, T>& b) const
		{
			bool res = false;
			for (int i = 0; i < N; ++i)
				res = res || (v[i] != b[i]);
			return res;
		}

		// Arithmetic operators
		SVec<N, T> operator=(const SVec<N, T>& b)
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

	template<int N, class T>
	std::ostream &operator<<(std::ostream &out, const SVec<N, T> &v)
	{
		out << '[' << v[0];
		for (int i = 1; i < N; ++i)
			out << ' ' << v[i];
		out << ']';
		return out;
	}

	template<int N, class T>
	std::istream &operator >> (std::istream &in, SVec<N, T> &v)
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

	//useful float point constants 
	const SimTK::Real Pi = SimTK::Pi;
	const SimTK::Real Infinity = numeric_limits<SimTK::Real>::max();
	const SimTK::Real Eps = SimTK::Eps;
	const SimTK::Real TinyReal = SimTK::TinyReal;
	constexpr size_t MaxSize_t = numeric_limits<size_t>::max();

	//vector with float point number
	using Vec2d = SimTK::Vec2;
	using Vec3d = SimTK::Vec3;

	//small matrix with float point number
	using Mat2d = SimTK::Mat22;
	using Mat3d = SimTK::Mat33;
	//small symmetric matrix with float point number
	using SymMat2d = SimTK::SymMat22;
	using SymMat3d = SimTK::SymMat33;

	//particle data type index
	const int indexScalar = 0;
	const int indexVector = 1;
	const int indexMatrix = 2;
	const int indexInteger = 3;
	const int indexBoolean = 4;

	Vec2d FirstAxisVector(Vec2d zero_vector);
	Vec3d FirstAxisVector(Vec3d zero_vector);
	Real getMinAbsoluteElement(Vec2d input);
	Real getMinAbsoluteElement(Vec3d input);
	Vec3d upgradeToVector3D(Real input);
	Vec3d upgradeToVector3D(Vec2d input);
	Vec3d upgradeToVector3D(Vec3d input);

	template<typename OutVectorType>
	OutVectorType upgradeVector(Real input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input;
		return out_vector;
	};
	template<typename OutVectorType>
	OutVectorType upgradeVector(Vec2d input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input[0];
		out_vector[1] = input[1];
		return out_vector;
	};
	template<typename OutVectorType>
	OutVectorType upgradeVector(Vec3d input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input[0];
		out_vector[1] = input[1];
		out_vector[2] = input[2];
		return out_vector;
	};

	Mat2d getInverse(Mat2d &A);
	Mat3d getInverse(Mat3d &A);
	Mat2d getAverageValue(Mat2d &A, Mat2d &B);
	Mat3d getAverageValue(Mat3d &A, Mat3d &B);
	Mat2d inverseCholeskyDecomposition(Mat2d &A);
	Mat3d inverseCholeskyDecomposition(Mat3d &A);

	/**
	 * @class Transform2d
	 * @brief Coordinate transfrom in 2D
	 */
	class Transform2d 
	{
		using Real = SimTK::Real;
		Real rotation_angle_;
		Vec2d translation_;
	public:
		Transform2d(SimTK::Real rotation_angle)
			: rotation_angle_(rotation_angle), translation_(0) {};
		Transform2d(SimTK::Real rotation_angle, Vec2d translation)
			: rotation_angle_(rotation_angle), translation_(translation) {};
		/** Forward tranformation. */
		Vec2d imposeTransform(Vec2d& origin) {
			Vec2d result(origin[0] * cos(rotation_angle_) - origin[1] * sin(rotation_angle_),
				origin[1] * cos(rotation_angle_) + origin[0] * sin(rotation_angle_));
				return result + translation_;
		};
		/** Inverse tranformation. */
		Vec2d imposeInverseTransform(Vec2d& result) {
			Vec2d origin(result[0] * cos(-rotation_angle_) - result[1] * sin(-rotation_angle_),
				result[1] * cos(-rotation_angle_) + result[0] * sin(-rotation_angle_));
			return origin - translation_;
		};
	};

	/** 
	* @function getVectorAfterRotation
	* @brief Each of these basic vector rotations appears counterclockwise 
	* @brief when the axis about which they occur points toward the observer, 
	* @brief and the coordinate system is right-handed. 
	*/
	Vec2d getVectorAfterRotation(Vec2d &initial_vector, Vec2d &rotation_angles);
	Vec3d getVectorAfterRotation(Vec3d &initial_vector, Vec3d &rotation_angles);

	/** Vector change rate after rotation. */
	Vec2d getVectorChangeRateAfterRotation(Vec2d &initial_vector, Vec2d &rotation_angles, Vec2d &angular_vel);
	Vec3d getVectorChangeRateAfterRotation(Vec3d &initial_vector, Vec3d &rotation_angles, Vec3d &angular_vel);
}

#endif //SPHINXSYS_BASE_SMALLVEC_H
