#ifndef SPHINXSYS_BASE_SMALLVEC_H
#define SPHINXSYS_BASE_SMALLVEC_H

#pragma once

#include <cmath>
#include <iostream>
#include <cassert>

#include "Simbody.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "scalar_functions.h"

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

	//vetor with integers
	using Vec1i = SVec<1, int>;
	using Vec2i = SVec<2, int>;
	using Vec3i = SVec<3, int>;

	//vector with unsigned int
	using Vec1u = SVec<1, size_t>;
	using Vec2u = SVec<2, size_t>;
	using Vec3u = SVec<3, size_t>;

	//vector with double float number
	using Vec1d = SimTK::Vec1;
	using Vec2d = SimTK::Vec2;
	using Vec3d = SimTK::Vec3;

	//small matrix with double float number
	using Mat1d = SimTK::Mat11;
	using Mat2d = SimTK::Mat22;
	using Mat3d = SimTK::Mat33;
	//small symmetric matrix with double float number
	using SymMat1d = SimTK::SymMat11;
	using SymMat2d = SimTK::SymMat22;
	using SymMat3d = SimTK::SymMat33;

	const SimTK::Real pi = SimTK::Pi;

	template<int N>
	SimTK::Vec<N> normalize(const SimTK::Vec<N> &r) { return r / (r.norm() + 1.0e-15); };

	Mat2d GeneralizedInverse(Mat2d &A);
	Mat3d GeneralizedInverse(Mat3d &A);
	Mat2d getInverse(Mat2d &A);
	Mat3d getInverse(Mat3d &A);
	Mat2d getAverageValue(Mat2d &A, Mat2d &B);
	Mat3d getAverageValue(Mat3d &A, Mat3d &B);
	Mat2d inverseCholeskyDecomposition(Mat2d &A);
	Mat3d inverseCholeskyDecomposition(Mat3d &A);

	Vec2d FisrtAxisVector(Vec2d zero_vector);
	Vec3d FisrtAxisVector(Vec3d zero_vector);

	SimTK::Real getMinAbslouteElement(Vec2d input);
	SimTK::Real getMinAbslouteElement(Vec3d input);

	SimTK::Real TensorDoubleDotProduct(Mat2d &A, Mat2d &B);
	SimTK::Real TensorDoubleDotProduct(Mat3d &A, Mat3d &B);

	/** Upgrade real to 3d vector. */
	Vec3d upgradeToVector3D(SimTK::Real input); 
	/** Upgrade 2d vector to 3d vector. */
	Vec3d upgradeToVector3D(Vec2d input);
	/** overload for 3d vector. */
	Vec3d upgradeToVector3D(Vec3d input);
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
		Vec2d ImposeTransform(Vec2d& origin) {
			Vec2d result(origin[0] * cos(rotation_angle_) - origin[1] * sin(rotation_angle_),
				origin[1] * cos(rotation_angle_) + origin[0] * sin(rotation_angle_));
				return result + translation_;
		};
		/** Inverse tranformation. */
		Vec2d ImposeInverseTransform(Vec2d& result) {
			Vec2d origin(result[0] * cos(-rotation_angle_) - result[1] * sin(-rotation_angle_),
				result[1] * cos(-rotation_angle_) + result[0] * sin(-rotation_angle_));
			return origin - translation_;
		};
	};
	/** CrossProduct computation. */
	Vec2d getCrossProduct(Vec2d &A, Vec2d &B);
	Vec3d getCrossProduct(Vec3d &A, Vec3d &B);
	/** User defined cross product for fiber calculation. */
	Mat2d getCrossProductMatrix(Vec2d &A);
	Mat3d getCrossProductMatrix(Vec3d &A);
}

#endif //SPHINXSYS_BASE_SMALLVEC_H
