/**
 * @file 	small_vectors.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "small_vectors.h"
//=================================================================================================//
namespace SPH {

	//=================================================================================================//
	Vec2d FirstAxisVector(Vec2d zero_vector)
	{
		return Vec2d(1.0, 0.0);
	}
	//=================================================================================================//
	Vec3d FirstAxisVector(Vec3d zero_vector)
	{
		return Vec3d(1.0, 0.0, 0.0);
	};
	//=================================================================================================//
	Real getMinAbsoluteElement(Vec2d input)
	{
		Real min = input.norm();
		for (int n = 0; n != input.size(); n++) SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	Real getMinAbsoluteElement(Vec3d input)
	{
		Real min = input.norm();
		for (int n = 0; n != input.size(); n++) SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(Real input) {
		return Vec3d(input, 0.0, 0.0);
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(Vec2d input) {
		return Vec3d(input[0], input[1], 0.0);
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(Vec3d input) {
		return input;
	}
	//=================================================================================================//
	Mat2d getInverse(Mat2d& A)
	{
		Mat2d minv(0);
		SimTK::Real det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
		SimTK::Real invdet = 1.0 / det;
		minv(0, 0) = A(1, 1) * invdet;
		minv(0, 1) = -A(0, 1) * invdet;
		minv(1, 0) = -A(1, 0) * invdet;
		minv(1, 1) = A(0, 0) * invdet;
		return minv;
	}
	//=================================================================================================//
	Mat3d getInverse(Mat3d& A)
	{
		SimTK::Real det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
			A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
			A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

		SimTK::Real invdet = 1 / det;
		Mat3d minv(0); // inverse of matrix m
		minv(0, 0) = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) * invdet;
		minv(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) * invdet;
		minv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * invdet;
		minv(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) * invdet;
		minv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * invdet;
		minv(1, 2) = (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * invdet;
		minv(2, 0) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)) * invdet;
		minv(2, 1) = (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1)) * invdet;
		minv(2, 2) = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * invdet;

		return minv;
	}
	//=================================================================================================//
	Mat2d getAverageValue(Mat2d& A, Mat2d& B)
	{
		Mat2d C(1.0);
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				C(i, j) = 2.0 * A(i, j) * B(i, j) / (A(i, j) + B(i, j) + TinyReal);
			}
		}
		return C;
	}
	//=================================================================================================//
	Mat3d getAverageValue(Mat3d& A, Mat3d& B)
	{
		Mat3d C(1.0);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				C(i, j) = 2.0 * A(i, j) * B(i, j) / (A(i, j) + B(i, j) + TinyReal);
			}
		}
		return C;
	}
	//=================================================================================================//
	Mat2d inverseCholeskyDecomposition(Mat2d& A)
	{
		Mat2d lower(0);
		int n = 2;
		/** Decomposing a matrix into Lower Triangular. */
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < (i + 1); j++)
			{
				double sum = 0;
				for (int k = 0; k < j; k++)
				{
					sum += lower(i, k) * lower(j, k);
				}
				if (i == j)
				{
					lower(i, j) = sqrt(A(i, i) - sum);
				}
				else
				{
					lower(i, j) = (1.0 / lower(j, j) * (A(i, j) - sum));
				}
			}
		}
		Mat2d inverse_lower = getInverse(lower);
		return inverse_lower;
	}
	//=================================================================================================//
	Mat3d inverseCholeskyDecomposition(Mat3d& A)
	{
		Mat3d lower(0);
		int n = 3;
		/** Decomposing a matrix into Lower Triangular. */
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < (i + 1); j++)
			{
				double sum = 0;
				for (int k = 0; k < j; k++)
				{
					sum += lower(i, k) * lower(j, k);
				}
				if (i == j)
				{
					lower(i, j) = sqrt(A(i, i) - sum);
				}
				else
				{
					lower(i, j) = (1.0 / lower(j, j) * (A(i, j) - sum));
				}
			}
		}
		Mat3d inverse_lower = getInverse(lower);
		return inverse_lower;
	}
	//=================================================================================================//
}
