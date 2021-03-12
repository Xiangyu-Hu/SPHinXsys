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
	Vec2d getVectorAfterRotation(Vec2d &initial_vector, Vec2d &rotation_angles)
	{
		/**The rotation matrix. */
		Mat2d rotation_matrix(0.0);
		rotation_matrix[0][0] = cos(rotation_angles[0]);
		rotation_matrix[0][1] = -sin(rotation_angles[0]);
		rotation_matrix[1][0] = -rotation_matrix[0][1];
		rotation_matrix[1][1] = rotation_matrix[0][0];

		return rotation_matrix * initial_vector;
	}
	//=================================================================================================//
	Vec3d getVectorAfterRotation(Vec3d &initial_vector, Vec3d &rotation_angles)
	{
		/**The rotation matrix about the X-axis. */
		Mat3d rotation_matrix_x(0.0);
		rotation_matrix_x[0][0] = 1.0;
		rotation_matrix_x[1][1] = cos(rotation_angles[0]);
		rotation_matrix_x[1][2] = -sin(rotation_angles[0]);
		rotation_matrix_x[2][1] = -rotation_matrix_x[1][2];
		rotation_matrix_x[2][2] = rotation_matrix_x[1][1];
		/**The rotation matrix about the Y-axis. */
		Mat3d rotation_matrix_y(0.0);
		rotation_matrix_y[0][0] = cos(rotation_angles[1]);
		rotation_matrix_y[0][2] = sin(rotation_angles[1]);
		rotation_matrix_y[1][1] = 1.0;
		rotation_matrix_y[2][0] = -rotation_matrix_y[0][2];
		rotation_matrix_y[2][2] = rotation_matrix_y[0][0];
		/**The rotation matrix about the Z-axis. */
		Mat3d rotation_matrix_z(0.0);
		rotation_matrix_z[0][0] = cos(rotation_angles[2]);
		rotation_matrix_z[0][1] = -sin(rotation_angles[2]);
		rotation_matrix_z[1][0] = -rotation_matrix_z[0][1];
		rotation_matrix_z[1][1] = rotation_matrix_z[0][0];
		rotation_matrix_z[2][2] = 1.0;

		return rotation_matrix_z * rotation_matrix_y * rotation_matrix_x * initial_vector;
	}
	//=================================================================================================//
	Vec2d getVectorChangeRateAfterRotation(Vec2d &initial_vector, Vec2d &rotation_angles, Vec2d &angular_vel)
	{
		/**The derivative of the rotation matrix. */
		Mat2d drotation_matrix_dt(0.0);
		drotation_matrix_dt[0][0] = -sin(rotation_angles[0]) * angular_vel[0];
		drotation_matrix_dt[0][1] = -cos(rotation_angles[0]) * angular_vel[0];
		drotation_matrix_dt[1][0] = -drotation_matrix_dt[0][1];
		drotation_matrix_dt[1][1] = drotation_matrix_dt[0][0];

		return drotation_matrix_dt * initial_vector;
	}
	//=================================================================================================//
	Vec3d getVectorChangeRateAfterRotation(Vec3d& initial_vector, Vec3d& rotation_angles, Vec3d& angular_vel)
	{
		/**The rotation matrix about the X-axis. */
		Mat3d rotation_matrix_x(0.0);
		rotation_matrix_x[0][0] = 1.0;
		rotation_matrix_x[1][1] = cos(rotation_angles[0]);
		rotation_matrix_x[1][2] = -sin(rotation_angles[0]);
		rotation_matrix_x[2][1] = -rotation_matrix_x[1][2];
		rotation_matrix_x[2][2] = rotation_matrix_x[1][1];
		/**The rotation matrix about the Y-axis. */
		Mat3d rotation_matrix_y(0.0);
		rotation_matrix_y[0][0] = cos(rotation_angles[1]);
		rotation_matrix_y[0][2] = sin(rotation_angles[1]);
		rotation_matrix_y[1][1] = 1.0;
		rotation_matrix_y[2][0] = -rotation_matrix_y[0][2];
		rotation_matrix_y[2][2] = rotation_matrix_y[0][0];
		/**The rotation matrix about the Z-axis. */
		Mat3d rotation_matrix_z(0.0);
		rotation_matrix_z[0][0] = cos(rotation_angles[2]);
		rotation_matrix_z[0][1] = -sin(rotation_angles[2]);
		rotation_matrix_z[1][0] = -rotation_matrix_z[0][1];
		rotation_matrix_z[1][1] = rotation_matrix_z[0][0];
		rotation_matrix_z[2][2] = 1.0;

		/**The derivative of the rotation matrix of the X-axis. */
		Mat3d drotation_matrix_x_dt(0.0);
		drotation_matrix_x_dt[1][1] = -sin(rotation_angles[0]) * angular_vel[0];
		drotation_matrix_x_dt[1][2] = -cos(rotation_angles[0]) * angular_vel[0];
		drotation_matrix_x_dt[2][1] = -drotation_matrix_x_dt[1][2];
		drotation_matrix_x_dt[2][2] = drotation_matrix_x_dt[1][1];
		/**The derivative of the rotation matrix of the Y-axis. */
		Mat3d drotation_matrix_y_dt(0.0);
		drotation_matrix_y_dt[0][0] = -sin(rotation_angles[1]) * angular_vel[1];
		drotation_matrix_y_dt[0][2] = cos(rotation_angles[1]) * angular_vel[1];
		drotation_matrix_y_dt[2][0] = -drotation_matrix_y_dt[0][2];
		drotation_matrix_y_dt[2][2] = drotation_matrix_y_dt[0][0];
		/**The derivative of the rotation matrix of the Z-axis. */
		Mat3d drotation_matrix_z_dt(0.0);
		drotation_matrix_z_dt[0][0] = -sin(rotation_angles[2]) * angular_vel[2];
		drotation_matrix_z_dt[0][1] = -cos(rotation_angles[2]) * angular_vel[2];
		drotation_matrix_z_dt[1][0] = -drotation_matrix_z_dt[0][1];
		drotation_matrix_z_dt[1][1] = drotation_matrix_z_dt[0][0];

		return (drotation_matrix_z_dt * rotation_matrix_y * rotation_matrix_x
			+ rotation_matrix_z * drotation_matrix_y_dt * rotation_matrix_x
			+ rotation_matrix_z * rotation_matrix_y * drotation_matrix_x_dt
			)* initial_vector;
	}
	//=================================================================================================//
}
