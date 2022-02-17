/**
 * @file 	vector_functions.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "vector_functions.h"
//=================================================================================================//
namespace SPH
{

	//=================================================================================================//
	Vec2d FirstAxisVector(const Vec2d &zero_vector)
	{
		return Vec2d(1.0, 0.0);
	}
	//=================================================================================================//
	Vec3d FirstAxisVector(const Vec3d &zero_vector)
	{
		return Vec3d(1.0, 0.0, 0.0);
	};
	//=================================================================================================//
	Real getMaxAbsoluteElement(const Vec2d &input)
	{
		Real max = 0.0;
		for (int n = 0; n != input.size(); n++)
			max = SMAX(fabs(input[n]), max);
		return max;
	}
	//=================================================================================================//
	Real getMaxAbsoluteElement(const Vec3d &input)
	{
		Real max = 0.0;
		for (int n = 0; n != input.size(); n++)
			max = SMAX(fabs(input[n]), max);
		return max;
	}
	//=================================================================================================//
	Real getMinAbsoluteElement(const Vec2d &input)
	{
		Real min = SimTK::Infinity;
		for (int n = 0; n != input.size(); n++)
			min = SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	Real getMinAbsoluteElement(const Vec3d &input)
	{
		Real min = SimTK::Infinity;
		for (int n = 0; n != input.size(); n++)
			min = SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(const Real &input)
	{
		return Vec3d(input, 0.0, 0.0);
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(const Vec2d &input)
	{
		return Vec3d(input[0], input[1], 0.0);
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(const Vec3d &input)
	{
		return input;
	}
	//=================================================================================================//
	Mat3d upgradeToMatrix3D(const Mat2d &input)
	{
		Mat3d output(0);
		output.col(0) = upgradeToVector3D(input.col(0));
		output.col(1) = upgradeToVector3D(input.col(1));
		return output;
	}
	//=================================================================================================//
	Mat3d upgradeToMatrix3D(const Mat3d &input)
	{
		return input;
	}
	//=================================================================================================//
	Mat2d getInverse(const Mat2d &A)
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
	Mat3d getInverse(const Mat3d &A)
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
	Mat2d getAverageValue(const Mat2d &A, const Mat2d &B)
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
	Mat3d getAverageValue(const Mat3d &A, const Mat3d &B)
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
	Mat2d inverseCholeskyDecomposition(const Mat2d &A)
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
	Mat3d inverseCholeskyDecomposition(const Mat3d &A)
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
	Mat2d getTransformationMatrix(const Vec2d &direction_of_y)
	{
		Mat2d transformation_matrix(0.0);
		transformation_matrix[0][0] = direction_of_y[1];
		transformation_matrix[0][1] = -direction_of_y[0];
		transformation_matrix[1][0] = direction_of_y[0];
		transformation_matrix[1][1] = direction_of_y[1];

		return transformation_matrix;
	}
	//=================================================================================================//
	Mat3d getTransformationMatrix(const Vec3d &direction_of_z)
	{
		Mat3d transformation_matrix(0.0);
		Real temp = 1.0 + direction_of_z[2];
		Real fraction = temp / (temp * temp + Eps);
		transformation_matrix[0][0] = direction_of_z[2] + powerN(direction_of_z[1], 2) * fraction;
		transformation_matrix[0][1] = -direction_of_z[0] * direction_of_z[1] * fraction;
		transformation_matrix[0][2] = -direction_of_z[0];
		transformation_matrix[1][0] = transformation_matrix[0][1];
		transformation_matrix[1][1] = direction_of_z[2] + powerN(direction_of_z[0], 2) * fraction;
		transformation_matrix[1][2] = -direction_of_z[1];
		transformation_matrix[2][0] = direction_of_z[0];
		transformation_matrix[2][1] = direction_of_z[1];
		transformation_matrix[2][2] = direction_of_z[2];

		return transformation_matrix;
	}
	//=================================================================================================//
	Mat2d getDiagonal(const Mat2d &A)
	{
		Mat2d diag(1.0);
		diag[0][0] = A[0][0];
		diag[1][1] = A[1][1];

		return diag;
	}
	Mat3d getDiagonal(const Mat3d &A)
	{
		Mat3d diag(1.0);
		diag[0][0] = A[0][0];
		diag[1][1] = A[1][1];
		diag[2][2] = A[2][2];

		return diag;
	}
	//=================================================================================================//
	Real getCosineOfAngleBetweenTwoVectors(const Vec2d &vector_1, const Vec2d &vector_2)
	{
		// returns the cosine of the angle between two vectors
		Real dot_product_1 = 0.0;
		for (int i = 0; i < vector_1.size(); i++)
		{
			dot_product_1 += vector_1[i] * vector_2[i];
		}
		Real cos_teta = dot_product_1 / (vector_1.norm() * vector_2.norm());

		return cos_teta;
	}
	//=================================================================================================//
	Real getCosineOfAngleBetweenTwoVectors(const Vec3d &vector_1, const Vec3d &vector_2)
	{
		// returns the cosine of the angle between two vectors
		Real dot_product_1 = 0.0;
		for (int i = 0; i < vector_1.size(); i++)
		{
			dot_product_1 += vector_1[i] * vector_2[i];
		}
		Real cos_teta = dot_product_1 / (vector_1.norm() * vector_2.norm());

		return cos_teta;
	}
	//=================================================================================================//
	Vec2d getVectorProjectionOfVector(const Vec2d &vector_1, const Vec2d &vector_2)
	{
		// get the projection of the vector_1 on vector 2, which is parallel to the vector_2, meaning it is the vector_2 * scalar
		Real dot_product_1 = 0.0;
		Real dot_product_2 = std::pow(vector_2.norm(), 2.0);
		for (int i = 0; i < vector_1.size(); i++)
		{
			dot_product_1 += vector_1[i] * vector_2[i];
		}
		// get scalar, which to multiply n_0 with
		Real lambda = dot_product_1 / dot_product_2;
		Vec2d proj_vector_1 = lambda * vector_2;

		return proj_vector_1;
	}
	//=================================================================================================//
	Vec3d getVectorProjectionOfVector(const Vec3d &vector_1, const Vec3d &vector_2)
	{
		// get the projection of the vector_1 on vector 2, which is parallel to the vector_2, meaning it is the vector_2 * scalar
		Real dot_product_1 = 0.0;
		Real dot_product_2 = std::pow(vector_2.norm(), 2.0);
		for (int i = 0; i < vector_1.size(); i++)
		{
			dot_product_1 += vector_1[i] * vector_2[i];
		}
		// get scalar, which to multiply n_0 with
		Real lambda = dot_product_1 / dot_product_2;
		Vec3d proj_vector_1 = lambda * vector_2;

		return proj_vector_1;
	}
	//=================================================================================================//
	bool checkIfPointInBoundingBox(const Vec2d &point, const BoundingBox2d &bbox)
	{
		return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
			   point[1] >= bbox.first[1] && point[1] <= bbox.second[1];
	}
	//=================================================================================================//
	bool checkIfPointInBoundingBox(const Vec3d &point, const BoundingBox3d &bbox)
	{
		return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
			   point[1] >= bbox.first[1] && point[1] <= bbox.second[1] &&
			   point[2] >= bbox.first[2] && point[2] <= bbox.second[2];
	}
	//=================================================================================================//
	Real MinimumDimension(const BoundingBox &bbox)
	{
		return getMinAbsoluteElement(bbox.second - bbox.first);
	}
	//=================================================================================================//
}
