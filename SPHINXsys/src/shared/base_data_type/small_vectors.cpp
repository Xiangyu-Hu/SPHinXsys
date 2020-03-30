/**
 * @file 	small_vectors.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "small_vectors.h"
//=================================================================================================//
namespace SPH {

	//=================================================================================================//
	Vec2d FisrtAxisVector(Vec2d zero_vector) 
	{ 
		return Vec2d(1.0, 0.0); 
	}
	//=================================================================================================//
	Vec3d FisrtAxisVector(Vec3d zero_vector) 
	{ 
		return Vec3d(1.0, 0.0, 0.0); 
	};
	//=================================================================================================//
	SimTK::Real getMinAbslouteElement(Vec2d input)
	{
		SimTK::Real min = input.norm();
		for (size_t n = 0; n != input.size(); n++) SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	SimTK::Real getMinAbslouteElement(Vec3d input)
	{
		SimTK::Real min = input.norm();
		for (size_t n = 0; n != input.size(); n++) SMIN(fabs(input[n]), min);
		return min;
	}
	//=================================================================================================//
	Vec3d upgradeToVector3D(SimTK::Real input) {
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
	/** 2x2 matrix,
	  * It is only works for kernel corection not reverse other matrices */
	Mat2d GeneralizedInverse(Mat2d &A)
	{
		// Find U such that U*A*A�*U� = diag
		Mat2d Su = A * ~A;
		SimTK::Real phi = 0.5*atan2(Su(0, 1) + Su(1, 0), Su(0, 0) - Su(1, 1));
		SimTK::Real Cphi = cos(phi);
		SimTK::Real Sphi = sin(phi);
		Mat2d U(Cphi, -Sphi, Sphi, Cphi);

		// Find W such that W�*A�*A*W = diag
		Mat2d Sw = ~A*A;
		SimTK::Real theta = 0.5*atan2(Sw(0, 1) + Sw(1, 0), Sw(0, 0) - Sw(1, 1));
		SimTK::Real Ctheta = cos(theta);
		SimTK::Real Stheta = sin(theta);
		Mat2d W(Ctheta, -Stheta, Stheta, Ctheta);

		// Find the singular values from U
		SimTK::Real SUsum = Su(0, 0) + Su(1, 1);
		SimTK::Real SUdif = sqrt((Su(0, 0) - Su(1, 1))*(Su(0, 0) - Su(1, 1)) + 4.0 * Su(0, 1)*Su(1, 0));
		Vec2d svals(sqrt((SUsum + SUdif) / 2), sqrt((SUsum - SUdif) / 2));
		Mat2d SIG(svals(0), 0.0, 0.0, svals(1));

		// Find the correction matrix for the right side
		Mat2d S = ~U*A*W;
		Mat2d C(SimTK::sign(S(0, 0)), 0.0, 0.0, SimTK::sign(S(1, 1)));
		Mat2d V = W * C;

		//regularization
		SimTK::Real strength = 0.1;
		Mat2d r_SIG(svals(0) / (svals(0)*svals(0) + strength *(svals(0) - 1.0)*(svals(0) - 1.0)),
			0.0, 0.0,
			svals(1) / (svals(1)*svals(1) + strength *(svals(1) - 1.0)*(svals(1) - 1.0)));

		return V * r_SIG*~U;
	}
	/** 3x3 matrix,
	  * It is only works for kernel corection not reverse other matrices */
//=================================================================================================//
	Mat3d GeneralizedInverse(Mat3d &A)
	{
		SimTK::Real det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
             A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
             A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

        SimTK::Real invdet = 1 / det;
        Mat3d minv; // inverse of matrix m
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
	SimTK::Real TensorDoubleDotProduct(Mat2d &A, Mat2d &B)
	{
		return  (A*~B).trace();
	}
//=================================================================================================//
	SimTK::Real TensorDoubleDotProduct(Mat3d &A, Mat3d &B)
	{
		return  (A*~B).trace();
	}
//=================================================================================================//
	Mat2d getInverse(Mat2d &A)
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
	Mat3d getInverse(Mat3d &A)
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
	Mat2d getAverageValue(Mat2d &A, Mat2d &B)
	{
		Mat2d C(1.0);
		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < 2; j++)
			{
				C(i,j) = 2.0 * A(i,j) * B(i,j) / (A(i,j) + B(i,j) + 1.0e-15);
			}
		}
		return C;
	}
//=================================================================================================//
	Mat3d getAverageValue(Mat3d &A, Mat3d &B)
	{
		Mat3d C(1.0);
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				C(i,j) = 2.0 * A(i,j) * B(i,j) / (A(i,j) + B(i,j) + 1.0e-15);
			}
		}
		return C;
	}
//=================================================================================================//
	Mat2d inverseCholeskyDecomposition(Mat2d &A)
	{
		Mat2d lower(0);
		int n = 2;
    	/** Decomposing a matrix into Lower Triangular. */
    	for (int i = 0; i < n; i++)
    	{
        	for (int j = 0; j < (i+1); j++) 
        	{
            	double sum = 0;
            	for (int k = 0; k < j; k++)
            	{
                	sum += lower(i,k) * lower(j,k);
            	}
            	if(i == j)
            	{
                	lower(i,j) = sqrt(A(i,i) - sum);
            	}else
            	{
                	lower(i,j) = (1.0 / lower(j,j) * (A(i,j) - sum));
            	}
        	}
    	}
		Mat2d inverse_lower = getInverse(lower);
		return inverse_lower;
	}
//=================================================================================================//
	Mat3d inverseCholeskyDecomposition(Mat3d &A)
	{
		Mat3d lower(0);
		int n = 3;
    	/** Decomposing a matrix into Lower Triangular. */
    	for (int i = 0; i < n; i++)
    	{
        	for (int j = 0; j < (i+1); j++) 
        	{
            	double sum = 0;
            	for (int k = 0; k < j; k++)
            	{
                	sum += lower(i,k) * lower(j,k);
            	}
            	if(i == j)
            	{
                	lower(i,j) = sqrt(A(i,i) - sum);
            	}else
            	{
                	lower(i,j) = (1.0 / lower(j,j) * (A(i,j) - sum));
            	}
        	}
    	}
		Mat3d inverse_lower = getInverse(lower);
		return inverse_lower;
	}
//=================================================================================================//
	Vec2d getCrossProduct(Vec2d &A, Vec2d &B)
	{
		return Vec2d(0.0);
	}
//=================================================================================================//
	Vec3d getCrossProduct(Vec3d &A, Vec3d &B)
	{
		SimTK::Real x_1 = A[1] * B[2] - A[2] * B[1];
		SimTK::Real x_2 = A[2] * B[0] - A[0] * B[2];
		SimTK::Real x_3 = A[0] * B[1] - A[1] * B[0];
		return Vec3d(x_1, x_2, x_3);
	}
//=================================================================================================//
	Mat2d getCrossProductMatrix(Vec2d &A)
	{
		return Mat2d(0.0);
	}
//=================================================================================================//
	Mat3d getCrossProductMatrix(Vec3d &A)
	{
		Mat3d cross_A(0.0);
		cross_A(1,1) =  0.0;
		cross_A(1,2) = -A[2];
		cross_A(1,3) =  A[1];

		cross_A(2,1) =  A[2];
		cross_A(2,2) =  0.0;
		cross_A(2,3) = -A[0];

		cross_A(3,1) = -A[1];
		cross_A(3,2) =  A[0];
		cross_A(3,3) =  0.0;
		
		return cross_A;
	}
//=================================================================================================//
}
