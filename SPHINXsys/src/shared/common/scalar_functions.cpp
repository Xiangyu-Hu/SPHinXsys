/**
 * @file 	scalar_functions.cpp
 * @author	Xiangyu Hu
 * @version	0.1
 */

#include "scalar_functions.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	int ThirdAxis(int axis_direction) {
		return SecondAxis(SecondAxis(axis_direction));
	}
	//=================================================================================================//
	double getLeftStateInWeno(double v_1, double v_2, double v_3, double v_4)
	{
		double v1 = v_1;
		double v2 = v_2;
		double v3 = v_3;
		double v4 = v_4;

		double f1 = 0.5*v2 + 0.5*v3;
		double f2 = (-0.5)*v1 + 1.5*v2;
		double f3 = v2 / 3.0 + 5.0*v3 / 6.0 - v4 / 6.0;

		double epsilon = 1.0e-6;
		double s1 = pow(v2 - v3, 2) + epsilon;
		double s2 = pow(v2 - v1, 2) + epsilon;
		double s3 = pow(3.0*v2 - 4.0*v3 + v4, 2) / 4.0 + 13.0*pow(v2 - 2.0*v3 + v4, 2) / 12.0 + epsilon;
		double s12 = 13.0*pow(v1 - 2.0 * v2 + v3, 2) / 12.0 + pow(v1 - v3, 2) / 4.0 + epsilon;
		double tau_4 = (v1*(547.0*v1 - 2522.0*v2 + 1922.0*v3 - 494.0*v4)
			+ v2 * (3443.0*v2 - 5966.0*v3 + 1602.0*v4)
			+ v3 * (2843.0*v3 - 1642.0*v4)
			+ 267.0*v4*v4) / 240.0;

		double alpha_1 = (1.0 + (tau_4 / s1)*(tau_4 / s12)) / 3.0;
		double alpha_2 = (1.0 + (tau_4 / s2)*(tau_4 / s12)) / 6.0;
		double alpha_3 = (1.0 + tau_4 / s3) / 2.0;
		double w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
		double w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
		double w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);
		double left_state = w_1 * f1 + w_2 * f2 + w_3 * f3;

		return left_state;
	}
	//=================================================================================================//
	double getRightStateInWeno(double v_1, double v_2, double v_3, double v_4)
	{
		double v1 = v_4;
		double v2 = v_3;
		double v3 = v_2;
		double v4 = v_1;

		double f1 = 0.5*v2 + 0.5*v3;
		double f2 = (-0.5)*v1 + 1.5*v2;
		double f3 = v2 / 3.0 + 5.0*v3 / 6.0 - v4 / 6.0;

		double epsilon = 1.0e-6;
		double s1 = pow(v2 - v3, 2) + epsilon;
		double s2 = pow(v2 - v1, 2) + epsilon;
		double s3 = pow(3.0*v2 - 4.0*v3 + v4, 2) / 4.0 + 13.0*pow(v2 - 2.0*v3 + v4, 2) / 12.0 + epsilon;
		double s12 = 13.0*pow(v1 - 2.0 * v2 + v3, 2) / 12.0 + pow(v1 - v3, 2) / 4.0 + epsilon;
		double tau_4 = (v1*(547.0*v1 - 2522.0*v2 + 1922.0*v3 - 494.0*v4)
			+ v2 * (3443.0*v2 - 5966.0*v3 + 1602.0*v4)
			+ v3 * (2843.0*v3 - 1642.0*v4)
			+ 267.0*v4*v4) / 240.0;

		double alpha_1 = (1.0 + (tau_4 / s1)*(tau_4 / s12)) / 3.0;
		double alpha_2 = (1.0 + (tau_4 / s2)*(tau_4 / s12)) / 6.0;
		double alpha_3 = (1.0 + tau_4 / s3) / 2.0;
		double w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
		double w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
		double w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);
		double right_state = w_1 * f1 + w_2 * f2 + w_3 * f3;

		return right_state;
	}
	//=================================================================================================//

}
