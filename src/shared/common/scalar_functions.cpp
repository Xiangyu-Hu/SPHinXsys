#include "scalar_functions.h"
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
int ThirdAxis(int first_axis)
{
    return SecondAxis(SecondAxis(first_axis));
}
//=================================================================================================//
Real getLeftStateInWeno(Real v_1, Real v_2, Real v_3, Real v_4)
{
    Real v1 = v_1;
    Real v2 = v_2;
    Real v3 = v_3;
    Real v4 = v_4;

    Real f1 = 0.5 * v2 + 0.5 * v3;
    Real f2 = (-0.5) * v1 + 1.5 * v2;
    Real f3 = v2 / 3.0 + 5.0 * v3 / 6.0 - v4 / 6.0;

    Real epsilon = 1.0e-6;
    Real s1 = pow(v2 - v3, 2) + epsilon;
    Real s2 = pow(v2 - v1, 2) + epsilon;
    Real s3 = pow(3.0 * v2 - 4.0 * v3 + v4, 2) / 4.0 + 13.0 * pow(v2 - 2.0 * v3 + v4, 2) / 12.0 + epsilon;
    Real s12 = 13.0 * pow(v1 - 2.0 * v2 + v3, 2) / 12.0 + pow(v1 - v3, 2) / 4.0 + epsilon;
    Real tau_4 = (v1 * (547.0 * v1 - 2522.0 * v2 + 1922.0 * v3 - 494.0 * v4) + v2 * (3443.0 * v2 - 5966.0 * v3 + 1602.0 * v4) + v3 * (2843.0 * v3 - 1642.0 * v4) + 267.0 * v4 * v4) / 240.0;

    Real alpha_1 = (1.0 + (tau_4 / s1) * (tau_4 / s12)) / 3.0;
    Real alpha_2 = (1.0 + (tau_4 / s2) * (tau_4 / s12)) / 6.0;
    Real alpha_3 = (1.0 + tau_4 / s3) / 2.0;
    Real w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
    Real w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
    Real w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);
    Real left_state = w_1 * f1 + w_2 * f2 + w_3 * f3;

    return left_state;
}
//=================================================================================================//
Real getRightStateInWeno(Real v_1, Real v_2, Real v_3, Real v_4)
{
    Real v1 = v_4;
    Real v2 = v_3;
    Real v3 = v_2;
    Real v4 = v_1;

    Real f1 = 0.5 * v2 + 0.5 * v3;
    Real f2 = (-0.5) * v1 + 1.5 * v2;
    Real f3 = v2 / 3.0 + 5.0 * v3 / 6.0 - v4 / 6.0;

    Real epsilon = 1.0e-6;
    Real s1 = pow(v2 - v3, 2) + epsilon;
    Real s2 = pow(v2 - v1, 2) + epsilon;
    Real s3 = pow(3.0 * v2 - 4.0 * v3 + v4, 2) / 4.0 + 13.0 * pow(v2 - 2.0 * v3 + v4, 2) / 12.0 + epsilon;
    Real s12 = 13.0 * pow(v1 - 2.0 * v2 + v3, 2) / 12.0 + pow(v1 - v3, 2) / 4.0 + epsilon;
    Real tau_4 = (v1 * (547.0 * v1 - 2522.0 * v2 + 1922.0 * v3 - 494.0 * v4) + v2 * (3443.0 * v2 - 5966.0 * v3 + 1602.0 * v4) + v3 * (2843.0 * v3 - 1642.0 * v4) + 267.0 * v4 * v4) / 240.0;

    Real alpha_1 = (1.0 + (tau_4 / s1) * (tau_4 / s12)) / 3.0;
    Real alpha_2 = (1.0 + (tau_4 / s2) * (tau_4 / s12)) / 6.0;
    Real alpha_3 = (1.0 + tau_4 / s3) / 2.0;
    Real w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
    Real w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
    Real w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);
    Real right_state = w_1 * f1 + w_2 * f2 + w_3 * f3;

    return right_state;
}
//=================================================================================================//
Real Heaviside(Real phi, Real half_width)
{
    Real normalized_phi = phi / half_width;
    return std::clamp(0.5 + 0.5 * normalized_phi, 0.0, 1.0);
}
//=================================================================================================//
} // namespace SPH
