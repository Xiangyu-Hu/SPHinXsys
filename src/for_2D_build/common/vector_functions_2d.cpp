#include "vector_functions.h"

namespace SPH
{
//=================================================================================================//
Vec2d degradeToVecd(const Vec3d &input)
{
    Vec2d output = Vec2d::Zero();
    for (int i = 0; i != Dimensions; i++)
        output[i] = input[i];
    return output;
}
//=================================================================================================//
Mat2d degradeToMatd(const Mat3d &input)
{
    Mat2d output = Mat2d::Zero();
    for (int i = 0; i != Dimensions; i++)
        for (int j = 0; j != Dimensions; j++)
            output(i, j) = input(i, j);
    return output;
}
//=================================================================================================//
} // namespace SPH
