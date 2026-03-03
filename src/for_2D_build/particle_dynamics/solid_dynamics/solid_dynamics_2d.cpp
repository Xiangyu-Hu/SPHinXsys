#include "elastic_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
//=========================================================================================//
Vec2d UpdateElasticNormalDirection::getRotatedNormalDirection(const Mat2d &F, const Vec2d &n0)
{
    // deformation tensor in 2D
    Real F00 = F(0, 0);
    Real F01 = F(0, 1);
    Real F10 = F(1, 0);
    Real F11 = F(1, 1);
    // polar decomposition
    Mat2d R{
        {F00 + F11, F01 - F10},
        {F10 - F01, F00 + F11},
    };

    return R * n0 / sqrt(pow(F00 + F11, 2) + pow(F01 - F10, 2));
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
