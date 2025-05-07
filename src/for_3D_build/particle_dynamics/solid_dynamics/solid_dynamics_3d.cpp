#include "elastic_dynamics.h"

#include "polar_decomposition_3x3.h"

using namespace polar;

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
Vec3d UpdateElasticNormalDirection::getRotatedNormalDirection(const Mat3d &F, const Vec3d &n0)
{
    Real Q[9], H[9], A[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i * 3 + j] = F(i, j);

    polar::polar_decomposition(Q, H, A);
    // this decomposition has the form A = Q*H, where Q is orthogonal and H is symmetric positive semi-definite.
    // Ref. "An algorithm to compute the polar decomposition of a 3*3 matrix, Nicholas J. Higham et al. Numer Algor(2016) "
    Mat3d R;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R(i, j) = Q[i * 3 + j];
    return R * n0;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
