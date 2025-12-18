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
    Eigen::JacobiSVD<Mat3d> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Mat3d U = svd.matrixU();
    Mat3d V = svd.matrixV();
    Mat3d R = U * V.transpose();
    return (R * n0).normalized();
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
