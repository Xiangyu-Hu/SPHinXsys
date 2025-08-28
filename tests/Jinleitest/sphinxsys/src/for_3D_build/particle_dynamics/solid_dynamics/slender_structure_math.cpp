#include "slender_structure_math.h"

namespace SPH
{
//=====================================================================================================//
namespace slender_structure_dynamics
{
//=================================================================================================//

Vec3d getVectorAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles)
{
    Real theta = sqrt(rotation_angles[0] * rotation_angles[0] +
                      rotation_angles[1] * rotation_angles[1] + rotation_angles[2] * rotation_angles[2]);
    Mat3d R1 = Mat3d::Zero();
    R1(0, 1) = -rotation_angles[2];
    R1(0, 2) = rotation_angles[1];
    R1(1, 0) = rotation_angles[2];
    R1(1, 2) = -rotation_angles[0];
    R1(2, 0) = -rotation_angles[1];
    R1(2, 1) = rotation_angles[0];

    Mat3d rotation_matrix = Mat3d::Identity() + sin(theta) / (theta + Eps) * R1 +
                            (1 - cos(theta)) / (theta * theta + Eps) * R1 * R1;

    return rotation_matrix * initial_vector;
}

//=================================================================================================//
Vec3d getVectorChangeRateAfterThinStructureRotation(
    const Vec3d &initial_vector, const Vec3d &rotation_angles, const Vec3d &angular_vel)
{
    Real sin_rotation_0 = sin(rotation_angles[0]);
    Real cos_rotation_0 = cos(rotation_angles[0]);

    Real sin_rotation_1 = sin(rotation_angles[1]);
    Real cos_rotation_1 = cos(rotation_angles[1]);

    Real dpseudo_n_dt_0 = -sin_rotation_0 * sin_rotation_1 * angular_vel[0] +
                          cos_rotation_0 * cos_rotation_1 * angular_vel[1];
    Real dpseudo_n_dt_1 = -cos_rotation_0 * angular_vel[0];
    Real dpseudo_n_dt_2 = -sin_rotation_0 * cos_rotation_1 * angular_vel[0] -
                          cos_rotation_0 * sin_rotation_1 * angular_vel[1];

    return Vec3d(dpseudo_n_dt_0, dpseudo_n_dt_1, dpseudo_n_dt_2);
}
//=================================================================================================//

Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
{
    Real sin_rotation_0 = sin(rotation[0]);
    Real cos_rotation_0 = cos(rotation[0]);
    Real sin_rotation_1 = sin(rotation[1]);
    Real cos_rotation_1 = cos(rotation[1]);

    Real rotation_0_a = -(dpseudo_n_d2t[2] * cos_rotation_1 + dpseudo_n_d2t[0] * sin_rotation_1 +
                          angular_vel[1] * angular_vel[1] * cos_rotation_0 + angular_vel[0] * angular_vel[0] * cos_rotation_0);
    Real rotation_0_b = sin_rotation_0 * angular_vel[0] * angular_vel[0] - dpseudo_n_d2t[1];
    Real angle_vel_dt_0 = sin_rotation_0 * rotation_0_a + cos_rotation_0 * rotation_0_b;

    Real rotation_1_a = dpseudo_n_d2t[0] * cos_rotation_1 - dpseudo_n_d2t[2] * sin_rotation_1 +
                        2.0 * angular_vel[1] * angular_vel[0] * sin_rotation_0;
    Real rotation_1_b1 = dpseudo_n_d2t[0] * cos_rotation_0 +
                         angular_vel[1] * angular_vel[1] * cos_rotation_0 * cos_rotation_0 * sin_rotation_1 +
                         angular_vel[0] * angular_vel[0] * sin_rotation_1 - dpseudo_n_d2t[1] * sin_rotation_1 * sin_rotation_0 +
                         2.0 * angular_vel[1] * angular_vel[0] * cos_rotation_1 * cos_rotation_0 * sin_rotation_0;
    Real rotation_1_b2 = -(dpseudo_n_d2t[2] * cos_rotation_0 +
                           angular_vel[1] * angular_vel[1] * cos_rotation_1 * cos_rotation_0 * cos_rotation_0 + angular_vel[0] * angular_vel[0] * cos_rotation_1 -
                           dpseudo_n_d2t[1] * cos_rotation_1 * sin_rotation_0 -
                           2.0 * angular_vel[1] * angular_vel[0] * cos_rotation_0 * sin_rotation_1 * sin_rotation_0);
    Real angle_vel_dt_1 = rotation_1_a * rotation_1_a * (rotation_1_b1 * cos_rotation_1 + rotation_1_b2 * sin_rotation_1) /
                          (rotation_1_b1 * rotation_1_b1 + rotation_1_b2 * rotation_1_b2 + Eps);

    return Vec3d(angle_vel_dt_0, angle_vel_dt_1, 0.0);
}
//=================================================================================================//

Vec3d getRotationFromPseudoNormalForSmallDeformation(
    const Vec3d &dpseudo_b_n_d2t, const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
{
    return Vec3d(0.0, dpseudo_n_d2t[0], 0.0);
}
//=================================================================================================//

Vec3d getRotationFromPseudoNormalForSmallDeformation_b(
    const Vec3d &dpseudo_b_n_d2t, const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
{
    return Vec3d(0.0, 0.0, dpseudo_b_n_d2t[0]);
}
//=================================================================================================//

Vec3d getNormalFromDeformationGradientTensor(const Mat3d &F)
{
    return F.col(0).cross(F.col(1)).normalized();
}
Vec3d getBinormalFromDeformationGradientTensor(const Mat3d &F)
{
    return (F.col(2).cross(F.col(0))).normalized();
}
//=================================================================================================//

Mat3d getCorrectedAlmansiStrain(const Mat3d &current_local_almansi_strain, const Real &nu_)
{
    Mat3d corrected_almansi_strain = current_local_almansi_strain;
    corrected_almansi_strain(2, 2) = -nu_ * (current_local_almansi_strain(0, 0) + current_local_almansi_strain(1, 1)) / (1.0 - nu_);
    return corrected_almansi_strain;
}
//=================================================================================================//

Mat3d getCorrectionMatrix(const Mat3d &local_deformation_part_one)
{
    Mat3d correction_matrix = Mat3d::Zero();
    correction_matrix.block<2, 2>(0, 0) = local_deformation_part_one.block<2, 2>(0, 0).inverse();
    return correction_matrix;
}

Mat3d getCorrectionMatrix_beam(const Mat3d &local_deformation_part_one)
{
    Mat3d correction_matrix = Mat3d::Zero();
    correction_matrix.block<1, 1>(0, 0) = local_deformation_part_one.block<1, 1>(0, 0).inverse();
    return correction_matrix;
}
//=================================================================================================//
} // namespace slender_structure_dynamics
} // namespace SPH
