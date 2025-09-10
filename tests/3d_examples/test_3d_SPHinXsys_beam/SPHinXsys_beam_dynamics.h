#include "sphinxsys.h"

using namespace SPH;

namespace SPH::slender_structure_dynamics
{
struct GaussianQuadrature
{
    size_t number_of_gaussian_points_ = 0;
    std::vector<Real> gaussian_points_ = {};
    std::vector<Real> gaussian_weights_ = {};

    void reset(int number_of_points)
    {
        switch (number_of_points)
        {
        case 1:
            number_of_gaussian_points_ = 1;
            gaussian_points_ = {0.0};
            gaussian_weights_ = {2.0};
            break;
        case 2:
            number_of_gaussian_points_ = 2;
            gaussian_points_ = []()
            {
                Real point = 1.0 / sqrt(3.0);
                return std::vector<Real>{-point, point};
            }();
            gaussian_weights_ = {1.0, 1.0};
            break;
        case 3:
            number_of_gaussian_points_ = 3;
            gaussian_points_ = []()
            {
                Real point = sqrt(3.0 / 5.0);
                return std::vector<Real>{-point, 0.0, point};
            }();
            gaussian_weights_ = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
            break;
        case 4:
            number_of_gaussian_points_ = 4;
            gaussian_points_ = []()
            {
                Real point_1 = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
                Real point_2 = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
                return std::vector<Real>{-point_1, -point_2, point_2, point_1};
            }();
            gaussian_weights_ = []()
            {
                Real weight_1 = (18.0 + sqrt(30.0)) / 36.0;
                Real weight_2 = (18.0 - sqrt(30.0)) / 36.0;
                return std::vector<Real>{weight_1, weight_2, weight_2, weight_1};
            }();
            break;
        case 5:
            number_of_gaussian_points_ = 5;
            gaussian_points_ = []()
            {
                Real point_1 = 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
                Real point_2 = 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
                return std::vector<Real>{-point_1, -point_2, 0.0, point_2, point_1};
            }();
            gaussian_weights_ = []()
            {
                Real weight_1 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
                Real weight_2 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
                return std::vector<Real>{weight_1, weight_2, 128.0 / 225.0, weight_2, weight_1};
            }();
            break;
        default:
            throw std::invalid_argument("Unsupported number of Gaussian points");
        }
    }
};

Vec3d get_angular_acc_from_pseudo_normal_acc(
    const Vec3d &local_dnb_dt2,
    const Vec3d &local_dn_dt2,
    const Vec3d &rotation,
    const Vec3d &angular_vel)
{
    Real cos1 = cos(rotation.x());
    Real sin1 = sin(rotation.x());
    Real cos2 = cos(rotation.y());
    Real sin2 = sin(rotation.y());
    Real cos3 = cos(rotation.z());
    Real sin3 = sin(rotation.z());
    Eigen::Matrix<Real, 6, 3> A{
        {cos3 * sin2 * cos1 + sin3 * sin1, cos3 * cos2 * sin1, -sin3 * sin2 * sin1 - cos3 * cos1},
        {sin3 * sin2 * cos1 - cos3 * sin1, sin3 * cos2 * sin1, cos3 * sin2 * sin1 - sin3 * cos1},
        {cos2 * cos1, -sin2 * sin1, 0},
        {-cos3 * sin2 * sin1 + sin3 * cos1, cos3 * cos2 * cos1, -sin3 * sin2 * cos1 + cos3 * sin1},
        {-sin3 * sin2 * sin1 - cos3 * cos1, sin3 * cos2 * cos1, cos3 * sin2 * cos1 + sin3 * sin1},
        {-cos2 * sin1, -sin2 * cos1, 0}};

    Real w_norm = angular_vel.squaredNorm();
    Real w1_2 = angular_vel.x() * angular_vel.x();
    Real w2_2 = angular_vel.y() * angular_vel.y();
    Real w3_2 = angular_vel.z() * angular_vel.z();
    Real w1_w2 = angular_vel.x() * angular_vel.y();
    Real w2_w3 = angular_vel.y() * angular_vel.z();
    Real w1_w3 = angular_vel.x() * angular_vel.z();
    Real w1_2pw3_2 = w1_2 + w3_2;

    Real b1 = local_dnb_dt2.x() + cos3 * sin2 * sin1 * w_norm + 2 * (sin3 * cos2 * sin1 * w2_w3 + sin3 * sin2 * cos1 * w1_w3 - cos3 * cos2 * cos1 * w1_w2) - (sin3 * cos1 * w1_2pw3_2 + 2 * cos3 * sin1 * w1_w3);
    Real b2 = local_dnb_dt2.y() + sin3 * sin2 * sin1 * w_norm - 2 * (cos3 * cos2 * sin1 * w2_w3 + cos3 * sin2 * cos1 * w1_w3 - sin3 * cos2 * cos1 * w1_w2) + (cos3 * cos1 * w1_2pw3_2 - 2 * sin3 * sin1 * w1_w3);
    Real b3 = local_dnb_dt2.z() + cos2 * sin1 * (w1_2 + w2_2) + 2 * sin2 * cos1 * w1_w2;

    Real n1 = local_dn_dt2.x() + cos3 * sin2 * cos1 * w_norm + 2 * (sin3 * cos2 * cos1 * w2_w3 - sin3 * sin2 * sin1 * w1_w3 + cos3 * cos2 * sin1 * w1_w2) + (sin3 * sin1 * w1_2pw3_2 - 2 * cos3 * cos1 * w1_w3);
    Real n2 = local_dn_dt2.y() + sin3 * sin2 * cos1 * w_norm - 2 * (cos3 * cos2 * cos1 * w2_w3 - cos3 * sin2 * sin1 * w1_w3 - sin3 * cos2 * sin1 * w1_w2) - (cos3 * sin1 * w1_2pw3_2 + 2 * sin3 * cos1 * w1_w3);
    Real n3 = local_dn_dt2.z() + cos2 * cos1 * (w1_2 + w2_2) - 2 * sin2 * sin1 * w1_w2;

    Eigen::Matrix<Real, 6, 1> b{b1, b2, b3, n1, n2, n3};

    return A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}

Vec3d get_dn_from_angular_acc(const Vec3d &rotation, const Vec3d &angular_vel)
{
    Real cos1 = cos(rotation.x());
    Real sin1 = sin(rotation.x());
    Real cos2 = cos(rotation.y());
    Real sin2 = sin(rotation.y());
    Real cos3 = cos(rotation.z());
    Real sin3 = sin(rotation.z());
    Real dn1 = -sin3 * sin2 * cos1 * angular_vel.z() + cos3 * cos2 * cos1 * angular_vel.y() - cos3 * sin2 * sin1 * angular_vel.x() + cos3 * sin1 * angular_vel.z() + sin3 * cos1 * angular_vel.x();
    Real dn2 = cos3 * sin2 * cos1 * angular_vel.z() + sin3 * cos2 * cos1 * angular_vel.y() - sin3 * sin2 * sin1 * angular_vel.x() + sin3 * sin1 * angular_vel.z() - cos3 * cos1 * angular_vel.x();
    Real dn3 = -sin2 * cos1 * angular_vel.y() - cos2 * sin1 * angular_vel.x();
    return Vec3d{dn1, dn2, dn3};
}

Vec3d get_db_from_angular_acc(const Vec3d &rotation, const Vec3d &angular_vel)
{
    Real cos1 = cos(rotation.x());
    Real sin1 = sin(rotation.x());
    Real cos2 = cos(rotation.y());
    Real sin2 = sin(rotation.y());
    Real cos3 = cos(rotation.z());
    Real sin3 = sin(rotation.z());
    Real db1 = -sin3 * sin2 * sin1 * angular_vel.z() + cos3 * cos2 * sin1 * angular_vel.y() + cos3 * sin2 * cos1 * angular_vel.x() - cos3 * cos1 * angular_vel.z() + sin3 * sin1 * angular_vel.x();
    Real db2 = cos3 * sin2 * sin1 * angular_vel.z() + sin3 * cos2 * sin1 * angular_vel.y() + sin3 * sin2 * cos1 * angular_vel.x() - sin3 * cos1 * angular_vel.z() - cos3 * sin1 * angular_vel.x();
    Real db3 = -sin2 * sin1 * angular_vel.y() + cos2 * cos1 * angular_vel.x();
    return Vec3d{db1, db2, db3};
}

std::pair<Vec3d, Vec3d>
get_torque_components(const Vec3d &m, const Vec3d &g1, const Vec3d &g2, const Vec3d &g3, Real alpha = 0.5)
{
    Real m1 = m.dot(g1);
    Real m2 = m.dot(g2);
    Real m3 = m.dot(g3);
    Vec3d my = -m3 * g1 + alpha * m1 * g3;
    Vec3d mz = m2 * g1 - (1 - alpha) * m1 * g2;
    return std::make_pair(my, mz);
}

class BeamStressRelaxationFirstHalf : public BaseBarRelaxation
{
  private:
    ElasticSolid &elastic_solid_;
    Real rho0_;
    Real inv_rho0_;
    Real smoothing_length_;
    Mat3d numerical_damping_scaling_matrix_;
    Real *rho_;
    Real *mass_;
    Vecd *b_n0_;
    Mat3d *global_stress_;
    Mat3d *global_moment_;
    Mat3d *global_b_moment_;
    Vec3d *global_shear_stress_;
    Vec3d *global_b_shear_stress_;
    Mat3d *global_F_;
    Mat3d *global_F_bending_;
    Mat3d *global_F_b_bending_;

    Real G0_;
    Real hourglass_control_factor_;
    bool hourglass_control_;
    const Real inv_W0_ = 1.0 / sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd);

    Vecd *m_prior_;             // the external torque per unit length in global coordinates
    Vecd *angular_acc_prior_L_; // the angular acceleration caused by external torque

    GaussianQuadrature gaussian_{};

  public:
    explicit BeamStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                           size_t number_of_gaussian_points = 3,
                                           bool hourglass_control = false,
                                           Real hourglass_control_factor = 0.005)
        : BaseBarRelaxation(inner_relation),
          elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial())),
          rho0_(elastic_solid_.ReferenceDensity()), inv_rho0_(1.0 / rho0_),
          smoothing_length_(sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
          numerical_damping_scaling_matrix_(Matd::Identity() * smoothing_length_),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          b_n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialBinormalDirection", "BinormalDirection")),
          global_stress_(particles_->registerStateVariableData<Mat3d>("GlobalStress")),
          global_moment_(particles_->registerStateVariableData<Mat3d>("GlobalMoment")),
          global_b_moment_(particles_->registerStateVariableData<Matd>("GlobalBinormalMoment")),
          global_shear_stress_(particles_->registerStateVariableData<Vec3d>("GlobalShearStress")),
          global_b_shear_stress_(particles_->registerStateVariableData<Vecd>("GlobalBinormalShearStress")),
          global_F_(particles_->registerStateVariableData<Mat3d>("GlobalDeformationGradient")),
          global_F_bending_(particles_->registerStateVariableData<Mat3d>("GlobalBendingDeformationGradient")),
          global_F_b_bending_(particles_->registerStateVariableData<Mat3d>("GlobalBinormalBendingDeformationGradient")),
          G0_(elastic_solid_.ShearModulus()),
          hourglass_control_(hourglass_control),
          m_prior_(particles_->registerStateVariableData<Vec3d>("ExternalTorquePerUnitLength")),
          angular_acc_prior_L_(particles_->registerStateVariableData<Vec3d>("PriorAngularAcceleration"))
    {
        gaussian_.reset(number_of_gaussian_points);
    }

    void initialization(size_t index_i, Real dt = 0.0)
    {
        // update position and rotation
        pos_[index_i] += vel_[index_i] * dt * 0.5;
        rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
        pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
        pseudo_b_n_[index_i] += dpseudo_b_n_dt_[index_i] * dt * 0.5;

        F_[index_i] += dF_dt_[index_i] * dt * 0.5;
        F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
        F_b_bending_[index_i] += dF_b_bending_dt_[index_i] * dt * 0.5;

        Real J = F_[index_i].determinant();
        rho_[index_i] = rho0_ / J;

        // compute global F
        const Mat3d &Q_0 = transformation_matrix0_[index_i];
        global_F_[index_i] = Q_0.transpose() * F_[index_i] * Q_0;
        global_F_bending_[index_i] = Q_0.transpose() * F_bending_[index_i] * Q_0;
        global_F_b_bending_[index_i] = Q_0.transpose() * F_b_bending_[index_i] * Q_0;

        // compute resultant stress and moment
        auto [resultant_stress_L, resultant_moment_L] = compute_resultant_stress(index_i);

        Mat3d stress_L = Mat3d::Zero();
        stress_L.col(xAxis) = resultant_stress_L.col(xAxis);
        global_stress_[index_i] = Q_0.transpose() * (stress_L * B_[index_i]) * Q_0;

        Mat3d moment_b_L = Mat3d::Zero();
        moment_b_L.col(xAxis) = resultant_moment_L.col(yAxis);
        global_b_moment_[index_i] = Q_0.transpose() * (moment_b_L * B_[index_i]) * Q_0;

        Mat3d moment_L = Mat3d::Zero();
        moment_L.col(xAxis) = resultant_moment_L.col(zAxis);
        global_moment_[index_i] = Q_0.transpose() * (moment_L * B_[index_i]) * Q_0;

        global_b_shear_stress_[index_i] = Q_0.transpose() * resultant_stress_L.col(yAxis);
        global_shear_stress_[index_i] = Q_0.transpose() * resultant_stress_L.col(zAxis);
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Mat3d Q_0_i = transformation_matrix0_[index_i];

        Vec3d df_ds = Vec3d::Zero();
        Vec3d dmy_ds = Vec3d::Zero();
        Vec3d dmz_ds = Vec3d::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vec3d grad_W_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            // use symmetry scheme for stress and moment
            df_ds += (global_stress_[index_i] + global_stress_[index_j]) * grad_W_ijV_j;
            dmy_ds += (global_b_moment_[index_i] + global_b_moment_[index_j]) * grad_W_ijV_j;
            dmz_ds += (global_moment_[index_i] + global_moment_[index_j]) * grad_W_ijV_j;
        }

        if (hourglass_control_)
        {
            auto [pos_hourglass, b_hourglass, n_hourglass] = compute_hourglass_control_force(index_i);
            df_ds += pos_hourglass;
            dmy_ds += b_hourglass;
            dmz_ds += n_hourglass;
        }

        force_[index_i] =
            mass_[index_i] * df_ds * inv_rho0_ / (thickness_[index_i] * width_[index_i]);
        Real Iy_rho_inv = 12.0 / width_[index_i] / pow(thickness_[index_i], 3) * inv_rho0_;
        Real Iz_rho_inv = 12.0 / thickness_[index_i] / pow(width_[index_i], 3) * inv_rho0_;
        dpseudo_n_d2t_[index_i] = Iy_rho_inv * (dmz_ds - global_shear_stress_[index_i]);
        dpseudo_b_n_d2t_[index_i] = Iz_rho_inv * (dmy_ds - global_b_shear_stress_[index_i]);

        Vec3d local_dpseudo_n_d2t = Q_0_i * dpseudo_n_d2t_[index_i];
        Vec3d local_dpseudo_b_n_d2t = Q_0_i * dpseudo_b_n_d2t_[index_i];

        dangular_vel_dt_[index_i] = get_angular_acc_from_pseudo_normal_acc(
            local_dpseudo_b_n_d2t, local_dpseudo_n_d2t, rotation_[index_i], angular_vel_[index_i]);

        // torque
        if (m_prior_[index_i].norm() < Eps)
            return;
        auto [my, mz] = get_torque_components(
            m_prior_[index_i],
            pseudo_b_n_[index_i].cross(pseudo_n_[index_i]),
            pseudo_b_n_[index_i],
            pseudo_n_[index_i]);
        Vec3d dpseudo_n_d2t_ext = Iy_rho_inv * mz;
        Vec3d dpseudo_b_n_d2t_ext = Iz_rho_inv * my;
        Vec3d local_dpseudo_n_d2t_ext = Q_0_i * dpseudo_n_d2t_ext;
        Vec3d local_dpseudo_b_n_d2t_ext = Q_0_i * dpseudo_b_n_d2t_ext;
        angular_acc_prior_L_[index_i] = get_angular_acc_from_pseudo_normal_acc(
            local_dpseudo_b_n_d2t_ext, local_dpseudo_n_d2t_ext, rotation_[index_i], angular_vel_[index_i]);
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
        angular_vel_[index_i] += (dangular_vel_dt_[index_i] + angular_acc_prior_L_[index_i]) * dt;
    }

  private:
    std::tuple<Vec3d, Vec3d, Vec3d> compute_hourglass_control_force(size_t index_i)
    {
        Mat3d Q_0_i = transformation_matrix0_[index_i];

        Vec3d pos_hourglass = Vec3d::Zero();
        Vec3d b_n_hourglass = Vec3d::Zero();
        Vec3d n_hourglass = Vec3d::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Mat3d Q_0_j = transformation_matrix0_[index_j];
            Vecd e_ij = inner_neighborhood.e_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
            Vecd pos_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij,
                r_ij,
                pos_[index_i],
                global_F_[index_i],
                pos_[index_j],
                global_F_[index_j]);
            Real limiter_pos = SMIN(2.0 * pos_jump.norm() / r_ij, 1.0);
            pos_hourglass += hourglass_control_factor_ * weight * G0_ * pos_jump * Dimensions *
                             inner_neighborhood.dW_ij_[n] * Vol_[index_j] * limiter_pos;

            Vecd pseudo_n_variation_i = pseudo_n_[index_i] - n0_[index_i];
            Vecd pseudo_n_variation_j = pseudo_n_[index_j] - n0_[index_j];
            Vecd pseudo_n_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij,
                r_ij,
                pseudo_n_variation_i,
                global_F_bending_[index_i],
                pseudo_n_variation_j,
                global_F_bending_[index_j]);
            Real limiter_pseudo_n = SMIN(
                2.0 * pseudo_n_jump.norm() / ((pseudo_n_variation_i - pseudo_n_variation_j).norm() + Eps), 1.0);
            n_hourglass += hourglass_control_factor_ * weight * G0_ * pseudo_n_jump * Dimensions *
                           inner_neighborhood.dW_ij_[n] * Vol_[index_j] * pow(thickness_[index_i], 3) * limiter_pseudo_n;

            Vecd pseudo_b_variation_i = pseudo_b_n_[index_i] - b_n0_[index_i];
            Vecd pseudo_b_variation_j = pseudo_b_n_[index_j] - b_n0_[index_j];
            Vecd pseudo_b_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij,
                r_ij,
                pseudo_b_variation_i,
                global_F_b_bending_[index_i],
                pseudo_b_variation_j,
                global_F_b_bending_[index_j]);
            Real limiter_pseudo_b = SMIN(
                2.0 * pseudo_b_jump.norm() / ((pseudo_b_variation_i - pseudo_b_variation_j).norm() + Eps), 1.0);
            b_n_hourglass += hourglass_control_factor_ * weight * G0_ * pseudo_b_jump * Dimensions *
                             inner_neighborhood.dW_ij_[n] * Vol_[index_j] * pow(width_[index_i], 3) * limiter_pseudo_b;
        }

        return std::make_tuple(pos_hourglass, b_n_hourglass, n_hourglass);
    }

    // For now we assume the bar has rectangular cross section.
    // Rectangular: 0.5 * width * point_y
    Real get_y_G(Real point_y, Real, Real width, Real) const { return 0.5 * width * point_y; }
    // Rectangular: 0.5 * thickness * point_z
    Real get_z_G(Real point_y, Real point_z, Real, Real thickness) const { return 0.5 * thickness * point_z; }
    // Rectangular: 0.25 * width * thickness
    Real get_J_G(Real point_y, Real point_z, Real width, Real thickness) const { return 0.25 * width * thickness; }
    std::pair<Mat3d, Mat3d> compute_resultant_stress(size_t index_i)
    {
        Real width = width_[index_i];
        Real thickness = thickness_[index_i];

        Mat3d resultant_stress_L = Mat3d::Zero();
        Mat3d resultant_moment_L = Mat3d::Zero();

        numerical_damping_scaling_matrix_(1, 1) = width_[index_i] < smoothing_length_ ? width_[index_i] : smoothing_length_;
        numerical_damping_scaling_matrix_(2, 2) = thickness_[index_i] < smoothing_length_ ? thickness_[index_i] : smoothing_length_;

        for (size_t i = 0; i != gaussian_.number_of_gaussian_points_; ++i)
        {
            Real weight_y = gaussian_.gaussian_weights_[i];
            Real point_y = gaussian_.gaussian_points_[i];

            for (size_t j = 0; j != gaussian_.number_of_gaussian_points_; ++j)
            {
                Real weight_z = gaussian_.gaussian_weights_[j];
                Real point_z = gaussian_.gaussian_points_[j];

                Real y_G = get_y_G(point_y, point_z, width, thickness);
                Real z_G = get_z_G(point_y, point_z, width, thickness);
                Real J_G = get_J_G(point_y, point_z, width, thickness);

                Matd F_gaussian_L = F_[index_i] + y_G * F_b_bending_[index_i] + z_G * F_bending_[index_i];
                // TODO: some material models have orientations, we have to use global F
                Mat3d P_gaussian_point_L = elastic_solid_.StressPK1(F_gaussian_L, index_i);

                // The numerical stress is S, need to convert it to P

                Matd dF_dt_gaussian_L =
                    dF_dt_[index_i] + y_G * dF_b_bending_dt_[index_i] + z_G * dF_bending_dt_[index_i];
                Mat3d damping_L =
                    F_gaussian_L * elastic_solid_.NumericalDampingRightCauchy(
                                       F_gaussian_L, dF_dt_gaussian_L, numerical_damping_scaling_matrix_, index_i);
                P_gaussian_point_L += damping_L;

                Vec3d T_L = P_gaussian_point_L.col(xAxis);

                /** Integrate Cauchy stress along thickness. */
                Real weight_yz = J_G * weight_y * weight_z;
                resultant_stress_L += weight_yz * P_gaussian_point_L;
                resultant_moment_L.col(yAxis) += weight_yz * y_G * T_L;
                resultant_moment_L.col(zAxis) += weight_yz * z_G * T_L;
            }
        }
        return {resultant_stress_L, resultant_moment_L};
    }
};

class BeamStressRelaxationSecondHalf : public BaseBarRelaxation
{
  public:
    explicit BeamStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation) {};

    void initialization(size_t index_i, Real dt = 0.0)
    {
        Mat3d Q0_t = transformation_matrix0_[index_i].transpose();
        pos_[index_i] += vel_[index_i] * dt * 0.5;
        rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
        Vec3d local_db_n_dt = get_db_from_angular_acc(rotation_[index_i], angular_vel_[index_i]);
        Vec3d local_dn_dt = get_dn_from_angular_acc(rotation_[index_i], angular_vel_[index_i]);
        dpseudo_n_dt_[index_i] = Q0_t * local_dn_dt;
        dpseudo_b_n_dt_[index_i] = Q0_t * local_db_n_dt;
        pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
        pseudo_b_n_[index_i] += dpseudo_b_n_dt_[index_i] * dt * 0.5;
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Matd Q0_i = transformation_matrix0_[index_i];

        Mat3d deformation_gradient_change_rate_part_one = Mat3d::Zero();
        Mat3d deformation_gradient_change_rate_part_two = Mat3d::Zero();
        Mat3d deformation_gradient_change_rate_part_three = Mat3d::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            deformation_gradient_change_rate_part_one -= (vel_[index_i] - vel_[index_j]) * gradW_ijV_j.transpose();
            deformation_gradient_change_rate_part_two -=
                (dpseudo_n_dt_[index_i] - dpseudo_n_dt_[index_j]) * gradW_ijV_j.transpose();
            deformation_gradient_change_rate_part_three -=
                (dpseudo_b_n_dt_[index_i] - dpseudo_b_n_dt_[index_j]) * gradW_ijV_j.transpose();
        }
        dF_dt_[index_i] = Q0_i * deformation_gradient_change_rate_part_one * Q0_i.transpose() * B_[index_i];
        dF_dt_[index_i].col(Dimensions - 1) = Q0_i * dpseudo_n_dt_[index_i];
        dF_dt_[index_i].col(Dimensions - 2) = Q0_i * dpseudo_b_n_dt_[index_i];

        dF_bending_dt_[index_i] = Q0_i * deformation_gradient_change_rate_part_two * Q0_i.transpose() * B_[index_i];
        dF_b_bending_dt_[index_i] = Q0_i * deformation_gradient_change_rate_part_three * Q0_i.transpose() * B_[index_i];
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        F_[index_i] += dF_dt_[index_i] * dt * 0.5;
        F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
        F_b_bending_[index_i] += dF_b_bending_dt_[index_i] * dt * 0.5;
    };
};
} // namespace SPH::slender_structure_dynamics