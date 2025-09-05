#include "sphinxsys.h"

using namespace SPH;

namespace SPH::slender_structure_dynamics
{
Mat3d get_skew_matrix(const Vec3d &v)
{
    return Mat3d{
        {0.0, -v.z(), v.y()},
        {v.z(), 0.0, -v.x()},
        {-v.y(), v.x(), 0.0}};
}

// Update the transformation matrix g_i = lambada * g0_i based on the incremental rotation
// lambda^n+1 = exp(dtheta) * lambda^n
Mat3d get_lambda(const Mat3d &lambda, const Vec3d &dtheta)
{
    Real min = std::numeric_limits<Real>::epsilon();
    Real dtheta_norm = dtheta.norm();
    Mat3d S_dtheta = get_skew_matrix(dtheta);
    double a1 = dtheta_norm > min ? sin(dtheta_norm) / dtheta_norm : 1.0;
    double a2 = dtheta_norm > min ? (1.0 - cos(dtheta_norm)) / dtheta_norm / dtheta_norm : 0.5;
    Mat3d exp_S_dtheta = Mat3d::Identity() + a1 * S_dtheta + a2 * S_dtheta * S_dtheta;
    return exp_S_dtheta * lambda;
}

// Get the T operator based on the incremental rotation
Mat3d get_T_operator(const Vec3d &dtheta)
{
    Real min = std::numeric_limits<Real>::epsilon();
    Real dtheta_norm = dtheta.norm();
    if (dtheta_norm < min)
        return Mat3d::Identity();
    double a1 = sin(dtheta_norm) / dtheta_norm;
    double a2 = (1.0 - cos(dtheta_norm)) / dtheta_norm / dtheta_norm;
    double a3 = (dtheta_norm - sin(dtheta_norm)) / (dtheta_norm * dtheta_norm * dtheta_norm);
    Mat3d S_dtheta = get_skew_matrix(dtheta);
    return a1 * Mat3d::Identity() + a2 * S_dtheta + a3 * S_dtheta * S_dtheta.transpose();
}

class SimoReissnerStressRelaxationFirstHalf : public BaseBarRelaxation
{
  private:
    Real G0_;
    bool hourglass_control_;
    Real hourglass_control_factor_;

    Real *mass_;
    Vec3d *b_n0_;
    Real inv_W0_;

    Mat3d *C_N_;               // diag[EA,GA2,GA3]
    Mat3d *C_M_;               // diag[GJ,EI2,EI3]
    Mat3d *I_rho_l_;           // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;              // rho0*A, the mass per unit length
    Vec3d *m_prior_;           // the external torque per unit length in global coordinates
    Vec3d *angular_acc_prior_; // the angular acceleration caused by external torque

    Vec3d *Tau_;   // material tensile-shear strain
    Vec3d *Kappa_; // material bending strain
    Vec3d *f_;     // resultant force per unit length in global coordinates
    Vec3d *m_;     // resultant moment per unit length in global coordinates
    Vec3d *dr_ds_; // dr_ds in global coordinates

  public:
    explicit SimoReissnerStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                                   bool hourglass_control = false,
                                                   Real hourglass_control_factor = 0.005)
        : BaseBarRelaxation(inner_relation),
          G0_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial()).ShearModulus()),
          hourglass_control_(hourglass_control),
          hourglass_control_factor_(hourglass_control_factor),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          b_n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialBinormalDirection", "BinormalDirection")),
          inv_W0_(1.0 / sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)),
          C_N_(particles_->registerStateVariableData<Matd>("TensileConstitutiveMatrix")),
          C_M_(particles_->registerStateVariableData<Matd>("BendingConstitutiveMatrix")),
          I_rho_l_(particles_->registerStateVariableData<Matd>("InertiaMatrix")),
          A_rho_(particles_->registerStateVariableData<Real>("RhoCrossSectionArea")),
          m_prior_(particles_->registerStateVariableData<Vecd>("ExternalTorquePerUnitLength")),
          angular_acc_prior_(particles_->registerStateVariableData<Vecd>("PriorAngularAcceleration")),
          Tau_(particles_->registerStateVariableData<Vecd>("MaterialTensileShearStrain")),
          Kappa_(particles_->registerStateVariableData<Vecd>("MaterialBendingStrain")),
          f_(particles_->registerStateVariableData<Vecd>("ResultantForce")),
          m_(particles_->registerStateVariableData<Vecd>("ResultantMoment")),
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")) {}

    void initialization(size_t index_i, Real dt = 0.0)
    {
        // compute spatial resultant force and moment in global coordinates
        Mat3d Q_t = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]).transpose();
        f_[index_i] = Q_t * C_N_[index_i] * Tau_[index_i];
        m_[index_i] = Q_t * C_M_[index_i] * Kappa_[index_i];
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Mat3d I_i_global = Q.transpose() * I_rho_l_[index_i] * Q;
        auto I_i_inverse = I_i_global.inverse();
        Vec3d h_i_global = I_i_global * angular_vel_[index_i];
        Vec3d w_cross_h = angular_vel_[index_i].cross(h_i_global);

        const Mat3d Q_0 = transformation_matrix0_[index_i];

        Mat3d grad_f = Mat3d::Zero();
        Mat3d grad_m = Mat3d::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vec3d grad_W_ijV_j = inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ij_[n] * Vol_[index_j];
            // use symmetry scheme for stress and moment
            grad_f += (f_[index_i] + f_[index_j]) * grad_W_ijV_j.transpose();
            grad_m += (m_[index_i] + m_[index_j]) * grad_W_ijV_j.transpose();
        }
        Vec3d df_ds = (grad_f * Q_0.transpose()).col(xAxis);
        Vec3d dm_ds = (grad_m * Q_0.transpose()).col(xAxis);

        if (hourglass_control_)
        {
            auto [pos_hourglass_force, angular_momentum_hourglass_torque] = compute_hourglass_control_force(index_i);
            df_ds += pos_hourglass_force;
            dm_ds += angular_momentum_hourglass_torque;
        }

        force_[index_i] = df_ds / A_rho_[index_i] * mass_[index_i];
        dangular_vel_dt_[index_i] = I_i_inverse * (dm_ds + dr_ds_[index_i].cross(f_[index_i]) - w_cross_h);

        angular_acc_prior_[index_i] = I_i_inverse * m_prior_[index_i];
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
        angular_vel_[index_i] += (dangular_vel_dt_[index_i] + angular_acc_prior_[index_i]) * dt;
    }

  private:
    std::pair<Vec3d, Vec3d> compute_hourglass_control_force(size_t index_i)
    {
        Vec3d pos_hourglass = Vec3d::Zero();
        Vec3d angular_momentum_hourglass = Vec3d::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd e_ij = inner_neighborhood.e_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
            Vecd pos_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij, r_ij, pos_[index_i], F_[index_i], pos_[index_j], F_[index_j]);
            Real limiter_pos = SMIN(2.0 * pos_jump.norm() / r_ij, 1.0);
            pos_hourglass += hourglass_control_factor_ * weight * G0_ * pos_jump * Dimensions *
                             inner_neighborhood.dW_ij_[n] * Vol_[index_j] * limiter_pos;

            Vecd pseudo_n_variation_i = pseudo_n_[index_i] - n0_[index_i];
            Vecd pseudo_n_variation_j = pseudo_n_[index_j] - n0_[index_j];
            Vecd pseudo_n_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij, r_ij, pseudo_n_variation_i, F_bending_[index_i], pseudo_n_variation_j, F_bending_[index_j]);
            Real limiter_pseudo_n = SMIN(
                2.0 * pseudo_n_jump.norm() / ((pseudo_n_variation_i - pseudo_n_variation_j).norm() + Eps), 1.0);
            Vec3d n_hourglass = hourglass_control_factor_ * weight * G0_ * pseudo_n_jump * Dimensions *
                                inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                                pow(thickness_[index_i], 3) * limiter_pseudo_n;

            Vecd pseudo_b_variation_i = pseudo_b_n_[index_i] - b_n0_[index_i];
            Vecd pseudo_b_variation_j = pseudo_b_n_[index_j] - b_n0_[index_j];
            Vecd pseudo_b_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij, r_ij, pseudo_b_variation_i, F_b_bending_[index_i], pseudo_b_variation_j, F_b_bending_[index_j]);
            Real limiter_pseudo_b = SMIN(
                2.0 * pseudo_b_jump.norm() / ((pseudo_b_variation_i - pseudo_b_variation_j).norm() + Eps), 1.0);
            Vec3d b_hourglass = hourglass_control_factor_ * weight * G0_ * pseudo_b_jump * Dimensions *
                                inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                                pow(width_[index_i], 3) * limiter_pseudo_b;

            angular_momentum_hourglass +=
                pseudo_n_[index_i].cross(n_hourglass) + pseudo_b_n_[index_i].cross(b_hourglass);
        }
        return std::make_pair(pos_hourglass, angular_momentum_hourglass);
    }
};

class SimoReissnerStressRelaxationSecondHalf : public BaseBarRelaxation
{
  private:
    Vec3d *b_n0_;
    Matd *lambda_;     // transformation matrix from initial local to current local coordinates
    Vec3d *Tau_;       // material tensile-shear strain
    Vec3d *Kappa_;     // material bending strain
    Vec3d *dkappa_dt_; // spatial bending strain rate
    Vec3d *dr_ds_;     // dr_ds in global coordinates
    Vec3d *dv_ds_;     // dv_ds in global coordinates
    Mat3d *initial_curvature_;

  public:
    explicit SimoReissnerStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation),
          b_n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialBinormalDirection", "BinormalDirection")),
          lambda_(particles_->registerStateVariableData<Mat3d>("TransformationFromInitialToCurrent", IdentityMatrix<Matd>::value)),
          Tau_(particles_->registerStateVariableData<Vec3d>("MaterialTensileShearStrain")),
          Kappa_(particles_->registerStateVariableData<Vec3d>("MaterialBendingStrain")),
          dkappa_dt_(particles_->registerStateVariableData<Vec3d>("SpatialBendingStrainRate")),
          dr_ds_(particles_->registerStateVariableData<Vec3d>("DrDs")),
          dv_ds_(particles_->registerStateVariableData<Vec3d>("DvDs")),
          initial_curvature_(particles_->registerStateVariableData<Mat3d>("InitialCurvature")) {};

    void initialization(size_t index_i, Real dt = 0.0)
    {
        pos_[index_i] += vel_[index_i] * dt;
        Vec3d dtheta = angular_vel_[index_i] * dt;
        rotation_[index_i] += dtheta;
        lambda_[index_i] = get_lambda(lambda_[index_i], dtheta);
        pseudo_b_n_[index_i] = lambda_[index_i] * b_n0_[index_i];
        pseudo_n_[index_i] = lambda_[index_i] * n0_[index_i];
    }
    void interaction(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];

        Mat3d grad_vel = Matd::Zero();
        Mat3d grad_dangular_vel = Matd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            grad_vel -= (vel_[index_i] - vel_[index_j]) * gradW_ijV_j.transpose();
            grad_dangular_vel -= (angular_vel_[index_i] - angular_vel_[index_j]) * gradW_ijV_j.transpose();
        }
        // in initial local coordinate
        grad_vel = Q0_i * grad_vel * Q0_i.transpose() * B_[index_i];
        grad_dangular_vel = Q0_i * grad_dangular_vel * Q0_i.transpose() * B_[index_i];

        // tranform back to global
        dv_ds_[index_i] = Q0_i.transpose() * grad_vel.col(xAxis);
        const Vec3d dw_ds = Q0_i.transpose() * grad_dangular_vel.col(xAxis);

        // compute dkappa_dt
        Vec3d dtheta = angular_vel_[index_i] * dt;
        Mat3d T_operator = get_T_operator(dtheta);
        dkappa_dt_[index_i] = T_operator.transpose() * dw_ds;
    }
    void update(size_t index_i, Real dt = 0.0)
    {
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        dr_ds_[index_i] += dv_ds_[index_i] * dt;
        Tau_[index_i] = Q * dr_ds_[index_i] - Vec3d::UnitX();
        Kappa_[index_i] += Q * dkappa_dt_[index_i] * dt;

        // update F and F_bending in global coordinates
        const auto &Q0_i = transformation_matrix0_[index_i];
        F_[index_i].col(xAxis) = dr_ds_[index_i];
        F_bending_[index_i].col(yAxis) = pseudo_b_n_[index_i];
        F_bending_[index_i].col(zAxis) = pseudo_n_[index_i];
        F_[index_i] = F_[index_i] * Q0_i;

        Mat3d S_K = get_skew_matrix(Kappa_[index_i]);
        Mat3d dlambda_ds = Q.transpose() * (S_K + initial_curvature_[index_i]);
        Vec3d dbn_ds = dlambda_ds.col(yAxis);
        Vec3d dn_ds = dlambda_ds.col(zAxis);
        F_b_bending_[index_i] = dbn_ds * Q0_i.row(xAxis);
        F_bending_[index_i] = dn_ds * Q0_i.row(xAxis);
    }
};

class BeamInitialGeometry : public LocalDynamics, public DataDelegateInner
{
  private:
    Vec3d *dr0_ds_;            // initial dr_ds (tangential direction)
    Vec3d *dr_ds_;             // current dr_ds (tangential direction)
    Mat3d *initial_curvature_; // lambda0.t * d(lambda0)/ds

    Vec3d *pos0_; // initial position
    Real *Vol_;
    Mat3d *transformation_matrix0_;
    Mat3d *B_;

  public:
    explicit BeamInitialGeometry(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          dr0_ds_(particles_->registerStateVariableData<Vec3d>("InitialDrDs")),
          dr_ds_(particles_->registerStateVariableData<Vec3d>("DrDs")),
          initial_curvature_(particles_->registerStateVariableData<Mat3d>("InitialCurvature")),
          pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          transformation_matrix0_(particles_->getVariableDataByName<Mat3d>("TransformationMatrix")),
          B_(particles_->getVariableDataByName<Mat3d>("LinearGradientCorrectionMatrix")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];
        Mat3d grad_r = Matd::Zero();
        Mat3d grad_g1 = Matd::Zero();
        Mat3d grad_b_n = Matd::Zero();
        Mat3d grad_n = Matd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            grad_r -= (pos0_[index_i] - pos0_[index_j]) * gradW_ijV_j.transpose();
            const auto &Q0_j = transformation_matrix0_[index_j];
            grad_g1 -= (Q0_i.row(xAxis).transpose() - Q0_j.row(xAxis).transpose()) * gradW_ijV_j.transpose();
            grad_b_n -= (Q0_i.row(yAxis).transpose() - Q0_j.row(yAxis).transpose()) * gradW_ijV_j.transpose();
            grad_n -= (Q0_i.row(zAxis).transpose() - Q0_j.row(zAxis).transpose()) * gradW_ijV_j.transpose();
        }

        grad_r = Q0_i * grad_r * Q0_i.transpose() * B_[index_i];
        grad_b_n = Q0_i * grad_b_n * Q0_i.transpose() * B_[index_i];
        grad_n = Q0_i * grad_n * Q0_i.transpose() * B_[index_i];

        // tranform back to global
        dr0_ds_[index_i] = Q0_i.transpose() * grad_r.col(xAxis);
        dr_ds_[index_i] = dr0_ds_[index_i];

        Vec3d dg1_ds = Q0_i.transpose() * grad_g1.col(xAxis);
        Vec3d db_n0_ds = Q0_i.transpose() * grad_b_n.col(xAxis);
        Vec3d dn0_ds = Q0_i.transpose() * grad_n.col(xAxis);
        Mat3d dlamba0_ds{};
        dlamba0_ds.col(xAxis) = dg1_ds;
        dlamba0_ds.col(yAxis) = db_n0_ds;
        dlamba0_ds.col(zAxis) = dn0_ds;
        initial_curvature_[index_i] = Q0_i * dlamba0_ds;
    }
};
} // namespace SPH::slender_structure_dynamics