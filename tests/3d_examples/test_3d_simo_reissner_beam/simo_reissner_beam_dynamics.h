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
    Real a1 = dtheta_norm > min ? sin(dtheta_norm) / dtheta_norm : 1.0;
    Real a2 = dtheta_norm > min ? (1.0 - cos(dtheta_norm)) / dtheta_norm / dtheta_norm : 0.5;
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
    Real a1 = sin(dtheta_norm) / dtheta_norm;
    Real a2 = (1.0 - cos(dtheta_norm)) / dtheta_norm / dtheta_norm;
    Real a3 = (dtheta_norm - sin(dtheta_norm)) / (dtheta_norm * dtheta_norm * dtheta_norm);
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

    Vec3d *dTau_;   // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_; // Kappa^n+1 - Kappa^n, used for numerical damping
    Real numerical_damping_factor_ = 0.5;

    Mat3d *C_N_;               // diag[EA,GA2,GA3]
    Mat3d *C_M_;               // diag[GJ,EI2,EI3]
    Mat3d *I_rho_l_;           // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;              // rho0*A, the mass per unit length
    Vec3d *m_prior_;           // the external torque per unit length in global coordinates
    Vec3d *angular_acc_prior_; // the angular acceleration caused by external torque

    Vec3d *Tau_;        // material tensile-shear strain
    Vec3d *Kappa_;      // material bending strain
    Vec3d *f_;          // resultant force per unit length in global coordinates
    Vec3d *m_;          // resultant moment per unit length in global coordinates
    Vec3d *dr_ds_;      // dr_ds in global coordinates
    Mat3d *grad_r_;     // nabla_r in global coordinates
    Mat3d *grad_theta_; // nabla_theta in global coordinates

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
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")),
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
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")),
          grad_r_(particles_->registerStateVariableData<Mat3d>("NablaR")),
          grad_theta_(particles_->registerStateVariableData<Mat3d>("NablaTheta")) {}

    void initialization(size_t index_i, Real dt = 0.0)
    {
        // compute spatial resultant force and moment in global coordinates
        Mat3d Q_t = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]).transpose();
        Vec3d f_m = C_N_[index_i] * (Tau_[index_i] + numerical_damping_factor_ * dTau_[index_i]);
        Vec3d m_m = C_M_[index_i] * (Kappa_[index_i] + numerical_damping_factor_ * dKappa_[index_i]);
        f_[index_i] = Q_t * f_m;
        m_[index_i] = Q_t * m_m;
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
            auto [pos_hourglass_force, theta_hourglass_torque] = compute_hourglass_control_force(index_i);
            df_ds += pos_hourglass_force;
            dm_ds += theta_hourglass_torque;
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
        Vec3d theta_hourglass = Vec3d::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd e_ij = inner_neighborhood.e_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
            Vecd pos_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij, r_ij, pos_[index_i], grad_r_[index_i], pos_[index_j], grad_r_[index_j]);
            Real limiter_pos = SMIN(2.0 * pos_jump.norm() / r_ij, 1.0);
            pos_hourglass += hourglass_control_factor_ * weight * G0_ * pos_jump * Dimensions *
                             inner_neighborhood.dW_ij_[n] * Vol_[index_j] * limiter_pos;

            Vecd theta_jump = thin_structure_dynamics::getLinearVariableJump(
                e_ij, r_ij, rotation_[index_i], grad_theta_[index_i], rotation_[index_j], grad_theta_[index_j]);
            Real limiter_theta = SMIN(
                2.0 * theta_jump.norm() / ((rotation_[index_i] - rotation_[index_j]).norm() + Eps), 1.0);
            theta_hourglass += hourglass_control_factor_ * weight * G0_ * theta_jump * Dimensions *
                               inner_neighborhood.dW_ij_[n] * Vol_[index_j] * limiter_theta;
        }
        return std::make_pair(pos_hourglass, theta_hourglass);
    }
};

class SimoReissnerStressRelaxationSecondHalf : public BaseBarRelaxation
{
  private:
    Vec3d *b_n0_;
    Matd *lambda_;            // transformation matrix from initial local to current local coordinates
    Vec3d *Tau_;              // material tensile-shear strain
    Vec3d *Kappa_;            // material bending strain
    Vec3d *dkappa_dt_;        // spatial bending strain rate
    Vec3d *dr_ds_;            // dr_ds in global coordinates
    Vec3d *dv_ds_;            // dv_ds in global coordinates
    Mat3d *grad_r_;           // nabla_r in global coordinates
    Mat3d *grad_angular_vel_; // nabla_w in global coordinates
    Mat3d *grad_theta_;       // nabla_theta in global coordinates
    Vec3d *dTau_;             // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_;           // Kappa^n+1 - Kappa^n, used for numerical damping

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
          grad_r_(particles_->registerStateVariableData<Mat3d>("NablaR")),
          grad_angular_vel_(particles_->registerStateVariableData<Mat3d>("NablaW")),
          grad_theta_(particles_->registerStateVariableData<Mat3d>("NablaTheta")),
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")) {};

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
        grad_angular_vel_[index_i] = Q0_i * grad_dangular_vel * Q0_i.transpose() * B_[index_i];

        // tranform back to global
        dv_ds_[index_i] = Q0_i.transpose() * grad_vel.col(xAxis);
        Vec3d dw_ds = Q0_i.transpose() * grad_angular_vel_[index_i].col(xAxis);

        // compute dkappa_dt
        Vec3d dtheta = angular_vel_[index_i] * dt;
        Mat3d T_operator = get_T_operator(dtheta);
        dkappa_dt_[index_i] = T_operator.transpose() * dw_ds;
    }
    void update(size_t index_i, Real dt = 0.0)
    {
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        dr_ds_[index_i] += dv_ds_[index_i] * dt;
        Vec3d Tau_n = Tau_[index_i];
        Vec3d Kappa_n = Kappa_[index_i];
        Tau_[index_i] = Q * dr_ds_[index_i] - Vec3d::UnitX();
        Kappa_[index_i] += Q * dkappa_dt_[index_i] * dt;
        dTau_[index_i] = Tau_[index_i] - Tau_n;
        dKappa_[index_i] = Kappa_[index_i] - Kappa_n;

        const auto &Q0_i = transformation_matrix0_[index_i];
        grad_r_[index_i] = dr_ds_[index_i] * Q0_i.row(xAxis);
        grad_theta_[index_i] += grad_angular_vel_[index_i] * dt;
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

class VelocityUpdate : public LocalDynamics
{
  private:
    Real *mass_;
    Vec3d *vel_;
    Vec3d *angular_vel_;
    Vec3d *pseudo_n_;
    Vec3d *pseudo_b_n_;

    Mat3d *lambda_;       // transformation matrix from initial local to current local coordinates
    Mat3d *I_rho_l_;      // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;         // rho0*A, the mass per unit length
    Vec3d *force_;        // internal moment
    Vec3d *force_prior_;  // the external torque per unit length in global coordinates
    Vec3d *moment_;       // internal moment
    Vec3d *moment_prior_; // the external torque per unit length in global coordinates

    Vec3d *Tau_;   // material tensile-shear strain
    Vec3d *Kappa_; // material bending strain
    Vec3d *f_;     // resultant force per unit length in global coordinates
    Vec3d *m_;     // resultant moment per unit length in global coordinates
    Vec3d *dr_ds_; // dr_ds in global coordinates

    Vec3d *linear_momentum_;  // linear momentum
    Vec3d *angular_momentum_; // skew symmetric angular momentum
    Vec3d *dtheta_;           // incremental rotation

    Vec3d *dTau_;   // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_; // Kappa^n+1 - Kappa^n, used for numerical damping

    Mat3d *grad_r_;      // nabla_r in global coordinates
    Mat3d *grad_dtheta_; // nabla_dtheta in global coordinates

    Vec3d *force_prior_discretized_;
    Vec3d *moment_prior_discretized_;

  public:
    explicit VelocityUpdate(SPHBody &body)
        : LocalDynamics(body),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          angular_vel_(particles_->getVariableDataByName<Vecd>("AngularVelocity")),
          pseudo_n_(particles_->registerStateVariableDataFrom<Vecd>("PseudoNormal", "NormalDirection")),
          pseudo_b_n_(particles_->registerStateVariableDataFrom<Vecd>("PseudoBinormal", "BinormalDirection")),
          lambda_(particles_->registerStateVariableData<Mat3d>("TransformationFromInitialToCurrent", IdentityMatrix<Matd>::value)),
          I_rho_l_(particles_->registerStateVariableData<Matd>("InertiaMatrix")),
          A_rho_(particles_->registerStateVariableData<Real>("RhoCrossSectionArea")),
          force_(particles_->getVariableDataByName<Vecd>("Force")),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          moment_(particles_->registerStateVariableData<Vecd>("InternalMoment")),
          moment_prior_(particles_->registerStateVariableData<Vecd>("ExternalTorquePerUnitLength")),
          Tau_(particles_->registerStateVariableData<Vecd>("MaterialTensileShearStrain")),
          Kappa_(particles_->registerStateVariableData<Vecd>("MaterialBendingStrain")),
          f_(particles_->registerStateVariableData<Vecd>("ResultantForce")),
          m_(particles_->registerStateVariableData<Vecd>("ResultantMoment")),
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")),
          linear_momentum_(particles_->registerStateVariableData<Vecd>("LinearMomentum")),
          angular_momentum_(particles_->registerStateVariableData<Vecd>("AngularMomentum")),
          dtheta_(particles_->registerStateVariableData<Vecd>("IncrementalRotation")),
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")),
          grad_r_(particles_->registerStateVariableData<Mat3d>("NablaR")),
          grad_dtheta_(particles_->registerStateVariableData<Mat3d>("NabladTheta")),
          force_prior_discretized_(particles_->registerStateVariableData<Vecd>("ForcePriorDiscretized")),
          moment_prior_discretized_(particles_->registerStateVariableData<Vecd>("MomentPriorDiscretized")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] += 0.5 * dt * (force_prior_discretized_[index_i] + force_[index_i]) / mass_[index_i];
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Mat3d I_i_global = Q.transpose() * I_rho_l_[index_i] * Q;
        angular_momentum_[index_i] += 0.5 * dt * (moment_[index_i] + moment_prior_discretized_[index_i]);
        angular_vel_[index_i] = 0.5 * I_i_global.inverse() * angular_momentum_[index_i];
    }
};

class SimoReissnerRelaxationFirstHalf_v2 : public BaseBarRelaxation
{
  private:
    Real *mass_;
    Vec3d *b_n0_;

    Mat3d *lambda_;       // transformation matrix from initial local to current local coordinates
    Mat3d *I_rho_l_;      // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;         // rho0*A, the mass per unit length
    Vec3d *moment_;       // internal moment
    Vec3d *moment_prior_; // the external torque per unit length in global coordinates

    Vec3d *Tau_;   // material tensile-shear strain
    Vec3d *Kappa_; // material bending strain
    Vec3d *f_;     // resultant force per unit length in global coordinates
    Vec3d *m_;     // resultant moment per unit length in global coordinates
    Vec3d *dr_ds_; // dr_ds in global coordinates

    Vec3d *linear_momentum_;  // linear momentum
    Vec3d *angular_momentum_; // skew symmetric angular momentum
    Vec3d *dtheta_;           // incremental rotation

    Vec3d *dTau_;   // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_; // Kappa^n+1 - Kappa^n, used for numerical damping

    Mat3d *grad_r_;      // nabla_r in global coordinates
    Mat3d *grad_dtheta_; // nabla_dtheta in global coordinates

    Vec3d *force_prev_;  // force at previous time step
    Vec3d *moment_prev_; // moment at previous time step
    Vec3d *force_prior_prev_;
    Vec3d *moment_prior_prev_;

  public:
    explicit SimoReissnerRelaxationFirstHalf_v2(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          b_n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialBinormalDirection", "BinormalDirection")),
          lambda_(particles_->registerStateVariableData<Mat3d>("TransformationFromInitialToCurrent", IdentityMatrix<Matd>::value)),
          I_rho_l_(particles_->registerStateVariableData<Matd>("InertiaMatrix")),
          A_rho_(particles_->registerStateVariableData<Real>("RhoCrossSectionArea")),
          moment_(particles_->registerStateVariableData<Vecd>("InternalMoment")),
          moment_prior_(particles_->registerStateVariableData<Vecd>("ExternalTorquePerUnitLength")),
          Tau_(particles_->registerStateVariableData<Vecd>("MaterialTensileShearStrain")),
          Kappa_(particles_->registerStateVariableData<Vecd>("MaterialBendingStrain")),
          f_(particles_->registerStateVariableData<Vecd>("ResultantForce")),
          m_(particles_->registerStateVariableData<Vecd>("ResultantMoment")),
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")),
          linear_momentum_(particles_->registerStateVariableData<Vecd>("LinearMomentum")),
          angular_momentum_(particles_->registerStateVariableData<Vecd>("AngularMomentum")),
          dtheta_(particles_->registerStateVariableData<Vecd>("IncrementalRotation")),
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")),
          grad_r_(particles_->registerStateVariableData<Mat3d>("NablaR")),
          grad_dtheta_(particles_->registerStateVariableData<Mat3d>("NabladTheta")),
          force_prev_(particles_->registerStateVariableData<Vecd>("ForcePrev")),
          moment_prev_(particles_->registerStateVariableData<Vecd>("MomentPrev")),
          force_prior_prev_(particles_->registerStateVariableData<Vecd>("ForcePriorPrev")),
          moment_prior_prev_(particles_->registerStateVariableData<Vecd>("MomentPriorPrev")) {}

    void initialization(size_t index_i, Real dt = 0.0)
    {
        // record n step force and moment
        force_prev_[index_i] = force_[index_i];
        moment_prev_[index_i] = moment_[index_i];
        force_prior_prev_[index_i] = force_prior_[index_i];
        moment_prior_prev_[index_i] = moment_prior_[index_i];

        // recompute linear momentum and angular momentum based on the last step velocity and angular velocity
        // to apply velocity constraint
        linear_momentum_[index_i] = vel_[index_i] * A_rho_[index_i];
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Mat3d I_i_global = Q.transpose() * I_rho_l_[index_i] * Q;
        angular_momentum_[index_i] = 2 * I_i_global * angular_vel_[index_i];

        // update position and rotation matrix based on the last step
        pos_[index_i] += dt * vel_[index_i];

        Vec3d Y = dt * angular_momentum_[index_i];
        // use Newton - Raphson method to find the incremental rotation
        auto get_r_and_grad_r = [&](Vec3d theta_vec) -> std::pair<Vec3d, Mat3d>
        {
            Real theta = theta_vec.norm();
            if (theta < TinyReal)
                return std::make_pair(-Y, I_i_global);
            Vec3d i_rho_theta = I_i_global * theta_vec;
            Real sin_t = sin(theta);
            Real cos_t = cos(theta);
            Vec3d residual = sin_t / theta * i_rho_theta + (1 - cos_t) / theta / theta * (theta_vec.cross(i_rho_theta)) - Y;
            Mat3d grad_residual = sin_t / theta * I_i_global + (theta * cos_t - sin_t) / pow(theta, 3) * i_rho_theta * theta_vec.transpose() +
                                  (1 - cos_t) / theta / theta * (-get_skew_matrix(i_rho_theta) + get_skew_matrix(theta_vec) * I_i_global) +
                                  (theta * sin_t - 2 * (1 - cos_t)) / pow(theta, 4) * (theta_vec.cross(i_rho_theta)) * theta_vec.transpose();
            return std::make_pair(residual, grad_residual);
        };
        constexpr int max_iterations = 20; // @WARN hard-coded maximum number of iterations
        constexpr auto infinity = std::numeric_limits<Real>::infinity();
        constexpr Real tolerance = 1e-6;
        Vec3d dtheta = dtheta_[index_i]; // use the last step value as the initial guess
        int it = 0;
        for (; it < max_iterations; ++it)
        {
            auto [r, grad_r] = get_r_and_grad_r(dtheta);
            if (r.norm() < tolerance)
                break;
            if (grad_r.determinant() < TinyReal)
            {
                std::cout << "grad r:" << grad_r << std::endl;
                throw std::runtime_error("Error: grad_r is singular or not finite in SimoReissnerRelaxationFirstHalf_v2");
            }
            dtheta -= grad_r.inverse() * r;
        }
        if (dtheta.allFinite() == false || it == max_iterations)
            throw std::runtime_error("Error: dtheta is not finite in SimoReissnerRelaxationFirstHalf_v2");
        dtheta_[index_i] = dtheta;
        // update rotation matrix
        lambda_[index_i] = get_lambda(lambda_[index_i], dtheta);
        pseudo_b_n_[index_i] = lambda_[index_i] * b_n0_[index_i];
        pseudo_n_[index_i] = lambda_[index_i] * n0_[index_i];
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];

        Mat3d grad_r = Matd::Zero();
        Mat3d grad_dtheta = Matd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            grad_r -= (pos_[index_i] - pos_[index_j]) * gradW_ijV_j.transpose();
            grad_dtheta -= (dtheta_[index_i] - dtheta_[index_j]) * gradW_ijV_j.transpose();
        }
        // in initial local coordinate
        grad_r_[index_i] = Q0_i * grad_r * Q0_i.transpose() * B_[index_i];
        grad_dtheta_[index_i] = Q0_i * grad_dtheta * Q0_i.transpose() * B_[index_i];
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];

        // tranform back to global
        dr_ds_[index_i] = Q0_i.transpose() * grad_r_[index_i].col(xAxis);
        Vec3d dtheta_ds_ = Q0_i.transpose() * grad_dtheta_[index_i].col(xAxis);

        // compute strain measures
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Vec3d Tau_n = Tau_[index_i];
        Vec3d Kappa_n = Kappa_[index_i];
        Tau_[index_i] = Q * dr_ds_[index_i] - Vec3d::UnitX();

        Mat3d T_operator = get_T_operator(dtheta_[index_i]);
        Kappa_[index_i] += Q * T_operator.transpose() * dtheta_ds_;

        dTau_[index_i] = Tau_[index_i] - Tau_n;
        dKappa_[index_i] = Kappa_[index_i] - Kappa_n;
    }
};

class SimoReissnerRelaxationSecondHalf_v2 : public BaseBarRelaxation
{
  private:
    Real *mass_;

    Mat3d *lambda_;       // transformation matrix from initial local to current local coordinates
    Mat3d *C_N_;          // diag[EA,GA2,GA3]
    Mat3d *C_M_;          // diag[GJ,EI2,EI3]
    Mat3d *I_rho_l_;      // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;         // rho0*A, the mass per unit length
    Vec3d *moment_;       // internal moment
    Vec3d *moment_prior_; // the external torque per unit length in global coordinates
    Vec3d *dr_ds_;        // dr_ds in global coordinates

    Vec3d *Tau_;   // material tensile-shear strain
    Vec3d *Kappa_; // material bending strain
    Vec3d *f_;     // resultant force per unit length in global coordinates
    Vec3d *m_;     // resultant moment per unit length in global coordinates

    Vec3d *linear_momentum_;  // linear momentum
    Vec3d *angular_momentum_; // skew symmetric angular momentum
    Vec3d *dtheta_;           // incremental rotation

    Vec3d *dTau_;   // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_; // Kappa^n+1 - Kappa^n, used for numerical damping

    Vec3d *force_prev_;  // force at previous time step
    Vec3d *moment_prev_; // moment at previous time step
    Vec3d *force_prior_prev_;
    Vec3d *moment_prior_prev_;

    Vec3d *force_prior_discretized_;
    Vec3d *moment_prior_discretized_;

    Real numerical_damping_factor_ = 0.5;
    const Real W0_;

  public:
    explicit SimoReissnerRelaxationSecondHalf_v2(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          lambda_(particles_->registerStateVariableData<Mat3d>("TransformationFromInitialToCurrent", IdentityMatrix<Matd>::value)),
          C_N_(particles_->registerStateVariableData<Matd>("TensileConstitutiveMatrix")),
          C_M_(particles_->registerStateVariableData<Matd>("BendingConstitutiveMatrix")),
          I_rho_l_(particles_->registerStateVariableData<Matd>("InertiaMatrix")),
          A_rho_(particles_->registerStateVariableData<Real>("RhoCrossSectionArea")),
          moment_(particles_->registerStateVariableData<Vecd>("InternalMoment")),
          moment_prior_(particles_->registerStateVariableData<Vecd>("ExternalTorquePerUnitLength")),
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")),
          Tau_(particles_->registerStateVariableData<Vecd>("MaterialTensileShearStrain")),
          Kappa_(particles_->registerStateVariableData<Vecd>("MaterialBendingStrain")),
          f_(particles_->registerStateVariableData<Vecd>("ResultantForce")),
          m_(particles_->registerStateVariableData<Vecd>("ResultantMoment")),
          linear_momentum_(particles_->registerStateVariableData<Vecd>("LinearMomentum")),
          angular_momentum_(particles_->registerStateVariableData<Vecd>("AngularMomentum")),
          dtheta_(particles_->registerStateVariableData<Vecd>("IncrementalRotation")),
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")),
          force_prev_(particles_->registerStateVariableData<Vecd>("ForcePrev")),
          moment_prev_(particles_->registerStateVariableData<Vecd>("MomentPrev")),
          force_prior_prev_(particles_->registerStateVariableData<Vecd>("ForcePriorPrev")),
          moment_prior_prev_(particles_->registerStateVariableData<Vecd>("MomentPriorPrev")),
          force_prior_discretized_(particles_->registerStateVariableData<Vecd>("ForcePriorDiscretized")),
          moment_prior_discretized_(particles_->registerStateVariableData<Vecd>("MomentPriorDiscretized")),
          W0_(sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)) {}

    void initialization(size_t index_i, Real dt = 0.0)
    {
        // update force and moment at n+1
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Mat3d Q_t = Q.transpose();
        Vec3d f_m = C_N_[index_i] * (Tau_[index_i] + numerical_damping_factor_ * dTau_[index_i]);
        Vec3d m_m = C_M_[index_i] * (Kappa_[index_i] + numerical_damping_factor_ * dKappa_[index_i]);
        f_[index_i] = Q_t * f_m;
        m_[index_i] = Q_t * m_m;
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];

        // compute dn_ds and dm_ds
        Mat3d grad_f = Mat3d::Zero();
        Mat3d grad_m = Mat3d::Zero();
        Vec3d dr_ds_cross_f = dr_ds_[index_i].cross(f_[index_i]) * W0_ * Vol_[index_i];
        Vec3d f_prior = force_prior_[index_i] * W0_ * Vol_[index_i];
        Vec3d m_prior = moment_prior_[index_i] * W0_ * Vol_[index_i];
        Real weight = W0_ * Vol_[index_i]; // self-contribution
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vec3d grad_W_ijV_j = inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ij_[n] * Vol_[index_j];
            // use symmetry scheme for stress and moment
            grad_f += (f_[index_i] + f_[index_j]) * grad_W_ijV_j.transpose();
            grad_m += (m_[index_i] + m_[index_j]) * grad_W_ijV_j.transpose();
            dr_ds_cross_f += dr_ds_[index_j].cross(f_[index_j]) * inner_neighborhood.W_ij_[n] * Vol_[index_j];
            f_prior += force_prior_[index_j] * inner_neighborhood.W_ij_[n] * Vol_[index_j];
            m_prior += moment_prior_[index_j] * inner_neighborhood.W_ij_[n] * Vol_[index_j];
            weight += inner_neighborhood.W_ij_[n] * Vol_[index_j];
        }
        Vec3d df_ds = (grad_f * Q0_i.transpose()).col(xAxis);
        Vec3d dm_ds = (grad_m * Q0_i.transpose()).col(xAxis);
        dr_ds_cross_f /= weight;
        force_prior_discretized_[index_i] = f_prior / weight;
        moment_prior_discretized_[index_i] = m_prior / weight;

        // if (hourglass_control_)
        // {
        //     auto [pos_hourglass_force, theta_hourglass_torque] = compute_hourglass_control_force(index_i);
        //     df_ds += pos_hourglass_force;
        //     dm_ds += theta_hourglass_torque;
        // }
        force_[index_i] = df_ds / A_rho_[index_i] * mass_[index_i];
        moment_[index_i] = dm_ds + dr_ds_cross_f;
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        linear_momentum_[index_i] += 0.5 * (force_[index_i] + force_prior_discretized_[index_i]) / mass_[index_i] * A_rho_[index_i] * dt;
        angular_momentum_[index_i] += 0.5 * (moment_[index_i] + moment_prior_discretized_[index_i]) * dt;

        // update velocity and angular velocity
        vel_[index_i] = linear_momentum_[index_i] / A_rho_[index_i];
        Mat3d Q = getTransformationMatrix(pseudo_n_[index_i], pseudo_b_n_[index_i]);
        Mat3d i_rho_global = Q.transpose() * I_rho_l_[index_i] * Q;
        angular_vel_[index_i] = 0.5 * i_rho_global.inverse() * angular_momentum_[index_i];
    }
};
} // namespace SPH::slender_structure_dynamics