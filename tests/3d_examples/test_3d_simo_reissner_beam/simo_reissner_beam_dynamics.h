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
    return a1 * Mat3d::Identity() + a2 * S_dtheta + a3 * dtheta * dtheta.transpose();
}

class BaseSimoRelaxation : public BaseBarRelaxation
{
  protected:
    Vec3d *b_n0_;
    Vec3d *dTau_;       // Tau^n+1 - Tau^n, used for numerical damping
    Vec3d *dKappa_;     // Kappa^n+1 - Kappa^n, used for numerical damping
    Vec3d *Tau_;        // material tensile-shear strain
    Vec3d *Kappa_;      // material bending strain
    Vec3d *dr_ds_;      // dr_ds in global coordinates
    Mat3d *grad_r_;     // nabla_r in global coordinates
    Mat3d *grad_theta_; // nabla_theta in global coordinates

  public:
    explicit BaseSimoRelaxation(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation),
          b_n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialBinormalDirection", "BinormalDirection")),
          dTau_(particles_->registerStateVariableData<Vec3d>("dTau")),
          dKappa_(particles_->registerStateVariableData<Vec3d>("dKappa")),
          Tau_(particles_->registerStateVariableData<Vecd>("Tau")),
          Kappa_(particles_->registerStateVariableData<Vecd>("Kappa")),
          dr_ds_(particles_->registerStateVariableData<Vecd>("DrDs")),
          grad_r_(particles_->registerStateVariableData<Mat3d>("NablaR")),
          grad_theta_(particles_->registerStateVariableData<Mat3d>("NablaTheta")) {}
};

class SimoReissnerStressRelaxationFirstHalf : public BaseSimoRelaxation
{
  private:
    Real G0_;
    bool hourglass_control_;
    Real hourglass_control_factor_;

    Real *mass_;
    Real inv_W0_;
    Real dp_;

    Real numerical_damping_factor_ = 0.5;

    Mat3d *C_N_;               // diag[EA,GA2,GA3]
    Mat3d *C_M_;               // diag[GJ,EI2,EI3]
    Mat3d *I_rho_l_;           // rho0*diag[J,I2,I3]the inertia matrix in current local coordinates
    Real *A_rho_;              // rho0*A, the mass per unit length
    Vec3d *M_prior_;           // the external torque on a particle in global coordinates
    Vec3d *angular_acc_prior_; // the angular acceleration caused by external torque

    Vec3d *f_; // resultant force per unit length in global coordinates
    Vec3d *m_; // resultant moment per unit length in global coordinates

  public:
    explicit SimoReissnerStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                                   bool hourglass_control = false,
                                                   Real hourglass_control_factor = 0.005)
        : BaseSimoRelaxation(inner_relation),
          G0_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial()).ShearModulus()),
          hourglass_control_(hourglass_control),
          hourglass_control_factor_(hourglass_control_factor),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          inv_W0_(1.0 / sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)),
          dp_(sph_body_.getSPHAdaptation().ReferenceSpacing()),
          C_N_(particles_->registerStateVariableData<Matd>("TensileConstitutiveMatrix")),
          C_M_(particles_->registerStateVariableData<Matd>("BendingConstitutiveMatrix")),
          I_rho_l_(particles_->registerStateVariableData<Matd>("InertiaMatrix")),
          A_rho_(particles_->registerStateVariableData<Real>("RhoCrossSectionArea")),
          M_prior_(particles_->registerStateVariableData<Vecd>("ExternalTorque")),
          angular_acc_prior_(particles_->registerStateVariableData<Vecd>("PriorAngularAcceleration")),
          f_(particles_->registerStateVariableData<Vecd>("ResultantForce")),
          m_(particles_->registerStateVariableData<Vecd>("ResultantMoment")) {}

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

        angular_acc_prior_[index_i] = I_i_inverse * M_prior_[index_i] / dp_;
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

class SimoReissnerStressRelaxationSecondHalf : public BaseSimoRelaxation
{
  private:
    Eigen::Quaternion<Real> *quaternion_; // quaternion

    Vec3d *dkappa_dt_;        // spatial bending strain rate
    Vec3d *dv_ds_;            // dv_ds in global coordinates
    Mat3d *grad_angular_vel_; // nabla_w in global coordinates

  public:
    explicit SimoReissnerStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseSimoRelaxation(inner_relation),
          quaternion_(particles_->addUniqueStateVariableData<Eigen::Quaternion<Real>>("Quaternion", Eigen::Quaternion<Real>::Identity())),
          dkappa_dt_(particles_->registerStateVariableData<Vec3d>("SpatialBendingStrainRate")),
          dv_ds_(particles_->registerStateVariableData<Vec3d>("DvDs")),
          grad_angular_vel_(particles_->registerStateVariableData<Mat3d>("NablaW")) {}

    void initialization(size_t index_i, Real dt = 0.0)
    {
        pos_[index_i] += vel_[index_i] * dt;
        rotation_[index_i] += angular_vel_[index_i] * dt;

        Eigen::Quaternion<Real> omega_q(0.0, angular_vel_[index_i].x(), angular_vel_[index_i].y(), angular_vel_[index_i].z());
        auto dq_dt = Eigen::Quaternion<Real>(omega_q * quaternion_[index_i]);
        quaternion_[index_i] = quaternion_[index_i].coeffs() + 0.5 * dt * dq_dt.coeffs();
        quaternion_[index_i].normalize();

        pseudo_b_n_[index_i] = quaternion_[index_i] * b_n0_[index_i];
        pseudo_n_[index_i] = quaternion_[index_i] * n0_[index_i];
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
    Vec3d *dr0_ds_; // initial dr_ds (tangential direction)
    Vec3d *dr_ds_;  // current dr_ds (tangential direction)

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
          pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          transformation_matrix0_(particles_->getVariableDataByName<Mat3d>("TransformationMatrix")),
          B_(particles_->getVariableDataByName<Mat3d>("LinearGradientCorrectionMatrix")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        const auto &Q0_i = transformation_matrix0_[index_i];
        Mat3d grad_r = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            grad_r -= (pos0_[index_i] - pos0_[index_j]) * gradW_ijV_j.transpose();
        }
        grad_r = Q0_i * grad_r * Q0_i.transpose() * B_[index_i];
        // tranform back to global
        dr0_ds_[index_i] = Q0_i.transpose() * grad_r.col(xAxis);
        dr_ds_[index_i] = dr0_ds_[index_i];
    }
};
} // namespace SPH::slender_structure_dynamics