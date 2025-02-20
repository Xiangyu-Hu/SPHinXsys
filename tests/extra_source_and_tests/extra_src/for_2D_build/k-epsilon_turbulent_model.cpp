// #pragma once
#include "k-epsilon_turbulent_model.hpp"
namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
using TurbuIntegration2ndHalfWithWallDissipativeRiemann = ComplexInteraction<Integration2ndHalf<Inner<>, Contact<Wall>>, DissipativeRiemannSolver>;
//=================================================================================================//
BaseTurbuClosureCoeff::BaseTurbuClosureCoeff()
    : Karman_(0.41), turbu_const_E_(9.8), C_mu_(0.09), turbulent_intensity_(5.0e-2),
      sigma_k_(1.0), C_l_(1.44), C_2_(1.92), sigma_E_(1.3), turbulent_length_ratio_for_epsilon_inlet_(0.07),
      start_time_laminar_(0.0), y_star_threshold_laminar_(11.225)
{
    C_mu_25_ = pow(C_mu_, 0.25);
    C_mu_75_ = pow(C_mu_, 0.75);
}
//=================================================================================================//
Real WallFunction::get_distance_from_P_to_wall(Real y_p_constant)
{
    return y_p_constant;
}
//=================================================================================================//
Real WallFunction::get_dimensionless_velocity(Real y_star, Real time)
{
    Real dimensionless_velocity = 0.0;
    if (y_star < y_star_threshold_laminar_ && time > start_time_laminar_)
    {
        dimensionless_velocity = laminar_law_wall_function(y_star);
    }
    else
    {
        dimensionless_velocity = log_law_wall_function(y_star);
    }
    return dimensionless_velocity;
}
//=================================================================================================//
Real WallFunction::get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity)
{
    Real vel_grad_mag = log_law_velocity_gradient(vel_fric_mag, denominator_log_law);
    return vel_grad_mag;
}
//=================================================================================================//
Real WallFunction::log_law_wall_function(Real y_star)
{
    Real u_star = abs(log(turbu_const_E_ * y_star) / Karman_);
    return u_star;
}
//=================================================================================================//
Real WallFunction::laminar_law_wall_function(Real y_star)
{
    Real u_star = y_star;
    return u_star;
}
//=================================================================================================//
Real WallFunction::log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law)
{
    return vel_fric_mag * vel_fric_mag / denominator_log_law;
}
//=================================================================================================//
Real WallFunction::laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity)
{
    return vel_fric_mag * vel_fric_mag / dynamic_viscosity;
}
//=================================================================================================//
GetVelocityGradient<Inner<>>::GetVelocityGradient(BaseInnerRelation &inner_relation, Real weight_sub)
    : GetVelocityGradient<DataDelegateInner>(inner_relation),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
      weight_sub_nearwall_(weight_sub)
{
    this->particles_->addEvolvingVariable<Matd>("TurbulentVelocityGradient");
}
//=================================================================================================//
void GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] != 1)
    {
        velocity_gradient_[index_i] = Matd::Zero();
        Vecd vel_i = vel_[index_i];
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];

            Real r_ij = inner_neighborhood.r_ij_[n];
            const Vecd &e_ij = inner_neighborhood.e_ij_[n];
            if (is_near_wall_P2_[index_i] == 10 && is_near_wall_P1_[index_j] == 1)
            {
                Matd P1 = -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
                Vecd vel_diff = velocity_gradient_[index_j] * r_ij * e_ij;
                Matd P2 = -vel_diff * nablaW_ijV_j.transpose();
                velocity_gradient_[index_i] += (1 - weight_sub_nearwall_) * P1 + weight_sub_nearwall_ * P2;
            }
            else
            {
                velocity_gradient_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
            }
        }
    }
}
//=================================================================================================//
void GetVelocityGradient<Inner<>>::update(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] != 1)
    {
        velocity_gradient_[index_i] *= turbu_B_[index_i];
    }
}
//=================================================================================================//
GetVelocityGradient<Contact<Wall>>::GetVelocityGradient(BaseContactRelation &contact_relation)
    : InteractionWithWall<GetVelocityGradient>(contact_relation),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient"))
{
    this->particles_->addEvolvingVariable<Matd>("Velocity_Gradient_Wall");
}
//=================================================================================================//
void GetVelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] != 1)
    {
        Vecd vel_i = vel_[index_i];
        for (size_t k = 0; k < DataDelegateContact::contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*DataDelegateContact::contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] * contact_neighborhood.e_ij_[n];
                velocity_gradient_[index_i] += -1.0 * (vel_i)*nablaW_ijV_j.transpose();
            }
        }
    }
}
//=================================================================================================//
K_TurbulentModelInner::K_TurbulentModelInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, bool is_STL)
    : BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      dk_dt_(particles_->registerStateVariable<Real>("ChangeRateOfTKE")),
      dk_dt_without_dissipation_(particles_->registerStateVariable<Real>("ChangeRateOfTKEWithoutDissipation")),
      k_production_(particles_->registerStateVariable<Real>("K_Production")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      turbu_k_(particles_->registerStateVariable<Real>("TurbulenceKineticEnergy", Real(initial_values[0]))),
      turbu_epsilon_(particles_->registerStateVariable<Real>("TurbulentDissipation", Real(initial_values[1]))),
      turbu_mu_(particles_->registerStateVariable<Real>("TurbulentViscosity", Real(initial_values[2]))),
      turbu_strain_rate_(particles_->getVariableDataByName<Matd>("TurbulentStrainRate")),
      is_extra_viscous_dissipation_(particles_->registerStateVariable<int>("TurbulentExtraViscousDissipation", is_extr_visc_dissipa)),
      is_STL_(is_STL),
      turbu_indicator_(particles_->registerStateVariable<int>("TurbulentIndicator")),
      k_diffusion_(particles_->registerStateVariable<Real>("K_Diffusion")),
      vel_x_(particles_->registerStateVariable<Real>("Velocity_X"))
{
    particles_->addEvolvingVariable<Real>("ChangeRateOfTKE");
    particles_->addEvolvingVariable<Real>("ChangeRateOfTKEWithoutDissipation");
    particles_->addEvolvingVariable<Real>("K_Production");
    particles_->addEvolvingVariable<Real>("TurbulenceKineticEnergy");
    particles_->addEvolvingVariable<Real>("TurbulentViscosity");
    particles_->addEvolvingVariable<Real>("TurbulentDissipation");
    particles_->addEvolvingVariable<Matd>("TurbulentStrainRate");
    particles_->addEvolvingVariable<Real>("K_Diffusion");
    particles_->addEvolvingVariable<Real>("Velocity_X");
    particles_->addEvolvingVariable<int>("TurbulentIndicator");
}
//=================================================================================================//
void K_TurbulentModelInner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_mu_i = turbu_mu_[index_i];
    Real turbu_k_i = turbu_k_[index_i];

    Real mu_eff_i = turbu_mu_[index_i] / sigma_k_ + mu_;

    dk_dt_[index_i] = 0.0;
    dk_dt_without_dissipation_[index_i] = 0.0;
    Real k_derivative(0.0);
    Real k_lap(0.0);
    Matd strain_rate = Matd::Zero();
    Matd Re_stress = Matd::Zero();

    Real k_production(0.0);
    Real k_dissipation(0.0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = turbu_mu_[index_j] / sigma_k_ + mu_;
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
    }
    strain_rate = 0.5 * (velocity_gradient_[index_i].transpose() + velocity_gradient_[index_i]);

    Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();

    Matd k_production_matrix = Re_stress.array() * velocity_gradient_[index_i].array();

    if (is_near_wall_P1_[index_i] != 1)
        k_production_[index_i] = k_production_matrix.sum();

    k_production = k_production_[index_i];
    k_dissipation = turbu_epsilon_[index_i];

    dk_dt_[index_i] = k_production - k_dissipation + k_lap;
    dk_dt_without_dissipation_[index_i] = k_production + k_lap;

    k_diffusion_[index_i] = k_lap;
    vel_x_[index_i] = vel_[index_i][0];
    turbu_strain_rate_[index_i] = strain_rate;
}
//=================================================================================================//
void K_TurbulentModelInner::update(size_t index_i, Real dt)
{
    if (is_STL_)
    {
        Real denominator = 1.0 + turbu_epsilon_[index_i] * dt / turbu_k_[index_i];
        turbu_k_[index_i] += dk_dt_without_dissipation_[index_i] * dt;
        turbu_k_[index_i] /= denominator;
    }
    else
    {
        turbu_k_[index_i] += dk_dt_[index_i] * dt;
    }
}
//=================================================================================================//
E_TurbulentModelInner::E_TurbulentModelInner(BaseInnerRelation &inner_relation, bool is_STL)
    : BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      depsilon_dt_(particles_->registerStateVariable<Real>("ChangeRateOfTDR")),
      depsilon_dt_without_dissipation_(particles_->registerStateVariable<Real>("ChangeRateOfTDRWithoutDissp")),
      ep_production(particles_->registerStateVariable<Real>("Ep_Production")),
      ep_dissipation_(particles_->registerStateVariable<Real>("Ep_Dissipation_")),
      ep_diffusion_(particles_->registerStateVariable<Real>("Ep_Diffusion_")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")),
      k_production_(particles_->getVariableDataByName<Real>("K_Production")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      is_STL_(is_STL)
{
    particles_->addEvolvingVariable<Real>("ChangeRateOfTDR");
    particles_->addEvolvingVariable<Real>("Ep_Production");
    particles_->addEvolvingVariable<Real>("Ep_Dissipation_");
    particles_->addEvolvingVariable<Real>("Ep_Diffusion_");
}
//=================================================================================================//
void E_TurbulentModelInner::
    interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_epsilon_i = turbu_epsilon_[index_i];

    Real mu_eff_i = turbu_mu_[index_i] / sigma_E_ + mu_;

    depsilon_dt_[index_i] = 0.0;
    depsilon_dt_without_dissipation_[index_i] = 0.0;
    Real epsilon_production(0.0);
    Real epsilon_derivative(0.0);
    Real epsilon_lap(0.0);
    Real epsilon_dissipation(0.0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = turbu_mu_[index_j] / sigma_E_ + mu_;
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        epsilon_derivative = (turbu_epsilon_i - turbu_epsilon_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        epsilon_lap += 2.0 * mu_harmo * epsilon_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
    }

    epsilon_production = C_l_ * turbu_epsilon_i * k_production_[index_i] / turbu_k_i;
    epsilon_dissipation = C_2_ * turbu_epsilon_i * turbu_epsilon_i / turbu_k_i;

    depsilon_dt_[index_i] = epsilon_production - epsilon_dissipation + epsilon_lap;
    depsilon_dt_without_dissipation_[index_i] = epsilon_production + epsilon_lap;

    ep_production[index_i] = epsilon_production;
    ep_dissipation_[index_i] = epsilon_dissipation;
    ep_diffusion_[index_i] = epsilon_lap;
}
//=================================================================================================//
void E_TurbulentModelInner::update(size_t index_i, Real dt)
{

    if (is_near_wall_P1_[index_i] != 1)
    {
        if (is_STL_)
        {
            Real denominator = 1.0 + C_2_ * turbu_epsilon_[index_i] * dt / turbu_k_[index_i];
            turbu_epsilon_[index_i] += depsilon_dt_without_dissipation_[index_i] * dt;
            turbu_epsilon_[index_i] /= denominator;
        }
        else
        {
            turbu_epsilon_[index_i] += depsilon_dt_[index_i] * dt;
        }
    }
}
//=================================================================================================//
TKEnergyForce<Inner<>>::TKEnergyForce(BaseInnerRelation &inner_relation)
    : TKEnergyForce<Base, DataDelegateInner>(inner_relation),
      test_k_grad_rslt_(this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}
//=================================================================================================//
void TKEnergyForce<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real turbu_k_i = turbu_k_[index_i];
    Vecd force = Vecd::Zero();
    Vecd k_gradient = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];

        k_gradient += (turbu_k_i + turbu_k_[index_j]) * nablaW_ijV_j;
    }
    force = -1.0 * (2.0 / 3.0) * k_gradient * mass_[index_i];
    force_[index_i] += force;

    test_k_grad_rslt_[index_i] = k_gradient;
}
//=================================================================================================//
TKEnergyForce<Contact<>>::TKEnergyForce(BaseContactRelation &contact_relation)
    : TKEnergyForce<Base, DataDelegateContact>(contact_relation),
      test_k_grad_rslt_(this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}
//=================================================================================================//
void TKEnergyForce<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real turbu_k_i = turbu_k_[index_i];
    Vecd force = Vecd::Zero();
    Vecd k_gradient = Vecd::Zero();
    for (size_t k = 0; k < DataDelegateContact::contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*DataDelegateContact::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] * contact_neighborhood.e_ij_[n];

            k_gradient += (turbu_k_i + turbu_k_i) * nablaW_ijV_j;
        }
    }
    force = -1.0 * (2.0 / 3.0) * k_gradient * mass_[index_i];
    force_[index_i] += force;

    test_k_grad_rslt_[index_i] += k_gradient;
}
//=================================================================================================//
TurbuViscousForce<Inner<>>::TurbuViscousForce(BaseInnerRelation &inner_relation)
    : TurbuViscousForce<DataDelegateInner>(inner_relation),
      turbu_indicator_(this->particles_->template getVariableDataByName<int>("TurbulentIndicator")),
      is_extra_viscous_dissipation_(this->particles_->template getVariableDataByName<int>("TurbulentExtraViscousDissipation")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}
//=================================================================================================//
void TurbuViscousForce<Inner<>>::interaction(size_t index_i, Real dt)
{
    turbu_indicator_[index_i] = 0;

    Real mu_eff_i = turbu_mu_[index_i] + molecular_viscosity_;

    Vecd force = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real mu_eff_j = turbu_mu_[index_j] + molecular_viscosity_;
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);

        Vecd shear_stress = mu_harmo * vel_derivative;
        Vecd shear_stress_eij = shear_stress.dot(e_ij) * e_ij;

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);

        Real dissipation = rho_[index_i] * smoothing_length_ * SMIN(Real(3.0) * SMAX(u_jump, Real(0)), c0_);
        Real dissipation_judge = dissipation;

        Vecd shear_stress_eij_corrected = shear_stress_eij;
        if (mu_harmo < dissipation_judge && is_extra_viscous_dissipation_[index_i] == 1)
        {
            shear_stress_eij_corrected = ((dissipation * vel_derivative).dot(e_ij)) * e_ij;
            turbu_indicator_[index_i]++;
        }
        shear_stress = (shear_stress - shear_stress_eij) + shear_stress_eij_corrected;

        Vecd force_j = 2.0 * mass_[index_i] * shear_stress * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];

        force += force_j;
    }
    viscous_force_[index_i] = force / rho_[index_i];
}
//=================================================================================================//
TurbuViscousForce<Contact<Wall>>::TurbuViscousForce(BaseContactRelation &wall_contact_relation)
    : BaseTurbuViscousForceWithWall(wall_contact_relation),
      wall_particle_spacing_(wall_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")) {}
//=================================================================================================//
void TurbuViscousForce<Contact<Wall>>::interaction(size_t index_i, Real dt)
{

    if (this->is_near_wall_P2_[index_i] != 10)
        return;

    Real current_time = *physical_time_;
    Real turbu_k_i = this->turbu_k_[index_i];
    Real turbu_k_i_05 = pow(turbu_k_i, 0.5);
    Real rho_i = this->rho_[index_i];
    const Vecd &vel_i = this->vel_[index_i];

    Real y_p_constant_i = this->y_p_[index_i];

    Vecd force = Vecd::Zero();
    Vecd e_j_n = Vecd::Zero();
    Vecd e_j_tau = Vecd::Zero();
    Matd WSS_j_tn = Matd::Zero();
    Matd WSS_j = Matd::Zero();
    Matd Q = Matd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        Vecd *n_k = this->wall_n_[k];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            e_j_n = n_k[index_j];
            Q = getTransformationMatrix(e_j_n);

            e_j_tau[0] = e_j_n[1];
            e_j_tau[1] = e_j_n[0] * (-1.0);

            Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));

            Real y_p_j = get_distance_from_P_to_wall(y_p_constant_i);
            Real y_star_j = rho_i * C_mu_25_ * turbu_k_i_05 * y_p_j / molecular_viscosity_;
            Real u_star_j = get_dimensionless_velocity(y_star_j, current_time);
            Real fric_vel_mag_j = sqrt(C_mu_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

            Real WSS_tn_mag_j = rho_i * fric_vel_mag_j * fric_vel_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));

            WSS_j_tn(0, 0) = 0.0;
            WSS_j_tn(0, 1) = WSS_tn_mag_j;
            WSS_j_tn(1, 0) = 0.0;
            WSS_j_tn(1, 1) = 0.0;

            WSS_j = Q.transpose() * WSS_j_tn * Q;
            Vecd force_j = 2.0 * mass_[index_i] * WSS_j * e_ij * contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;

            force += force_j;
        }
    }
    viscous_force_[index_i] += force;
}
//=================================================================================================//
TurbulentEddyViscosity::
    TurbulentEddyViscosity(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")),
      wall_Y_plus_(particles_->getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(particles_->getVariableDataByName<Real>("WallYstar")),
      viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())),
      mu_(viscosity_.ReferenceViscosity()) {}
//=================================================================================================//
void TurbulentEddyViscosity::update(size_t index_i, Real dt)
{
    turbu_mu_[index_i] = rho_[index_i] * C_mu_ * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
}
//=================================================================================================//
TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_turbu_(U_max), advectionCFL_(advectionCFL),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = viscosity_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
    speed_ref_turbu_ = SMAX(viscous_speed, speed_ref_turbu_);
}
//=================================================================================================//
Real TurbulentAdvectionTimeStepSize::reduce(size_t index_i, Real dt)
{
    Real turbu_viscous_speed = (viscosity_.ReferenceViscosity() + turbu_mu_[index_i]) / fluid_.ReferenceDensity() / smoothing_length_min_;
    Real turbu_viscous_speed_squire = turbu_viscous_speed * turbu_viscous_speed;
    Real vel_n_squire = vel_[index_i].squaredNorm();
    Real vel_bigger = SMAX(turbu_viscous_speed_squire, vel_n_squire);

    return vel_bigger;
}
//=================================================================================================//
Real TurbulentAdvectionTimeStepSize::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * smoothing_length_min_ / (SMAX(speed_max, speed_ref_turbu_) + TinyReal);
}
//=================================================================================================//
InflowTurbulentCondition::InflowTurbulentCondition(BodyPartByCell &body_part, Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet)
    : BaseFlowBoundaryCondition(body_part), type_turbu_inlet_(type_turbu_inlet),
      relaxation_rate_(relaxation_rate),
      CharacteristicLength_(CharacteristicLength),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation"))
{
    TurbulentLength_ = turbulent_length_ratio_for_epsilon_inlet_ * CharacteristicLength_;
}
//=================================================================================================//
void InflowTurbulentCondition::update(size_t index_i, Real dt)
{
    Real target_in_turbu_k = getTurbulentInflowK(pos_[index_i], vel_[index_i], turbu_k_[index_i]);
    turbu_k_[index_i] += relaxation_rate_ * (target_in_turbu_k - turbu_k_[index_i]);
    Real target_in_turbu_E = getTurbulentInflowE(pos_[index_i], turbu_k_[index_i], turbu_epsilon_[index_i]);
    turbu_epsilon_[index_i] += relaxation_rate_ * (target_in_turbu_E - turbu_epsilon_[index_i]);
}
//=================================================================================================//
Real InflowTurbulentCondition::getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k)
{
    Real u = velocity[0];
    Real temp_in_turbu_k = 1.5 * pow((turbulent_intensity_ * u), 2);
    Real turbu_k_original = turbu_k;
    if (type_turbu_inlet_ == 1)
    {
        Real channel_height = CharacteristicLength_;

        Real Y = 0.0;
        if (position[1] < channel_height / 2.0)
        {
            Y = position[1];
        }
        else if (position[1] > channel_height / 2.0)
        {
            Y = channel_height - position[1];
        }

        int polynomial_order = 8;
        int num_coefficient = polynomial_order + 1;

        Real coeff[] = {
            1.159981e-02, -4.662944e-02, 2.837400e-01,
            -1.193955e+00, 3.034851e+00, -4.766077e+00,
            4.529136e+00, -2.380854e+00, 5.307586e-01};
        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }
        temp_in_turbu_k = polynomial_value;
    }
    if (position[0] < 0.0)
    {
        turbu_k_original = temp_in_turbu_k;
    }
    return turbu_k_original;
}
//=================================================================================================//
Real InflowTurbulentCondition::getTurbulentInflowE(Vecd &position, Real &turbu_k, Real &turbu_E)
{
    Real temp_in_turbu_E = C_mu_75_ * pow(turbu_k, 1.5) / TurbulentLength_;
    Real turbu_E_original = turbu_E;
    if (type_turbu_inlet_ == 1)
    {
        Real channel_height = CharacteristicLength_;

        Real Y = 0.0;
        if (position[1] < channel_height / 2.0)
        {
            Y = position[1];
        }
        else if (position[1] > channel_height / 2.0)
        {
            Y = channel_height - position[1];
        }

        int polynomial_order = 8;
        int num_coefficient = polynomial_order + 1;

        Real coeff[] = {
            1.428191e-02, -1.766636e-01, 1.153107e+00,
            -4.515606e+00, 1.103752e+01, -1.694146e+01,
            1.584534e+01, -8.241577e+00, 1.825421e+00};

        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }
        temp_in_turbu_E = polynomial_value;
    }
    if (position[0] < 0.0)
    {
        turbu_E_original = temp_in_turbu_E;
    }
    return turbu_E_original;
}
//=================================================================================================//
JudgeIsNearWall::
    JudgeIsNearWall(BaseInnerRelation &inner_relation,
                    BaseContactRelation &contact_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      distance_to_dummy_interface_(particles_->registerStateVariable<Real>("DistanceToDummyInterface")),
      distance_to_dummy_interface_up_average_(particles_->registerStateVariable<Real>("DistanceToDummyInterfaceUpAver")),
      is_near_wall_P1_(particles_->registerStateVariable<int>("IsNearWallP1")),
      is_near_wall_P2_(particles_->registerStateVariable<int>("IsNearWallP2")),
      index_nearest_(particles_->registerStateVariable<int>("NearestIndex")),
      e_nearest_tau_(particles_->registerStateVariable<Vecd>("WallNearestTangentialUnitVector")),
      e_nearest_normal_(particles_->registerStateVariable<Vecd>("WallNearestNormalUnitVector")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      dimension_(2),
      fluid_particle_spacing_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
      wall_particle_spacing_(contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
      distance_from_wall_(particles_->getVariableDataByName<Vecd>("DistanceFromWall"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }

    particles_->addEvolvingVariable<Real>("DistanceToDummyInterfaceUpAver");
    particles_->addEvolvingVariable<Real>("DistanceToDummyInterface");
    particles_->addEvolvingVariable<int>("NearestIndex");
    particles_->addEvolvingVariable<int>("IsNearWallP1");
    particles_->addEvolvingVariable<int>("IsNearWallP2");
    particles_->addEvolvingVariable<Vecd>("WallNearestTangentialUnitVector");
    particles_->addEvolvingVariable<Vecd>("WallNearestNormalUnitVector");
};
//=================================================================================================//
void JudgeIsNearWall::interaction(size_t index_i, Real dt)
{
    is_near_wall_P2_[index_i] = 0;
    index_nearest_[index_i] = 0;
    distance_to_dummy_interface_[index_i] = 0.0;
    distance_to_dummy_interface_up_average_[index_i] = 0.0;
    e_nearest_normal_[index_i] = Vecd::Zero();
    e_nearest_tau_[index_i] = Vecd::Zero();

    int id_nearest_j = 0;
    Real r_dummy_normal = 0.0;
    Real r_dummy_normal_j = 0.0;
    Real r_min = 1.0e3;
    Vecd e_i_nearest_tau = Vecd::Zero();
    Vecd e_i_nearest_n = Vecd::Zero();
    Real ttl_weight(0);
    Real r_dmy_itfc_n_sum(0);
    int count_average(0);
    int is_near_contact = 0;

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = (contact_Vol_[k]);
        Vecd *n_k = (contact_n_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        if (contact_neighborhood.current_size_ != 0)
            is_near_contact++;

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Vecd &n_k_j = n_k[index_j];

            r_dummy_normal_j = abs(n_k_j.dot(r_ij * e_ij)) - 0.5 * wall_particle_spacing_;

            if (r_ij < r_min && r_dummy_normal_j > 0.0 + TinyReal)
            {
                r_min = r_ij;

                r_dummy_normal = r_dummy_normal_j;
                e_i_nearest_n = n_k[index_j];
                id_nearest_j = index_j;
            }

            if (r_dummy_normal_j - r_dummy_normal > (-1.0 * TinyReal))
            {
                count_average++;

                r_dmy_itfc_n_sum += weight_j * r_dummy_normal_j;
                ttl_weight += weight_j;
            }
        }
    }

    if (is_near_contact > 0)
    {
        is_near_wall_P2_[index_i] = 10;

        if (dimension_ == 2)
        {
            e_i_nearest_tau[0] = e_i_nearest_n[1];
            e_i_nearest_tau[1] = e_i_nearest_n[0] * (-1.0);
        }

        distance_to_dummy_interface_up_average_[index_i] = r_dmy_itfc_n_sum / ttl_weight;

        index_nearest_[index_i] = id_nearest_j;
        e_nearest_normal_[index_i] = e_i_nearest_n;
        e_nearest_tau_[index_i] = e_i_nearest_tau;
        distance_to_dummy_interface_[index_i] = r_dummy_normal;
    }
}
//=================================================================================================//
void JudgeIsNearWall::update(size_t index_i, Real dt)
{
    is_near_wall_P1_[index_i] = 0;
    if (is_near_wall_P2_[index_i] == 10)
    {

        Real distance = distance_to_dummy_interface_[index_i];

        if (distance < 1.0 * fluid_particle_spacing_)
            is_near_wall_P1_[index_i] = 1;
    }
}
//=================================================================================================//
StandardWallFunctionCorrection::
    StandardWallFunctionCorrection(BaseInnerRelation &inner_relation,
                                   BaseContactRelation &contact_relation, Real y_p_constant)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      y_p_(particles_->registerStateVariable<Real>("Y_P", y_p_constant)),
      wall_Y_plus_(particles_->registerStateVariable<Real>("WallYplus")),
      wall_Y_star_(particles_->registerStateVariable<Real>("WallYstar")),
      velo_tan_(particles_->registerStateVariable<Real>("TangentialVelocity")),
      velo_friction_(particles_->registerStateVariable<Vecd>("FrictionVelocity")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")), rho_(particles_->getVariableDataByName<Real>("Density")),
      viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())),
      molecular_viscosity_(viscosity_.ReferenceViscosity()),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(particles_->getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      k_production_(particles_->getVariableDataByName<Real>("K_Production")),
      distance_to_dummy_interface_(particles_->getVariableDataByName<Real>("DistanceToDummyInterface")),
      distance_to_dummy_interface_up_average_(particles_->getVariableDataByName<Real>("DistanceToDummyInterfaceUpAver")),
      index_nearest(particles_->getVariableDataByName<int>("NearestIndex")),
      e_nearest_tau_(particles_->getVariableDataByName<Vecd>("WallNearestTangentialUnitVector")),
      e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")),
      physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }

    particles_->addEvolvingVariable<Real>("Y_P");
    particles_->addEvolvingVariable<Real>("WallYplus");
    particles_->addEvolvingVariable<Real>("WallYstar");
    particles_->addEvolvingVariable<Real>("TangentialVelocity");
    particles_->addEvolvingVariable<Vecd>("FrictionVelocity");
};
//=================================================================================================//
void StandardWallFunctionCorrection::interaction(size_t index_i, Real dt)
{
    velo_tan_[index_i] = 0.0;
    velo_friction_[index_i] = Vecd::Zero();
    wall_Y_plus_[index_i] = 0.0;
    wall_Y_star_[index_i] = 0.0;
    Real current_time = *physical_time_;

    if (is_near_wall_P2_[index_i] == 10)
    {
        Real y_p_constant_i = y_p_[index_i];

        Real turbu_k_i_05 = pow(turbu_k_[index_i], 0.5);
        Real turbu_k_i_15 = pow(turbu_k_[index_i], 1.5);

        Vecd e_i_nearest_tau = e_nearest_tau_[index_i];
        Vecd e_i_nearest_n = e_nearest_normal_[index_i];
        const Vecd &vel_i = vel_[index_i];
        Real rho_i = rho_[index_i];
        Real nu_i = molecular_viscosity_ / rho_i;

        wall_Y_star_[index_i] = y_p_constant_i * C_mu_25_ * turbu_k_i_05 / nu_i;

        Real velo_fric_mag = 0.0;
        Real velo_tan_mag = 0.0;

        velo_tan_mag = abs(e_i_nearest_tau.dot(vel_i));
        velo_tan_[index_i] = velo_tan_mag;

        Real u_star = get_dimensionless_velocity(wall_Y_star_[index_i], current_time);
        velo_fric_mag = sqrt(C_mu_25_ * turbu_k_i_05 * velo_tan_mag / u_star);
        velo_friction_[index_i] = velo_fric_mag * e_i_nearest_tau;
        if (vel_i.dot(velo_friction_[index_i]) < 0.0)
            velo_friction_[index_i] = -1.0 * velo_friction_[index_i];

        wall_Y_plus_[index_i] = y_p_constant_i * velo_fric_mag / nu_i;

        if (is_near_wall_P1_[index_i] == 1)
        {
            Matd vel_grad_i_tn = Matd::Zero();
            Matd Q = Matd::Zero();
            Real total_weight = 0.0;

            Real epsilon_p_weighted_sum = 0.0;
            Real dudn_p_weighted_sum = 0.0;
            Real G_k_p_weighted_sum = 0.0;

            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Real *Vol_k = (contact_Vol_[k]);
                Vecd *n_k = (contact_n_[k]);
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    Real epsilon_p_j = 0.0;
                    Real dudn_p_j = 0.0;
                    Real G_k_p_j = 0.0;

                    Real y_p_j = 0.0;

                    Vecd e_j_tau = Vecd::Zero();

                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd e_j_n = n_k[index_j];

                    e_j_tau[0] = e_j_n[1];
                    e_j_tau[1] = e_j_n[0] * (-1.0);

                    y_p_j = get_distance_from_P_to_wall(y_p_constant_i);

                    Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                    total_weight += weight_j;

                    Real denominator_log_law_j = C_mu_25_ * turbu_k_i_05 * Karman_ * y_p_j;

                    Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
                    Real y_star_j = C_mu_25_ * turbu_k_i_05 * y_p_j / nu_i;
                    Real u_star_j = get_dimensionless_velocity(y_star_j, current_time);
                    Real fric_vel_mag_j = sqrt(C_mu_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

                    Real dudn_p_mag_j = get_near_wall_velocity_gradient_magnitude(y_star_j, fric_vel_mag_j, denominator_log_law_j, nu_i);
                    dudn_p_j = dudn_p_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));
                    dudn_p_weighted_sum += weight_j * dudn_p_j;

                    if (y_star_j < y_star_threshold_laminar_ && current_time > start_time_laminar_)
                    {
                        epsilon_p_j = 2.0 * turbu_k_[index_i] * nu_i / (y_p_j * y_p_j);
                        G_k_p_j = 0.0;
                    }
                    else
                    {
                        epsilon_p_j = C_mu_75_ * turbu_k_i_15 / (Karman_ * y_p_j);
                        G_k_p_j = rho_i * fric_vel_mag_j * fric_vel_mag_j * dudn_p_mag_j;
                    }

                    epsilon_p_weighted_sum += weight_j * epsilon_p_j;
                    G_k_p_weighted_sum += weight_j * G_k_p_j;
                }
            }
            turbu_epsilon_[index_i] = epsilon_p_weighted_sum / total_weight;

            vel_grad_i_tn(0, 0) = 0.0;
            vel_grad_i_tn(0, 1) = dudn_p_weighted_sum / total_weight;
            vel_grad_i_tn(1, 0) = 0.0;
            vel_grad_i_tn(1, 1) = 0.0;

            Q = getTransformationMatrix(e_i_nearest_n);

            velocity_gradient_[index_i] = Q.transpose() * vel_grad_i_tn * Q;

            k_production_[index_i] = G_k_p_weighted_sum / total_weight;
        }
    }
}
//=================================================================================================//
ConstrainNormalVelocityInRegionP::
    ConstrainNormalVelocityInRegionP(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")) {}
//=================================================================================================//
void ConstrainNormalVelocityInRegionP::update(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] == 1)
    {
        vel_[index_i] = vel_[index_i] - (vel_[index_i].dot(e_nearest_normal_[index_i])) * e_nearest_normal_[index_i];
    }
}
//=================================================================================================//
void TurbulentLinearGradientCorrectionMatrix<Inner<>>::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity();

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ij.transpose();
    }
    turbu_B_[index_i] = local_configuration;
}
//=================================================================================================//
void TurbulentLinearGradientCorrectionMatrix<Inner<>>::update(size_t index_i, Real dt)
{
    Real det_sqr = SMAX(turbu_alpha_ - turbu_B_[index_i].determinant(), Real(0));
    Matd inverse = turbu_B_[index_i].inverse();
    Real weight1_ = turbu_B_[index_i].determinant() / (turbu_B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (turbu_B_[index_i].determinant() + det_sqr);
    turbu_B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
}
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//