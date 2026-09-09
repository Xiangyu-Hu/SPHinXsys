
#include "udf_k-omega_turbulent_model.hpp"
namespace SPH
{

namespace fluid_dynamics
{
namespace udf
{

kOmega_BaseTurbuClosureCoeff::kOmega_BaseTurbuClosureCoeff()
    : std_kw_beta_star_(0.09), std_kw_sigma_star_(0.6),
      std_kw_alpha_(0.52), std_kw_sigma_(0.5), std_kw_f_beta_(1.0), std_kw_beta_0_(0.0708),
      std_kw_sigma_do_(0.125), std_kw_C_lim_(0.875), std_kw_beta_i_(0.075),
      turbulent_intensity_for_k_inlet_(0.05), turbulent_length_ratio_for_omega_inlet_(0.07),
      C_mu_75_for_omega_inlet_(pow(0.09, 0.75))
{
    std_kw_beta_ = std_kw_beta_0_ * std_kw_f_beta_;
    std_kw_beta_star_25_ = pow(std_kw_beta_star_, 0.25);
    std_kw_beta_star_5_ = pow(std_kw_beta_star_, 0.5);
}

kOmega_GetVelocityGradient<Inner<>>::kOmega_GetVelocityGradient(BaseInnerRelation &inner_relation)
    : kOmega_GetVelocityGradient<DataDelegateInner>(inner_relation),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
      turbu_strain_rate_(particles_->registerStateVariableData<Matd>("TurbulentStrainRate")),
      turbu_strain_rate_magnitude_(particles_->registerStateVariableData<Real>("TurbulentStrainRateMagnitude")),
      turbu_strain_rate_traceless_magnitude_(particles_->registerStateVariableData<Real>("TurbulentStrainRateTracelessMagnitude")){}

void kOmega_GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
{
    velocity_gradient_[index_i] = Matd::Zero();
    Vecd vel_i = vel_[index_i];
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
        velocity_gradient_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
    }
}

void kOmega_GetVelocityGradient<Inner<>>::update(size_t index_i, Real dt)
{
    velocity_gradient_[index_i] *= B_[index_i];
    turbu_strain_rate_[index_i] = 0.5 * (velocity_gradient_[index_i].transpose() + velocity_gradient_[index_i]);
    Real strain_rate_trace = turbu_strain_rate_[index_i].trace();
    Matd strain_rate_traceless = turbu_strain_rate_[index_i] - (1.0 / Dimensions) * strain_rate_trace * Matd::Identity();

    Real strain_rate_squire = (turbu_strain_rate_[index_i].array() * turbu_strain_rate_[index_i].array()).sum();
    turbu_strain_rate_magnitude_[index_i] = sqrt(2.0 * strain_rate_squire);
    Real strain_rate_traceless_squire = (strain_rate_traceless.array() * strain_rate_traceless.array()).sum();
    turbu_strain_rate_traceless_magnitude_[index_i] = sqrt(2.0 * strain_rate_traceless_squire);
}

kOmega_GetVelocityGradient<Contact<Wall>>::kOmega_GetVelocityGradient(BaseContactRelation &contact_relation)
    : InteractionWithWall<kOmega_GetVelocityGradient>(contact_relation),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")) {}

void kOmega_GetVelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Matd vel_grad = Matd::Zero();

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Real *Vol_k = wall_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            vel_grad += -1.0 * 2.0 * (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
        }
    }
    velocity_gradient_[index_i] += vel_grad;
}

kOmegaTurbulentEddyViscosity::
    kOmegaTurbulentEddyViscosity(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      wall_Y_plus_(particles_->getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(particles_->getVariableDataByName<Real>("WallYstar")),
      turbu_strain_rate_traceless_magnitude_(particles_->getVariableDataByName<Real>("TurbulentStrainRateTracelessMagnitude")),
      viscosity_(sph_body_->getMaterialProperty<Viscosity>()),
      mu_(viscosity_.ReferenceViscosity()) {}

void kOmegaTurbulentEddyViscosity::update(size_t index_i, Real dt)
{
    Real limited_omega = std_kw_C_lim_ * turbu_strain_rate_traceless_magnitude_[index_i] / sqrt(std_kw_beta_star_);
    Real turbu_omega_tilde_ = SMAX(turbu_omega_[index_i], limited_omega);
    turbu_mu_[index_i] = rho_[index_i] * turbu_k_[index_i] / turbu_omega_tilde_;
}

kOmega_WallFunctionCorrection::
    kOmega_WallFunctionCorrection(BaseInnerRelation &inner_relation,
                                  BaseContactRelation &contact_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      y_p_(particles_->getVariableDataByName<Real>("Y_P")),
      wall_Y_plus_(particles_->registerStateVariableData<Real>("WallYplus")),
      wall_Y_star_(particles_->registerStateVariableData<Real>("WallYstar")),
      velo_tan_(particles_->registerStateVariableData<Real>("TangentialVelocity")),
      velo_friction_(particles_->registerStateVariableData<Vecd>("FrictionVelocity")),
      wall_shear_stress_(particles_->registerStateVariableData<Real>("WallShearStress")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")), rho_(particles_->getVariableDataByName<Real>("Density")),
      viscosity_(sph_body_->getMaterialProperty<Viscosity>()),
      molecular_viscosity_(viscosity_.ReferenceViscosity()),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
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
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")),
      is_blended_(particles_->getVariableDataByName<int>("TurbulentWallTreatmentType")),
      turbu_strain_rate_magnitude_(particles_->getVariableDataByName<Real>("TurbulentStrainRateMagnitude")),
      laminar_fraction_for_blend_(particles_->registerStateVariableData<Real>("LaminarFractionForBlend"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
};

void kOmega_WallFunctionCorrection::interaction(size_t index_i, Real dt)
{
    velo_tan_[index_i] = 0.0;
    Real vel_fric_mag_previous = velo_friction_[index_i].norm();
    velo_friction_[index_i] = Vecd::Zero();
    wall_Y_plus_[index_i] = 0.0;
    wall_Y_star_[index_i] = 0.0;
    Real current_time = *physical_time_;

    wall_shear_stress_[index_i] = 0.0;
    laminar_fraction_for_blend_[index_i] = 0.0;
    if (is_near_wall_P2_[index_i] == 10)
    {
        Real y_p_constant_i = y_p_[index_i];
        Real turbu_k_i = turbu_k_[index_i];
        Real turbu_k_i_05 = pow(turbu_k_i, 0.5);
        Vecd e_i_nearest_tau = e_nearest_tau_[index_i];
        Vecd e_i_nearest_n = e_nearest_normal_[index_i];
        const Vecd &vel_i = vel_[index_i];
        Real rho_i = rho_[index_i];
        Real nu_i = molecular_viscosity_ / rho_i;
        wall_Y_star_[index_i] = y_p_constant_i * C_mu_wf_25_ * turbu_k_i_05 / nu_i;
        Real velo_fric_mag = 0.0;
        Real velo_tan_mag = 0.0;

        velo_tan_mag = abs(e_i_nearest_tau.dot(vel_i));
        velo_tan_[index_i] = velo_tan_mag;

        Real u_star_previous = velo_tan_mag / vel_fric_mag_previous;
        if ((u_star_previous > 100.0 || u_star_previous <= TinyReal) && current_time > start_time_laminar_)
        {
            u_star_previous = wall_Y_star_[index_i] + 10.0 * TinyReal;
        }

        if (wall_Y_star_[index_i] != static_cast<Real>(wall_Y_star_[index_i]))
        {
            std::cout << "y* is not a real value, please check" << std::endl;
            std::cin.get();
        }

        Real u_star = get_dimensionless_velocity(wall_Y_star_[index_i], current_time, u_star_previous, is_blended_[index_i]);
        velo_fric_mag = sqrt(C_mu_wf_25_ * turbu_k_i_05 * velo_tan_mag / u_star);

        if (velo_fric_mag != static_cast<Real>(velo_fric_mag))
        {
            std::cout << "friction velocity is not a real, please check" << std::endl;
            std::cout << "velo_fric=" << velo_fric_mag << std::endl
                      << "velo_tan_mag=" << velo_tan_mag << std::endl;
            std::cout << "turbu_k_=" << pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "sum=" << (Karman_ * velo_tan_mag * C_mu_wf_25_ * pow(turbu_k_[index_i], 0.5) / log(turbu_const_E_ * C_mu_wf_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_)) << std::endl;
            std::cout << "numerator=" << Karman_ * velo_tan_mag * C_mu_wf_25_ * pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "denominator=" << log(turbu_const_E_ * C_mu_wf_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_) << std::endl;
            Real temp = C_mu_wf_25_ * pow(turbu_k_[index_i], 0.5) * velo_tan_mag / u_star;

            std::cout << "temp =" << temp << std::endl;

            std::cout << "pow(turbu_k_[index_i], 0.5) =" << pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "velo_tan_mag / u_star =" << velo_tan_mag / u_star << std::endl;
            std::cout << "velo_tan_mag =" << velo_tan_mag << std::endl;
            std::cout << " u_star =" << u_star << std::endl;
            std::cin.get();
        }

        velo_friction_[index_i] = velo_fric_mag * e_i_nearest_tau;
        if (vel_i.dot(velo_friction_[index_i]) < 0.0)
            velo_friction_[index_i] = -1.0 * velo_friction_[index_i];
        wall_Y_plus_[index_i] = y_p_constant_i * velo_fric_mag / nu_i;
        if (is_near_wall_P1_[index_i] == 1)
        {
            wall_shear_stress_[index_i] = rho_i * velo_fric_mag * velo_fric_mag;

            Matd vel_grad_i_tn = Matd::Zero();
            Matd Q = Matd::Zero();
            Real total_weight = 0.0;

            Real dudn_p_weighted_sum = 0.0;
            Real G_k_p_weighted_sum = 0.0;
            Real omega_p_weighted_sum = 0.0;

            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Real *Vol_k = (contact_Vol_[k]);
                Vecd *n_k = (contact_n_[k]);
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    Real dudn_p_j = 0.0;
                    Real G_k_p_j = 0.0;
                    Real omega_p_j = 0.0;

                    Real y_p_j = 0.0;

                    Vecd e_j_tau = Vecd::Zero();

                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd e_j_n = n_k[index_j];
                    e_j_tau[0] = e_j_n[1];
                    e_j_tau[1] = e_j_n[0] * (-1.0);

                    y_p_j = get_distance_from_P_to_wall(y_p_constant_i);

                    Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                    total_weight += weight_j;

                    Real denominator_log_law_j = C_mu_wf_25_ * turbu_k_i_05 * Karman_ * y_p_j;

                    Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
                    Real y_star_j = C_mu_wf_25_ * turbu_k_i_05 * y_p_j / nu_i;
                    Real u_star_j = get_dimensionless_velocity(y_star_j, current_time, u_star_previous, is_blended_[index_i]);
                    Real fric_vel_mag_j = sqrt(C_mu_wf_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

                    Real dudn_p_mag_j = get_near_wall_velocity_gradient_magnitude(y_star_j, fric_vel_mag_j, denominator_log_law_j, nu_i);
                    dudn_p_j = dudn_p_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));
                    dudn_p_weighted_sum += weight_j * dudn_p_j;

                    if (is_blended_[index_i])
                    {
                        Real local_Re = y_p_j * turbu_k_i_05 / nu_i;
                        Real lam_frac = std::exp(-1.0 * local_Re / 11.0);
                        Real turbu_frac = 1.0 - lam_frac;
                        Real u_tau_temp_square_laminar = nu_i * vel_i_tau_mag / y_p_j;
                        Real u_tau_temp_square_turbulent = std_kw_beta_star_5_ * turbu_k_i;
                        Real u_tau_temp = sqrt(lam_frac * u_tau_temp_square_laminar + turbu_frac * u_tau_temp_square_turbulent);

                        Real omega_p_j_lam = 6.0 * nu_i / (std_kw_beta_i_ * y_p_j * y_p_j);

                        Real omega_p_j_turbu = u_tau_temp / (std_kw_beta_star_5_ * Karman_ * y_p_j);
                        Real G_k_p_j_turbu = (u_tau_temp * vel_i_tau_mag / u_star_j) * (u_tau_temp * vel_i_tau_mag / u_star_j) / (nu_i * Karman_ * y_star_j);

                        Real G_lam_k_p = turbu_mu_[index_i] / rho_i * turbu_strain_rate_magnitude_[index_i] * turbu_strain_rate_magnitude_[index_i];

                        omega_p_j = lam_frac * omega_p_j_lam + turbu_frac * omega_p_j_turbu;
                        G_k_p_j = lam_frac * G_lam_k_p + turbu_frac * G_k_p_j_turbu;
                    }
                    else
                    {
                        if (y_star_j < y_star_threshold_laminar_ && current_time > start_time_laminar_)
                        {
                            G_k_p_j = 0.0;
                            omega_p_j = 6.0 * nu_i / (std_kw_beta_i_ * y_p_j * y_p_j);
                        }
                        else
                        {
                            Real u_tau_temp = C_mu_wf_25_ * turbu_k_i_05;
                            G_k_p_j = (u_tau_temp * vel_i_tau_mag / u_star_j) * (u_tau_temp * vel_i_tau_mag / u_star_j) / (nu_i * Karman_ * y_star_j);
                            omega_p_j = u_tau_temp / (std_kw_beta_star_5_ * Karman_ * y_p_j);
                        }
                    }
                    G_k_p_weighted_sum += weight_j * G_k_p_j;
                    omega_p_weighted_sum += weight_j * omega_p_j;
                }
            }
            turbu_omega_[index_i] = omega_p_weighted_sum / total_weight;
            vel_grad_i_tn(0, 0) = 0.0;
            vel_grad_i_tn(0, 1) = dudn_p_weighted_sum / total_weight;
            vel_grad_i_tn(1, 0) = 0.0;
            vel_grad_i_tn(1, 1) = 0.0;
            Q = getTransformationMatrix(e_i_nearest_n);
            k_production_[index_i] = G_k_p_weighted_sum / total_weight;
            laminar_fraction_for_blend_[index_i] = std::exp(-1.0 * y_p_constant_i * turbu_k_i_05 / nu_i / 11.0);
        }
    }
}

kOmega_kTransportEquationInner::kOmega_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, int is_blended)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      dk_dt_(particles_->registerStateVariableData<Real>("ChangeRateOfTKE")),
      dk_dt_without_dissipation_(particles_->registerStateVariableData<Real>("ChangeRateOfTKEWithoutDissipation")),
      k_production_(particles_->registerStateVariableData<Real>("K_Production")),
      turbu_k_(particles_->registerStateVariableData<Real>("TurbulenceKineticEnergy", Real(initial_values[0]))),
      turbu_omega_(particles_->registerStateVariableData<Real>("TurbulentSpecificDissipation", Real(initial_values[1]))),
      turbu_mu_(particles_->registerStateVariableData<Real>("TurbulentViscosity", Real(initial_values[2]))),
      is_extra_viscous_dissipation_(particles_->registerStateVariableData<int>("TurbulentExtraViscousDissipation", is_extr_visc_dissipa)),
      is_blended_(particles_->registerStateVariableData<int>("TurbulentWallTreatmentType", is_blended)),
      turbu_indicator_(particles_->registerStateVariableData<int>("TurbulentIndicator")),
      k_diffusion_(particles_->registerStateVariableData<Real>("K_Diffusion")),
      turbu_strain_rate_(particles_->getVariableDataByName<Matd>("TurbulentStrainRate")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient"))
{
    particles_->addEvolvingVariable<Real>("TurbulenceKineticEnergy");
    particles_->addEvolvingVariable<Real>("TurbulentSpecificDissipation");
}

void kOmega_kTransportEquationInner::update(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_mu_i = turbu_mu_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];
    Real k_diffusion = k_diffusion_[index_i];
    Matd strain_rate = turbu_strain_rate_[index_i];

    Real k_dissipation = std_kw_beta_star_ * turbu_k_i * turbu_omega_i;

    dk_dt_[index_i] = 0.0;
    dk_dt_without_dissipation_[index_i] = 0.0;

    Real k_production(0.0);
    Matd Re_stress = Matd::Zero();

    Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
    Matd k_production_matrix = Re_stress.array() * velocity_gradient_[index_i].array();
    if (is_near_wall_P1_[index_i] != 1)
        k_production_[index_i] = k_production_matrix.sum();

    k_production = k_production_[index_i];

    dk_dt_[index_i] = k_production - k_dissipation + k_diffusion;
    dk_dt_without_dissipation_[index_i] = k_production + k_diffusion;

    turbu_k_[index_i] += dk_dt_[index_i] * dt;
}

kOmega_TKE_Diffusion::kOmega_TKE_Diffusion(BaseInnerRelation &inner_relation)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      k_diffusion_(particles_->getVariableDataByName<Real>("K_Diffusion")) {}

void kOmega_TKE_Diffusion::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];

    Real mu_eff_i = mu_ + std_kw_sigma_star_ * rho_i * turbu_k_i / turbu_omega_i;

    Real k_derivative(0.0);
    Real k_lap(0.0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = mu_ + std_kw_sigma_star_ * rho_i * turbu_k_[index_j] / turbu_omega_[index_j];
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
    }
    k_diffusion_[index_i] = k_lap / rho_i;
}

kOmega_omegaTransportEquationInner::kOmega_omegaTransportEquationInner(BaseInnerRelation &inner_relation)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      domega_dt_(particles_->registerStateVariableData<Real>("ChangeRateOfTDR")),
      domega_dt_without_dissipation_(particles_->registerStateVariableData<Real>("ChangeRateOfTDRWithoutDissp")),
      omega_production_(particles_->registerStateVariableData<Real>("omega_Production")),
      omega_dissipation_(particles_->registerStateVariableData<Real>("omega_Dissipation")),
      omega_diffusion_(particles_->registerStateVariableData<Real>("omega_Diffusion")),
      gradient_dot_k_omega_(particles_->registerStateVariableData<Real>("GradientDotKOmega")),
      omega_cross_diffusion_(particles_->registerStateVariableData<Real>("omega_Cross_Diffusion")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      k_production_(particles_->getVariableDataByName<Real>("K_Production")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")){}

void kOmega_omegaTransportEquationInner::update(size_t index_i, Real dt)
{
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];
    Real omega_diffusion = omega_diffusion_[index_i];
    Real grad_dot_k_omega = gradient_dot_k_omega_[index_i];

    domega_dt_[index_i] = 0.0;
    domega_dt_without_dissipation_[index_i] = 0.0;
    Real omega_production(0.0);
    Real omega_dissipation(0.0);
    Real omega_cross_diffusion(0.0);

    omega_production = std_kw_alpha_ * turbu_omega_i * k_production_[index_i] / turbu_k_i;
    omega_dissipation = std_kw_beta_ * turbu_omega_i * turbu_omega_i;

    std_kw_sigma_d_ = grad_dot_k_omega > 0.0 ? std_kw_sigma_do_ : 0.0;

    omega_cross_diffusion = std_kw_sigma_d_ / turbu_omega_i * grad_dot_k_omega;

    domega_dt_[index_i] = omega_production - omega_dissipation + omega_diffusion + omega_cross_diffusion;
    domega_dt_without_dissipation_[index_i] = omega_production + omega_diffusion + omega_cross_diffusion;
    if (is_near_wall_P1_[index_i] != 1)
    {
        turbu_omega_[index_i] += domega_dt_[index_i] * dt;
    }
    omega_production_[index_i] = omega_production;
    omega_dissipation_[index_i] = omega_dissipation;
    omega_cross_diffusion_[index_i] = omega_cross_diffusion;
}

kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner::kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner(BaseInnerRelation &inner_relation)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      gradient_dot_k_omega_(particles_->getVariableDataByName<Real>("GradientDotKOmega")),
      omega_diffusion_(particles_->getVariableDataByName<Real>("omega_Diffusion")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}

void kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];

    Real mu_eff_i = mu_ + std_kw_sigma_ * rho_i * turbu_k_i / turbu_omega_i;

    Real omega_derivative(0.0);
    Real omega_lap(0.0);
    Vecd k_gradient = Vecd::Zero();
    Vecd omega_gradient = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = mu_ + std_kw_sigma_ * rho_i * turbu_k_[index_j] / turbu_omega_[index_j];
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        omega_derivative = (turbu_omega_i - turbu_omega_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        omega_lap += 2.0 * mu_harmo * omega_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
        k_gradient += -1.0 * (turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
        omega_gradient += -1.0 * (turbu_omega_i - turbu_omega_[index_j]) * nablaW_ijV_j;
    }
    k_gradient = B_[index_i] * k_gradient;
    omega_gradient = B_[index_i] * omega_gradient;

    omega_diffusion_[index_i] = omega_lap / rho_i;
    gradient_dot_k_omega_[index_i] = k_gradient.dot(omega_gradient);
}

kOmega_InflowTurbulentCondition::kOmega_InflowTurbulentCondition(BodyPartByCell &body_part, Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet_omega, int type_turbu_inlet_k)
    : BaseFlowBoundaryCondition(body_part), type_turbu_inlet_omega_(type_turbu_inlet_omega), type_turbu_inlet_k_(type_turbu_inlet_k),
      relaxation_rate_(relaxation_rate),
      CharacteristicLength_(CharacteristicLength),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation"))
{
    TurbulentLength_ = turbulent_length_ratio_for_omega_inlet_ * CharacteristicLength_;
}

void kOmega_InflowTurbulentCondition::update(size_t index_i, Real dt)
{
    Real target_inflow_turbu_k = getTurbulentInflowK(pos_[index_i], vel_[index_i], turbu_k_[index_i]);
    turbu_k_[index_i] += relaxation_rate_ * (target_inflow_turbu_k - turbu_k_[index_i]);

    if (type_turbu_inlet_omega_ == 1 || type_turbu_inlet_omega_ == 0)
    {
        Real target_inflow_temp_turbu_E = getTurbulentInflowTemporaryEpsilon(pos_[index_i], turbu_k_[index_i], 0.0);
        Real inflow_turbu_omega_temp = target_inflow_temp_turbu_E / std_kw_beta_star_ / target_inflow_turbu_k;
        Real target_inflow_turbu_omega = (pos_[index_i][0] < 0.0) ? inflow_turbu_omega_temp : turbu_omega_[index_i];
        turbu_omega_[index_i] += relaxation_rate_ * (target_inflow_turbu_omega - turbu_omega_[index_i]);
    }
    else if (type_turbu_inlet_omega_ == 2)
    {
        Real target_inflow_turbu_omega = getTurbulentInflowOmega(pos_[index_i], vel_[index_i], turbu_k_[index_i]);;
        turbu_omega_[index_i] += relaxation_rate_ * (target_inflow_turbu_omega - turbu_omega_[index_i]);
    }
}

Real kOmega_InflowTurbulentCondition::getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k)
{
    Real u = velocity[0];
    Real temp_in_turbu_k = 1.5 * pow((turbulent_intensity_for_k_inlet_ * u), 2);
    Real turbu_k_original = turbu_k;
    if (type_turbu_inlet_k_ == 1)
    {
        Real channel_height = CharacteristicLength_;
        Real Y = (position[1] < channel_height / 2.0) ? position[1] : channel_height - position[1];

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

        if (Y > channel_height / 2.0 || Y < 0.0)
        {
            std::cout << "position[1]=" << position[1] << std::endl;
            std::cout << "Y=" << Y << std::endl;
            std::cout << "polynomial_value=" << polynomial_value << std::endl;
            std::cout << "Stop" << std::endl;
            std::cout << "=================" << std::endl;
            std::cin.get();
        }

        temp_in_turbu_k = polynomial_value;
    }
    if (type_turbu_inlet_k_ == 2)
    {
        const Real channel_height = CharacteristicLength_;
        Real Y = (position[1] < channel_height / 2.0) ? position[1] : channel_height - position[1];
        const Real y1 = 0.2;

        Real polynomial_value = 0.0;
        if (Y > 0.0 && Y <= y1)
        {
            static const std::vector<Real> coeff1 = {
            2.372066e-03, -2.184474e-03, -5.428215e-02,
            2.932582e+00, -9.519733e+01, 2.027762e+03,
            -2.936207e+04, 2.936498e+05, -2.027765e+06,
            9.488074e+06, -2.870793e+07, 5.066263e+07,
            -3.959145e+07,
            };
            polynomial_value = polyEval(coeff1, Y);
        }
        else
        {
            static const std::vector<Real> coeff2 = {
            2.373419e-03, -2.741400e-03, 1.415229e-03,
            -3.237809e-03, 7.040568e-03, -9.918119e-03,
            8.545716e-03, -1.622058e-05, -8.788399e-03,
            9.634837e-03, -4.033092e-03, -2.133763e-05,
            3.376056e-04,
            };
            polynomial_value = polyEval(coeff2, Y);
        }

        if (Y > channel_height / 2.0 || Y < 0.0)
        {
            std::cout << "position[1]=" << position[1] << std::endl;
            std::cout << "Y=" << Y << std::endl;
            std::cout << "polynomial_value=" << polynomial_value << std::endl;
            std::cout << "Stop" << std::endl;
            std::cout << "=================" << std::endl;
            std::cin.get();
        }

        temp_in_turbu_k = polynomial_value;
    }
    return (position[0] < 0.0) ? temp_in_turbu_k : turbu_k_original;
}

Real kOmega_InflowTurbulentCondition::getTurbulentInflowTemporaryEpsilon(Vecd &position, Real &turbu_k, Real turbu_E)
{
    Real temp_in_turbu_E = C_mu_75_for_omega_inlet_ * pow(turbu_k, 1.5) / TurbulentLength_;
    Real turbu_E_original = turbu_E;
    if (type_turbu_inlet_omega_ == 1)
    {
        Real channel_height = CharacteristicLength_;
        Real Y = (position[1] < channel_height / 2.0) ? position[1] : channel_height - position[1];

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

    return (position[0] < 0.0) ? temp_in_turbu_E : turbu_E_original;
}

Real kOmega_InflowTurbulentCondition::getTurbulentInflowOmega(Vecd& position, Vecd& velocity, Real& turbu_omega)
{
    Real temp_in_turbu_omega = 0.0;
    Real turbu_omega_original = turbu_omega;
    if (type_turbu_inlet_omega_ == 2)
    {
        const Real channel_height = CharacteristicLength_;
        Real Y = (position[1] <= channel_height * 0.5)
            ? position[1]
            : channel_height - position[1];
        const Real y1 = 0.2;

        Real polynomial_value = 0.0;
        if (Y > 0.0 && Y <= y1)
        {;
            polynomial_value = 2.2337739892e-01 / (Y + -9.5734226803e-07) + 4.9871467182e-02;
        }
        else
        {
            static const std::vector<Real> coeff2 = {
             5.610110e+00, -5.994682e+01, 3.670464e+02,
    -        1.398382e+03, 3.311412e+03, -4.223060e+03,
             -9.007083e+01, 1.109756e+04, -2.188842e+04,
             2.287992e+04, -1.420357e+04, 4.953297e+03,
             -7.511954e+02,
            };
            polynomial_value = polyEval(coeff2, Y);
        }

        if (Y < 0.0 || Y > channel_height * 0.5)
        {
            std::cout << "position[1] = " << position[1] << '\n'
                << "Y = " << Y << '\n'
                << "omega = " << polynomial_value << '\n'
                << "=================\n";
            std::cin.get();
        }

        temp_in_turbu_omega = polynomial_value;
    }
    return (position[0] < 0.0) ? temp_in_turbu_omega : turbu_omega_original;
}

Real kOmega_InflowTurbulentCondition::polyEval(const std::vector<Real>& a, Real x)
{
    const int n = static_cast<int>(a.size());
    if (n == 0)
    {
        std::cout << "size of coefficient = 0 " << '\n'
            << "=================\n";
        std::cin.get();
    }
    Real result = a[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        result = result * x + a[i];
    }
    return result;
}
}
}

}
