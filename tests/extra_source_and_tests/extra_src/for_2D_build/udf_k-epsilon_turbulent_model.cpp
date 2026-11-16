
#include "udf_k-epsilon_turbulent_model.hpp"
namespace SPH
{

namespace fluid_dynamics
{
namespace udf
{
using TurbuIntegration2ndHalfWithWallDissipativeRiemann = ComplexInteraction<Integration2ndHalf<Inner<>, Contact<Wall>>, DissipativeRiemannSolver>;

kEpsilon_TurbulentClosureCoefficient::kEpsilon_TurbulentClosureCoefficient()
    : Karman_(0.41), turbu_const_E_(9.8), C_mu_(0.09), turbulent_intensity_(5.0e-2),
      sigma_k_(1.0), C_l_(1.44), C_2_(1.92), sigma_E_(1.3), turbulent_length_ratio_for_epsilon_inlet_(0.07)
{
    inv_turbu_E_ = 1.0 / turbu_const_E_;
    C_mu_25_ = pow(C_mu_, 0.25);
    C_mu_75_ = pow(C_mu_, 0.75);
}

kEpsilon_kTransportEquationInner::kEpsilon_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, bool is_STL)
    : kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      dk_dt_(particles_->registerStateVariableData<Real>("ChangeRateOfTKE")),
      dk_dt_without_dissipation_(particles_->registerStateVariableData<Real>("ChangeRateOfTKEWithoutDissipation")),
      k_production_(particles_->registerStateVariableData<Real>("K_Production")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      turbu_k_(particles_->registerStateVariableData<Real>("TurbulenceKineticEnergy", Real(initial_values[0]))),
      turbu_epsilon_(particles_->registerStateVariableData<Real>("TurbulentDissipation", Real(initial_values[1]))),
      turbu_mu_(particles_->registerStateVariableData<Real>("TurbulentViscosity", Real(initial_values[2]))),
      turbu_strain_rate_(particles_->getVariableDataByName<Matd>("TurbulentStrainRate")),
      is_extra_viscous_dissipation_(particles_->registerStateVariableData<int>("TurbulentExtraViscousDissipation", is_extr_visc_dissipa)),
      is_STL_(is_STL),
      turbu_indicator_(particles_->registerStateVariableData<int>("TurbulentIndicator")),
      k_diffusion_(particles_->registerStateVariableData<Real>("K_Diffusion"))
{
    particles_->addEvolvingVariable<Real>("ChangeRateOfTKE");
    particles_->addEvolvingVariable<Real>("ChangeRateOfTKEWithoutDissipation");

    particles_->addEvolvingVariable<Real>("K_Production");
    particles_->addVariableToWrite<Real>("K_Production");

    particles_->addEvolvingVariable<Real>("TurbulenceKineticEnergy");
    particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

    particles_->addEvolvingVariable<Real>("TurbulentViscosity");
    particles_->addVariableToWrite<Real>("TurbulentViscosity");

    particles_->addEvolvingVariable<Real>("TurbulentDissipation");
    particles_->addVariableToWrite<Real>("TurbulentDissipation");

    particles_->addEvolvingVariable<Matd>("TurbulentStrainRate");
    particles_->addVariableToWrite<Matd>("TurbulentStrainRate");

    particles_->addEvolvingVariable<Real>("K_Diffusion");
    particles_->addVariableToWrite<Real>("K_Diffusion");

    particles_->addVariableToWrite<Real>("ChangeRateOfTKE");

    particles_->addEvolvingVariable<int>("TurbulentIndicator");
    particles_->addVariableToWrite<int>("TurbulentIndicator");

}

void kEpsilon_kTransportEquationInner::update(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_mu_i = turbu_mu_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real k_diffusion = k_diffusion_[index_i];
    Real k_dissipation = turbu_epsilon_[index_i];
    Matd vel_grad_i = velocity_gradient_[index_i];

    dk_dt_[index_i] = 0.0;
    dk_dt_without_dissipation_[index_i] = 0.0;
    Matd strain_rate = Matd::Zero();
    Matd Re_stress = Matd::Zero();

    Real k_production(0.0);

    strain_rate = 0.5 * (vel_grad_i.transpose() + vel_grad_i);

    Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();

    Matd k_production_matrix = Re_stress.array() * vel_grad_i.array();
    k_production = k_production_matrix.sum();

    dk_dt_[index_i] = k_production - k_dissipation + k_diffusion;
    dk_dt_without_dissipation_[index_i] = k_production + k_diffusion;

    if (is_near_wall_P1_[index_i] != 1)
        k_production_[index_i] = k_production;

    turbu_strain_rate_[index_i] = strain_rate;

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

kEpsilon_TKE_Diffusion::kEpsilon_TKE_Diffusion(BaseInnerRelation &inner_relation)
    : kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      k_diffusion_(particles_->getVariableDataByName<Real>("K_Diffusion")) {}

void kEpsilon_TKE_Diffusion::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_k_i = turbu_k_[index_i];

    Real mu_eff_i = turbu_mu_[index_i] / sigma_k_ + mu_;

    Real k_derivative(0.0);
    Real k_lap(0.0);

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = turbu_mu_[index_j] / sigma_k_ + mu_;
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
    }
    k_diffusion_[index_i] = k_lap;
}

kEpsilon_epsilonTransportEquationInner::kEpsilon_epsilonTransportEquationInner(BaseInnerRelation &inner_relation, bool is_STL)
    : kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      depsilon_dt_(particles_->registerStateVariableData<Real>("ChangeRateOfTDR")),
      depsilon_dt_without_dissipation_(particles_->registerStateVariableData<Real>("ChangeRateOfTDRWithoutDissp")),
      ep_production(particles_->registerStateVariableData<Real>("Ep_Production")),
      ep_dissipation_(particles_->registerStateVariableData<Real>("Ep_Dissipation_")),
      ep_diffusion_(particles_->registerStateVariableData<Real>("Ep_Diffusion_")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")),
      k_production_(particles_->getVariableDataByName<Real>("K_Production")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      is_STL_(is_STL)
{
    particles_->addEvolvingVariable<Real>("ChangeRateOfTDR");
    particles_->addVariableToWrite<Real>("ChangeRateOfTDR");

    particles_->addEvolvingVariable<Real>("Ep_Production");
    particles_->addVariableToWrite<Real>("Ep_Production");

    particles_->addEvolvingVariable<Real>("Ep_Dissipation_");
    particles_->addVariableToWrite<Real>("Ep_Dissipation_");

    particles_->addEvolvingVariable<Real>("Ep_Diffusion_");
    particles_->addVariableToWrite<Real>("Ep_Diffusion_");
}

void kEpsilon_epsilonTransportEquationInner::update(size_t index_i, Real dt)
{
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_epsilon_i = turbu_epsilon_[index_i];
    Real epsilon_diffusion = ep_diffusion_[index_i];
    Real k_production_i = k_production_[index_i];

    depsilon_dt_[index_i] = 0.0;
    depsilon_dt_without_dissipation_[index_i] = 0.0;
    Real epsilon_production(0.0);
    Real epsilon_dissipation(0.0);

    epsilon_production = C_l_ * turbu_epsilon_i * k_production_i / turbu_k_i;
    epsilon_dissipation = C_2_ * turbu_epsilon_i * turbu_epsilon_i / turbu_k_i;

    depsilon_dt_[index_i] = epsilon_production - epsilon_dissipation + epsilon_diffusion;
    depsilon_dt_without_dissipation_[index_i] = epsilon_production + epsilon_diffusion;

    ep_production[index_i] = epsilon_production;
    ep_dissipation_[index_i] = epsilon_dissipation;

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

kEpsilon_TDR_Diffusion::kEpsilon_TDR_Diffusion(BaseInnerRelation &inner_relation)
    : kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      ep_diffusion_(particles_->getVariableDataByName<Real>("Ep_Diffusion_")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")) {}

void kEpsilon_TDR_Diffusion::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Real turbu_epsilon_i = turbu_epsilon_[index_i];

    Real mu_eff_i = turbu_mu_[index_i] / sigma_E_ + mu_;

    Real epsilon_derivative(0.0);
    Real epsilon_lap(0.0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = turbu_mu_[index_j] / sigma_E_ + mu_;
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        epsilon_derivative = (turbu_epsilon_i - turbu_epsilon_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        epsilon_lap += 2.0 * mu_harmo * epsilon_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
    }
    ep_diffusion_[index_i] = epsilon_lap;
}

kEpsilon_TurbulentEddyViscosity::
    kEpsilon_TurbulentEddyViscosity(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation")),
      wall_Y_plus_(particles_->getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(particles_->getVariableDataByName<Real>("WallYstar")),
      viscosity_(sph_body_->getMaterialProperty<Viscosity>()),
      mu_(viscosity_.ReferenceViscosity()) {}

void kEpsilon_TurbulentEddyViscosity::update(size_t index_i, Real dt)
{
    turbu_mu_[index_i] = rho_[index_i] * C_mu_ * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
}

kEpsilon_InflowTurbulentCondition::kEpsilon_InflowTurbulentCondition(BodyPartByCell &body_part, Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet)
    : BaseFlowBoundaryCondition(body_part), type_turbu_inlet_(type_turbu_inlet),
      relaxation_rate_(relaxation_rate),
      CharacteristicLength_(CharacteristicLength),
      turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_epsilon_(particles_->getVariableDataByName<Real>("TurbulentDissipation"))
{
    TurbulentLength_ = turbulent_length_ratio_for_epsilon_inlet_ * CharacteristicLength_;
}

void kEpsilon_InflowTurbulentCondition::update(size_t index_i, Real dt)
{
    Real target_in_turbu_k = getTurbulentInflowK(pos_[index_i], vel_[index_i], turbu_k_[index_i]);
    turbu_k_[index_i] += relaxation_rate_ * (target_in_turbu_k - turbu_k_[index_i]);
    Real target_in_turbu_E = getTurbulentInflowE(pos_[index_i], turbu_k_[index_i], turbu_epsilon_[index_i]);
    turbu_epsilon_[index_i] += relaxation_rate_ * (target_in_turbu_E - turbu_epsilon_[index_i]);
}

Real kEpsilon_InflowTurbulentCondition::getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k)
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
    if (position[0] < 0.0)
    {
        turbu_k_original = temp_in_turbu_k;
    }
    return turbu_k_original;
}

Real kEpsilon_InflowTurbulentCondition::getTurbulentInflowE(Vecd &position, Real &turbu_k, Real &turbu_E)
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

        if (Y > channel_height / 2.0 || Y < 0.0)
        {
            std::cout << "position[1]=" << position[1] << std::endl;
            std::cout << "Y=" << Y << std::endl;
            std::cout << "polynomial_value=" << polynomial_value << std::endl;
            std::cout << "Stop" << std::endl;
            std::cout << "=================" << std::endl;
            std::cin.get();
        }

        temp_in_turbu_E = polynomial_value;
    }
    if (position[0] < 0.0)
    {
        turbu_E_original = temp_in_turbu_E;
    }
    return turbu_E_original;
}

kEpsilon_StandardWallFunctionCorrection::
    kEpsilon_StandardWallFunctionCorrection(BaseInnerRelation &inner_relation,
                                            BaseContactRelation &contact_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      wall_Y_plus_(particles_->registerStateVariableData<Real>("WallYplus")),
      wall_Y_star_(particles_->registerStateVariableData<Real>("WallYstar")),
      velo_tan_(particles_->registerStateVariableData<Real>("TangentialVelocity")),
      velo_friction_(particles_->registerStateVariableData<Vecd>("FrictionVelocity")),
      y_p_(particles_->getVariableDataByName<Real>("Y_P")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")), rho_(particles_->getVariableDataByName<Real>("Density")),
      viscosity_(sph_body_->getMaterialProperty<Viscosity>()),
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
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }

    particles_->addEvolvingVariable<Real>("Y_P");
    particles_->addVariableToWrite<Real>("Y_P");

    particles_->addEvolvingVariable<Real>("WallYplus");
    particles_->addVariableToWrite<Real>("WallYplus");

    particles_->addEvolvingVariable<Real>("WallYstar");
    particles_->addVariableToWrite<Real>("WallYstar");

    particles_->addEvolvingVariable<Real>("TangentialVelocity");
    particles_->addVariableToWrite<Real>("TangentialVelocity");

    particles_->addEvolvingVariable<Vecd>("FrictionVelocity");
    particles_->addVariableToWrite<Vecd>("FrictionVelocity");
};

void kEpsilon_StandardWallFunctionCorrection::interaction(size_t index_i, Real dt)
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

        wall_Y_star_[index_i] = y_p_constant_i * C_mu_wf_25_ * turbu_k_i_05 / nu_i;

        Real velo_fric_mag = 0.0;
        Real velo_tan_mag = 0.0;

        velo_tan_mag = abs(e_i_nearest_tau.dot(vel_i));
        velo_tan_[index_i] = velo_tan_mag;

        if (wall_Y_star_[index_i] != static_cast<Real>(wall_Y_star_[index_i]))
        {
            std::cout << "y* is not a real value, please check" << std::endl;
            std::cin.get();
        }

        Real u_star = get_dimensionless_velocity(wall_Y_star_[index_i], current_time, 0.0, 0);
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

                    Real denominator_log_law_j = C_mu_wf_25_ * turbu_k_i_05 * Karman_ * y_p_j;

                    Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
                    Real y_star_j = C_mu_wf_25_ * turbu_k_i_05 * y_p_j / nu_i;
                    Real u_star_j = get_dimensionless_velocity(y_star_j, current_time, 0.0, 0);
                    Real fric_vel_mag_j = sqrt(C_mu_wf_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

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
                        epsilon_p_j = C_mu_wf_75_ * turbu_k_i_15 / (Karman_ * y_p_j);
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

}

}

}
