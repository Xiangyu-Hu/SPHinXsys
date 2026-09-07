
#include "udf_common_turbulence_model.hpp"
namespace SPH
{

namespace fluid_dynamics
{

namespace udf
{

WallFunctionCoefficient::WallFunctionCoefficient()
    : Karman_(0.41), turbu_const_E_(9.8), C_mu_wf_(0.09),
      start_time_laminar_(2.0), y_star_threshold_laminar_(11.225)
{
    C_mu_wf_25_ = pow(C_mu_wf_, 0.25);
    C_mu_wf_75_ = pow(C_mu_wf_, 0.75);
    inv_turbu_E_ = 1.0 / turbu_const_E_;
}

Real WallFunction::get_distance_from_P_to_wall(Real y_p_constant)
{
    return y_p_constant;
}

Real WallFunction::get_dimensionless_velocity(Real y_star, Real time, Real u_star_previous, int is_blended)
{
    Real dimensionless_velocity = 0.0;
    if (is_blended && time > start_time_laminar_)
    {
        dimensionless_velocity = Spalding_wall_function(y_star, u_star_previous);
    }
    else
    {

        if (y_star < y_star_threshold_laminar_ && time > start_time_laminar_)
        {
            dimensionless_velocity = laminar_law_wall_function(y_star);
        }
        else
        {
            if (std::isnan(dimensionless_velocity) || std::isinf(dimensionless_velocity))
            {
                std::cout << "u_star=" << dimensionless_velocity << std::endl;
                std::cout << "y_star=" << y_star << std::endl;
                std::cout << "system pause" << std::endl;

                std::cin.get();
            }
            dimensionless_velocity = log_law_wall_function(y_star);
        }
    }
    if (std::isnan(dimensionless_velocity) || std::isinf(dimensionless_velocity))
    {
        std::cout << "u_star=" << dimensionless_velocity << std::endl;
        std::cout << "y_star=" << y_star << std::endl;
    }

    return dimensionless_velocity;
}

Real WallFunction::get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity)
{
    Real vel_grad_mag = log_law_velocity_gradient(vel_fric_mag, denominator_log_law);
    return vel_grad_mag;
}

Real WallFunction::log_law_wall_function(Real y_star)
{
    Real u_star = abs(log(turbu_const_E_ * y_star) / Karman_);
    return u_star;
}

Real WallFunction::laminar_law_wall_function(Real y_star)
{
    Real u_star = y_star;
    return u_star;
}

Real WallFunction::Spalding_wall_function(Real y_star, Real u_star_guess)
{
    if (u_star_guess > 500.0 || u_star_guess <= TinyReal)
    {
        std::cout << "u_star_guess > 500.0 || u_star_guess <= TinyReal, please check." << std::endl;
        std::cout << "u_star_guess=" << u_star_guess << std::endl;
    }
    Real u_star = u_star_guess;
    int max_iter = 10;
    Real tolerance = 0.01;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        Real Karman_u_star = SMIN(Karman_ * u_star, 50.0);
        Real f = u_star + inv_turbu_E_ * (std::exp(Karman_u_star) - 1.0 - Karman_u_star - 0.5 * pow(Karman_u_star, 2) - 1.0 / 6.0 * pow(Karman_u_star, 3)) - y_star;
        Real df = 1.0 + inv_turbu_E_ * (Karman_ * std::exp(Karman_u_star) - Karman_ - Karman_ * Karman_u_star - 0.5 * Karman_ * pow(Karman_u_star, 2));
        u_star -= f / df;
        Real residue = std::abs(f / df);
        if (residue <= tolerance)
            break;
    }
    if (u_star > 500.0)
    {
        std::cout << "u+ is larger than 500, please check." << std::endl;
        std::cout << "u_star=" << u_star << std::endl;
    }
    return u_star;
}

Real WallFunction::log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law)
{
    return vel_fric_mag * vel_fric_mag / denominator_log_law;
}

Real WallFunction::laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity)
{
    return vel_fric_mag * vel_fric_mag / dynamic_viscosity;
}

kEpsilon_GetVelocityGradient<Inner<>>::kEpsilon_GetVelocityGradient(BaseInnerRelation &inner_relation, Real weight_sub)
    : kEpsilon_GetVelocityGradient<DataDelegateInner>(inner_relation),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
      weight_sub_nearwall_(weight_sub) {}

void kEpsilon_GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
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

void kEpsilon_GetVelocityGradient<Inner<>>::update(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] != 1)
    {
        velocity_gradient_[index_i] *= turbu_B_[index_i];
    }
}

TKEnergyForce<Inner<>>::TKEnergyForce(BaseInnerRelation &inner_relation)
    : TKEnergyForce<Base, DataDelegateInner>(inner_relation),
      test_k_grad_rslt_(this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}

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

TKEnergyForce<Contact<>>::TKEnergyForce(BaseContactRelation &contact_relation)
    : TKEnergyForce<Base, DataDelegateContact>(contact_relation),
      test_k_grad_rslt_(this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}

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

TurbuViscousForce<Inner<>>::TurbuViscousForce(BaseInnerRelation &inner_relation)
    : TurbuViscousForce<DataDelegateInner>(inner_relation),
      turbu_indicator_(this->particles_->template getVariableDataByName<int>("TurbulentIndicator")),
      is_extra_viscous_dissipation_(this->particles_->template getVariableDataByName<int>("TurbulentExtraViscousDissipation")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      dUdn_P_sublayer_(particles_->getVariableDataByName<Matd>("dUdnFromSublayer")),
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {}

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
        if (this->is_near_wall_P2_[index_i] == 10)
        {
            shear_stress_eij_corrected = shear_stress_eij;
        }
        shear_stress = (shear_stress - shear_stress_eij) + shear_stress_eij_corrected;

        Vecd force_j = 2.0 * mass_[index_i] * shear_stress * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];

        bool is_inner_SS_correction = true;
        if(is_inner_SS_correction)
        {
            if (is_near_wall_P1_[index_i] == 1)
            {
                Matd shear_stress_sublayer = mu_eff_i * (dUdn_P_sublayer_[index_i] + dUdn_P_sublayer_[index_i].transpose());
                Vecd kernel_gradient = (e_ij * inner_neighborhood.dW_ij_[n]);
                force_j = 2.0 * mass_[index_i] * (shear_stress_sublayer * kernel_gradient) * this->Vol_[index_j];
            }
            if (is_near_wall_P2_[index_i] == 10 && is_near_wall_P1_[index_i] != 1)
            {
                if (is_near_wall_P1_[index_j] == 1)
                {
                    Matd shear_stress_sublayer = mu_eff_j * (dUdn_P_sublayer_[index_j] + dUdn_P_sublayer_[index_j].transpose());
                    Vecd kernel_gradient = (e_ij * inner_neighborhood.dW_ij_[n]);
                    force_j = 2.0 * mass_[index_i] * (shear_stress_sublayer * kernel_gradient) * this->Vol_[index_j];
                }
            }
        }
        force += force_j;
    }
    viscous_force_[index_i] = force / rho_[index_i];
}

TurbuViscousForce<Contact<Wall>>::TurbuViscousForce(BaseContactRelation &wall_contact_relation)
    : BaseTurbuViscousForceWithWall(wall_contact_relation),
      wall_particle_spacing_(wall_contact_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing()),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")),
      is_blended_(particles_->getVariableDataByName<int>("TurbulentWallTreatmentType")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      friction_velocity_from_sublayer_(particles_->getVariableDataByName<Real>("FrictionVelocityFromSublayer")),
      turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
      B_only_wall_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrixOnlyWall")) {}

void TurbuViscousForce<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    if (this->is_near_wall_P1_[index_i] != 1)
        return;

    Real vel_fric_mag_previous = velo_friction_[index_i].norm();
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

            Real u_star_previous = vel_i_tau_mag / vel_fric_mag_previous;
            if ((u_star_previous > 100.0 || u_star_previous <= TinyReal) && current_time > start_time_laminar_)
            {
                u_star_previous = wall_Y_star_[index_i] + 10.0 * TinyReal;
            }

            Real y_p_j = get_distance_from_P_to_wall(y_p_constant_i);
            Real y_star_j = rho_i * C_mu_wf_25_ * turbu_k_i_05 * y_p_j / molecular_viscosity_;
            Real u_star_j = get_dimensionless_velocity(y_star_j, current_time, u_star_previous, is_blended_[index_i]);
            Real fric_vel_mag_j = sqrt(C_mu_wf_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);
            if (is_near_wall_P1_[index_i] == 1)
            {
                fric_vel_mag_j = friction_velocity_from_sublayer_[index_i];
            }
            Real WSS_tn_mag_j = rho_i * fric_vel_mag_j * fric_vel_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));

            WSS_j_tn(0, 0) = 0.0;
            WSS_j_tn(0, 1) = WSS_tn_mag_j;
            WSS_j_tn(1, 0) = 0.0;
            WSS_j_tn(1, 1) = 0.0;
            WSS_j = Q.transpose() * WSS_j_tn * Q;

            bool is_add_B_only_wall_to_WSS_P1 = false;
            Vecd force_j = Vecd::Zero();
            if (is_add_B_only_wall_to_WSS_P1)
            {
                if (is_near_wall_P1_[index_i] == 1)
                {
                    Vecd kernel_gradient = (e_ij * contact_neighborhood.dW_ij_[n]);
                    force_j = 2.0 * mass_[index_i] * (WSS_j * kernel_gradient) * this->Vol_[index_j] / rho_i;
                }
                else
                {
                    force_j = 2.0 * mass_[index_i] * WSS_j * e_ij * contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
                }
            }
            else
            {
                force_j = 2.0 * mass_[index_i] * WSS_j * e_ij * contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;
            }

            force += force_j;
        }
    }
    viscous_force_[index_i] += force;
}

TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax<Real>>(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      speed_ref_turbu_(U_max), advectionCFL_(advectionCFL),
      turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      fluid_(DynamicCast<Fluid>(this, sph_body_->getMatterMaterial())),
      viscosity_(sph_body_->getMaterialProperty<Viscosity>())
{
    Real viscous_speed = viscosity_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
    speed_ref_turbu_ = SMAX(viscous_speed, speed_ref_turbu_);
}

Real TurbulentAdvectionTimeStepSize::reduce(size_t index_i, Real dt)
{
    Real turbu_viscous_speed = (viscosity_.ReferenceViscosity() + turbu_mu_[index_i]) / fluid_.ReferenceDensity() / smoothing_length_min_;
    Real turbu_viscous_speed_squire = turbu_viscous_speed * turbu_viscous_speed;
    Real vel_n_squire = vel_[index_i].squaredNorm();
    Real vel_bigger = SMAX(turbu_viscous_speed_squire, vel_n_squire);

    return vel_bigger;
}

Real TurbulentAdvectionTimeStepSize::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * smoothing_length_min_ / (SMAX(speed_max, speed_ref_turbu_) + TinyReal);
}

JudgeIsNearWall::
    JudgeIsNearWall(BaseInnerRelation &inner_relation,
                    BaseContactRelation &contact_relation, Real constant_y_p)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      distance_to_dummy_interface_(particles_->registerStateVariableData<Real>("DistanceToDummyInterface")),
      distance_to_dummy_interface_up_average_(particles_->registerStateVariableData<Real>("DistanceToDummyInterfaceUpAver")),
      is_near_wall_P1_(particles_->registerStateVariableData<int>("IsNearWallP1")),
      is_near_wall_P2_(particles_->registerStateVariableData<int>("IsNearWallP2")),
      index_nearest_(particles_->registerStateVariableData<int>("NearestIndex")),
      e_nearest_tau_(particles_->registerStateVariableData<Vecd>("WallNearestTangentialUnitVector")),
      e_nearest_normal_(particles_->registerStateVariableData<Vecd>("WallNearestNormalUnitVector")),
      y_p_(particles_->registerStateVariableData<Real>("Y_P")),
      constant_y_p_(constant_y_p),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      fluid_particle_spacing_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing()),
      wall_particle_spacing_(contact_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
};

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

        if (Dimensions == 2)
        {
            e_i_nearest_tau[0] = e_i_nearest_n[1];
            e_i_nearest_tau[1] = e_i_nearest_n[0] * (-1.0);
        }

        if (r_dmy_itfc_n_sum <= 0.0)
        {
            std::cout << "r_dmy_itfc_n_sum is almost zero" << std::endl;
            std::cout << "count=" << count_average << std::endl;
            std::cin.get();
        }

        distance_to_dummy_interface_up_average_[index_i] = r_dmy_itfc_n_sum / ttl_weight;

        index_nearest_[index_i] = id_nearest_j;
        e_nearest_normal_[index_i] = e_i_nearest_n;
        e_nearest_tau_[index_i] = e_i_nearest_tau;
        distance_to_dummy_interface_[index_i] = r_dummy_normal;
    }
}

void JudgeIsNearWall::update(size_t index_i, Real dt)
{
    is_near_wall_P1_[index_i] = 0;
    y_p_[index_i] = 0.0;
    if (is_near_wall_P2_[index_i] == 10)
    {
        y_p_[index_i] = constant_y_p_ + fluid_particle_spacing_ / 2.0;
        Real distance = distance_to_dummy_interface_[index_i];
        if (distance < 1.0 * fluid_particle_spacing_)
        {
            is_near_wall_P1_[index_i] = 1;
            y_p_[index_i] = constant_y_p_;
        }
    }
}

ConstrainNormalVelocityInRegionP::
    ConstrainNormalVelocityInRegionP(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
      e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")),
      wall_Y_star_(particles_->getVariableDataByName<Real>("WallYstar")) {}

void ConstrainNormalVelocityInRegionP::update(size_t index_i, Real dt)
{
    if (is_near_wall_P1_[index_i] == 1 && wall_Y_star_[index_i] >= y_star_threshold_laminar_)
    {
        vel_[index_i] = vel_[index_i] - (vel_[index_i].dot(e_nearest_normal_[index_i])) * e_nearest_normal_[index_i];
    }
}

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

void TurbulentLinearGradientCorrectionMatrix<Inner<>>::update(size_t index_i, Real dt)
{
    Real det_sqr = SMAX(turbu_alpha_ - turbu_B_[index_i].determinant(), Real(0));
    Matd inverse = turbu_B_[index_i].inverse();
    Real weight1_ = turbu_B_[index_i].determinant() / (turbu_B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (turbu_B_[index_i].determinant() + det_sqr);
    turbu_B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();

    if (is_near_wall_P1_[index_i] == 1)
    {
        det_sqr = SMAX(turbu_alpha_ - B_only_wall_[index_i].determinant(), Real(0));
        inverse = B_only_wall_[index_i].inverse();
        weight1_ = B_only_wall_[index_i].determinant() / (B_only_wall_[index_i].determinant() + det_sqr);
        weight2_ = det_sqr / (B_only_wall_[index_i].determinant() + det_sqr);
        B_only_wall_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
    }

}

TurbulentLinearGradientCorrectionMatrix<Contact<>>::
TurbulentLinearGradientCorrectionMatrix(BaseContactRelation& contact_relation)
    : TurbulentLinearGradientCorrectionMatrix<DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
}

void TurbulentLinearGradientCorrectionMatrix<Contact<>>::interaction(size_t index_i, Real dt)
{
    B_only_wall_[index_i] = ZeroData<Matd>::value;
    Matd local_configuration = Eps * Matd::Identity();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real* Vol_k = contact_Vol_[k];
        Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
    }
    B_only_wall_[index_i] += local_configuration;
}

GetLimiterOfTransportVelocityCorrection::
    GetLimiterOfTransportVelocityCorrection(SPHBody &sph_body, Real slope)
    : LocalDynamics(sph_body),
      h_ref_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()),
      zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("KernelGradientIntegral")),
      slope_(slope),
      limiter_tvc_(particles_->registerStateVariableData<Real>("LimiterOfTVC")){}

void GetLimiterOfTransportVelocityCorrection::update(size_t index_i, Real dt)
{
    Real squared_norm = zero_gradient_residue_[index_i].squaredNorm();
    limiter_tvc_[index_i] = SMIN(slope_ * squared_norm * h_ref_ * h_ref_, Real(1));
}

NonDimensionalisePressure::
    NonDimensionalisePressure(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      p_dimensionless_(particles_->registerStateVariableData<Real>("PressureDimensionless")) {}

void NonDimensionalisePressure::update(size_t index_i, Real dt)
{
    p_dimensionless_[index_i] = p_[index_i] / ((rho_[index_i] - 1.0));
}

}

}

}
