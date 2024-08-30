//#pragma once
#include "k-epsilon_turbulent_model.hpp"
namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
	using TurbuIntegration2ndHalfWithWallDissipativeRieman = ComplexInteraction<Integration2ndHalf<Inner<>, Contact<Wall>>, DissipativeRiemannSolver>;
//=================================================================================================//
	BaseTurbuClosureCoeff::BaseTurbuClosureCoeff()
		: Karman_(0.41), turbu_const_E_(9.8), C_mu_(0.09), turbulent_intensity_(5.0e-2), 
		sigma_k_(1.0),  C_l_(1.44), C_2_(1.92), sigma_E_(1.3), turbulent_length_ratio_for_epsilon_inlet_(0.07),
		start_time_laminar_(0.0), y_star_threshold_laminar_(11.225)
	{
		C_mu_25_ = pow(C_mu_, 0.25);
		C_mu_75_ = pow(C_mu_, 0.75);
	}
//=================================================================================================//
	Real WallFunction::get_distance_from_P_to_wall(Real y_p_constant)
	{
		//** Check the distance. *  
		//if (y_p_constant < 0.05 * dp_wall)
		//{
		//	std::cout << "y_p_j < 0.05 * wall_particle_spacing_" << std::endl;
		//	std::cin.get();
		//}
		//y_p_j = abs(e_j_n.dot(r_ij * e_ij)) - 0.5 * wall_particle_spacing_;
		
		//** Use the constant y_p strategy. *  
		return y_p_constant;
	}
	//=================================================================================================//
	Real WallFunction:: get_dimensionless_velocity(Real y_star)
	{
		Real dimensionless_velocity = 0.0;
		if (y_star < y_star_threshold_laminar_ && GlobalStaticVariables::physical_time_ > start_time_laminar_)
		{
			dimensionless_velocity = laminar_law_wall_functon(y_star);
		}
		else
		{
			dimensionless_velocity = log_law_wall_functon(y_star);
		}
		if (std::isnan(dimensionless_velocity) || std::isinf(dimensionless_velocity))
		{
			std::cout << "u_star=" << dimensionless_velocity << std::endl;
			std::cout << "y_star=" << y_star << std::endl;
			std::cin.get();
		}
		//if (dimensionless_velocity<0.0)
		//{
		//	std::cout << "dimensionless_velocity<0.0" << dimensionless_velocity << std::endl;
		//	std::cin.get();
		//}

		return dimensionless_velocity;
	}
	//=================================================================================================//
	Real WallFunction:: get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity)
	{
		Real vel_grad_mag = log_law_velocity_gradient(vel_fric_mag, denominator_log_law);
		return vel_grad_mag;
	}
	//=================================================================================================//
	Real WallFunction::log_law_wall_functon(Real y_star)
	{
		//** u_star should larger than 0 *
		Real u_star = abs(log(turbu_const_E_ * y_star ) / Karman_);
		return u_star;
	}
	//=================================================================================================//
	Real WallFunction::laminar_law_wall_functon(Real y_star)
	{
		Real u_star = y_star ;
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
	GetVelocityGradient<Inner<>>::GetVelocityGradient(BaseInnerRelation& inner_relation, Real weight_sub)
		: GetVelocityGradient<DataDelegateInner>(inner_relation),
		velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
		B_(*particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
		turbu_B_(*particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
		weight_sub_nearwall_(weight_sub)
	{
		this->particles_->addVariableToSort<Matd>("TurbulentVelocityGradient");
		this->particles_->addVariableToWrite<Matd>("TurbulentVelocityGradient");
	}
	//=================================================================================================//
	void GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
    {
		//** The near wall velo grad is updated in wall function part *
		if (is_near_wall_P1_[index_i] != 1)
		{
			velocity_gradient_[index_i] = Matd::Zero();
			Vecd vel_i = vel_[index_i];
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
				
				Real r_ij = inner_neighborhood.r_ij_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];
				if (is_near_wall_P2_[index_i] == 10 && is_near_wall_P1_[index_j] == 1)
				{
					//Vecd vel_ps = vel_[index_j] + 0.5 * r_ij * velocity_gradient_[index_j] * e_ij ;
					//Vecd vel_ps =  0.5 * vel_[index_j] + 0.5 * vel_[index_i] ;
                    //velocity_gradient_[index_i] += -2.0 * (vel_i - vel_ps) * nablaW_ijV_j.transpose();
					//std::cout<<"Calculate velocity gradient according to interpolation!"<<std::endl;
					//std::cin.get();
					//Vecd vel_diff = velocity_gradient_[index_j] * r_ij * e_ij;
					//Vecd vel_diff = velocity_gradient_[index_j] * e_ij;
					
					//Real factor_B = B_[index_i].norm() / turbu_B_[index_i].norm();
					Matd P1 = - (vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
					Vecd vel_diff = velocity_gradient_[index_j] * r_ij * e_ij;
					Matd P2 = - vel_diff * nablaW_ijV_j.transpose();
					
					//velocity_gradient_[index_i] += - weight_sub_nearwall_ * vel_diff * nablaW_ijV_j.transpose();
					velocity_gradient_[index_i] +=  ( 1 - weight_sub_nearwall_ ) * P1 +  weight_sub_nearwall_* P2;
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
			//velocity_gradient_[index_i] *= B_[index_i];
		    velocity_gradient_[index_i] *= turbu_B_[index_i];
		}
	    
	}
	//=================================================================================================//
	GetVelocityGradient<Contact<Wall>>::GetVelocityGradient(BaseContactRelation& contact_relation)
		: InteractionWithWall<GetVelocityGradient>(contact_relation),
		velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient"))
	{
		this->particles_->addVariableToSort<Matd>("Velocity_Gradient_Wall");
		this->particles_->addVariableToWrite<Matd>("Velocity_Gradient_Wall");
	}
	//=================================================================================================//
	void GetVelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
	{
		//** The near wall velo grad is updated in wall function part *
		if (is_near_wall_P1_[index_i] != 1)
		{
			Vecd vel_i = vel_[index_i];
			for (size_t k = 0; k < DataDelegateContact::contact_configuration_.size(); ++k)
			{
				//StdLargeVec<Vecd>& vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*DataDelegateContact::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n]* this->Vol_[index_j] * contact_neighborhood.e_ij_[n];

					//velocity_gradient_[index_i] += -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
					velocity_gradient_[index_i] += -1.0 * (vel_i ) * nablaW_ijV_j.transpose();
				}
			}

		}
	}
//=================================================================================================//
	TransferVelocityGradient::
		TransferVelocityGradient(SPHBody& sph_body)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body), 
		is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
		velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
		vel_grad_(*this->particles_->template registerSharedVariable<Matd>("VelocityGradient")){}
	//=================================================================================================//
	void TransferVelocityGradient::update(size_t index_i, Real dt)
	{
		velocity_gradient_[index_i] = Matd::Zero();
		if (is_near_wall_P1_[index_i] != 1)
		{
			velocity_gradient_[index_i]=vel_grad_[index_i]; //** Transfer value, but exclude P region */
		}
	}
//=================================================================================================//
	K_TurtbulentModelInner::K_TurtbulentModelInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values, int is_extr_visc_dissipa)
		: BaseTurtbulentModel<Base, DataDelegateInner>(inner_relation),
		dk_dt_(*particles_->registerSharedVariable<Real>("ChangeRateOfTKE")),
		dk_dt_without_dissipation_(*particles_->registerSharedVariable<Real>("ChangeRateOfTKEWithoutDissipation")),
		k_production_(*particles_->registerSharedVariable<Real>("K_Production")),
		is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
		velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation")),
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		turbu_strain_rate_(*particles_->getVariableDataByName<Matd>("TurbulentStrainRate")),
		is_extra_viscous_dissipation_(*particles_->registerSharedVariable<int>("TurbulentExtraViscousDissipation")),
		turbu_indicator_(*particles_->registerSharedVariable<int>("TurbulentIndicator")),
		k_diffusion_(*particles_->registerSharedVariable<Real>("K_Diffusion")),
		vel_x_(*particles_->registerSharedVariable<Real>("Velocity_X"))
	{
		//particles_->registerSharedVariable(dk_dt_, "ChangeRateOfTKE");
		particles_->addVariableToSort<Real>("ChangeRateOfTKE");
		particles_->addVariableToSort<Real>("ChangeRateOfTKEWithoutDissipation");

		//particles_->registerSharedVariable(k_production_, "K_Production");
		particles_->addVariableToSort<Real>("K_Production");
		particles_->addVariableToWrite<Real>("K_Production");

		particles_->addVariableToSort<Real>("TurbulenceKineticEnergy");
		particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

		particles_->addVariableToSort<Real>("TurbulentViscosity");
		particles_->addVariableToWrite<Real>("TurbulentViscosity");

		particles_->addVariableToSort<Real>("TurbulentDissipation");
		particles_->addVariableToWrite<Real>("TurbulentDissipation");

		particles_->addVariableToSort<Matd>("TurbulentStrainRate");
		particles_->addVariableToWrite<Matd>("TurbulentStrainRate");


		//** Obtain Initial values for transport equations *		
		std::fill(turbu_k_.begin(), turbu_k_.end(), initial_values[0]);
		std::fill(turbu_epsilon_.begin(), turbu_epsilon_.end(), initial_values[1]);
		std::fill(turbu_mu_.begin(), turbu_mu_.end(), initial_values[2]);

		//** for test */
		//particles_->registerSharedVariable(k_diffusion_, "K_Diffusion");
		particles_->addVariableToSort<Real>("K_Diffusion");
		particles_->addVariableToWrite<Real>("K_Diffusion");

		particles_->addVariableToWrite<Real>("ChangeRateOfTKE");
		
		//particles_->registerSharedVariable(vel_x_, "Velocity_X");
		particles_->addVariableToSort<Real>("Velocity_X");

		//particles_->registerSharedVariable(turbu_indicator_, "TurbulentIndicator");
		particles_->addVariableToSort<int>("TurbulentIndicator");
		particles_->addVariableToWrite<int>("TurbulentIndicator");

		//particles_->registerSharedVariable(is_extra_viscous_dissipation_, "TurbulentExtraViscousDissipation");
		std::fill(is_extra_viscous_dissipation_.begin(), is_extra_viscous_dissipation_.end(), is_extr_visc_dissipa);
	}
	//=================================================================================================//
	void K_TurtbulentModelInner::interaction(size_t index_i, Real dt)
	{
		//Vecd vel_i = vel_[index_i];
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
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real mu_eff_j = turbu_mu_[index_j] / sigma_k_ + mu_;
			Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
			k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
			k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ij_[n]* this->Vol_[index_j] / rho_i;
		}
		strain_rate = 0.5 * (velocity_gradient_[index_i].transpose() + velocity_gradient_[index_i]);

		Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
		//Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i;

		Matd k_production_matrix = Re_stress.array() * velocity_gradient_[index_i].array();
		//** The near wall k production is updated in wall function part *
		if (is_near_wall_P1_[index_i] != 1)
			k_production_[index_i] = k_production_matrix.sum();
			
		k_production = k_production_[index_i];
		k_dissipation = turbu_epsilon_[index_i] ;

		//** Linearize the source term *
		//k_production = k_production_prior_[index_i];
		//k_dissipation = ( turbu_epsilon_prior_[index_i] / turbu_k_prior_[index_i] ) * turbu_k_[index_i];
		
		dk_dt_[index_i] = k_production - k_dissipation + k_lap;
		dk_dt_without_dissipation_[index_i] = k_production + k_lap;

		//** for test */
		k_diffusion_[index_i] = k_lap;
		vel_x_[index_i] = vel_[index_i][0];
		turbu_strain_rate_[index_i] = strain_rate;

	}
	//=================================================================================================//
	void K_TurtbulentModelInner::update(size_t index_i, Real dt)
	{
		turbu_k_[index_i] += dk_dt_[index_i] * dt;
		//** If use source term linearisation *
		//Real denominator = 1.0 + turbu_epsilon_[index_i] * dt / turbu_k_[index_i];
		//turbu_k_[index_i] +=  dk_dt_without_dissipation_[index_i] * dt;
		//turbu_k_[index_i] /=  denominator ;
	}
	//=================================================================================================//
	// void K_TurtbulentModelInner::update_prior_turbulent_value()
	// {
	// 	k_production_prior_ = k_production_;
	// 	turbu_k_prior_ = turbu_k_;
	// 	turbu_epsilon_prior_ = turbu_epsilon_;
	// }
//=================================================================================================//
	E_TurtbulentModelInner::E_TurtbulentModelInner(BaseInnerRelation& inner_relation)
		: BaseTurtbulentModel<Base, DataDelegateInner>(inner_relation),
		depsilon_dt_(*particles_->registerSharedVariable<Real>("ChangeRateOfTDR")),
		depsilon_dt_without_disspation_(*particles_->registerSharedVariable<Real>("ChangeRateOfTDRWithoutDissp")),
		ep_production(*particles_->registerSharedVariable<Real>("Ep_Production")),
		ep_dissipation_(*particles_->registerSharedVariable<Real>("Ep_Dissipation_")),
		ep_diffusion_(*particles_->registerSharedVariable<Real>("Ep_Diffusion_")),
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation")),
		k_production_(*particles_->getVariableDataByName<Real>("K_Production")),		
		is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1"))
	{
		//particles_->registerSharedVariable(depsilon_dt_, "ChangeRateOfTDR");
		particles_->addVariableToSort<Real>("ChangeRateOfTDR");
		particles_->addVariableToWrite<Real>("ChangeRateOfTDR");

		//particles_->registerSharedVariable(ep_production, "Ep_Production");
		particles_->addVariableToSort<Real>("Ep_Production");
		particles_->addVariableToWrite<Real>("Ep_Production");

		//particles_->registerSharedVariable(ep_dissipation_, "Ep_Dissipation_");
		particles_->addVariableToSort<Real>("Ep_Dissipation_");
		particles_->addVariableToWrite<Real>("Ep_Dissipation_");

		//particles_->registerSharedVariable(ep_diffusion_, "Ep_Diffusion_");
		particles_->addVariableToSort<Real>("Ep_Diffusion_");
		particles_->addVariableToWrite<Real>("Ep_Diffusion_");
	}
	//=================================================================================================//
	void E_TurtbulentModelInner::
		interaction(size_t index_i, Real dt)
	{
		Real rho_i = rho_[index_i];
		//Real turbu_mu_i = turbu_mu_[index_i];
		Real turbu_k_i = turbu_k_[index_i];
		Real turbu_epsilon_i = turbu_epsilon_[index_i];

		Real mu_eff_i = turbu_mu_[index_i] / sigma_E_ + mu_;

		depsilon_dt_[index_i] = 0.0;
		depsilon_dt_without_disspation_[index_i] = 0.0;
		Real epsilon_production(0.0);
		Real epsilon_derivative(0.0);
		Real epsilon_lap(0.0);
		Real epsilon_dissipation(0.0);
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real mu_eff_j = turbu_mu_[index_j] / sigma_E_ + mu_;
			Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
			epsilon_derivative = (turbu_epsilon_i - turbu_epsilon_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
			epsilon_lap += 2.0 * mu_harmo * epsilon_derivative * inner_neighborhood.dW_ij_[n]* this->Vol_[index_j] / rho_i;
		}

		epsilon_production = C_l_ * turbu_epsilon_i * k_production_[index_i] / turbu_k_i;
		epsilon_dissipation = C_2_ * turbu_epsilon_i * turbu_epsilon_i / turbu_k_i;
		
		//** Linearize the source term *
		//epsilon_production = C_l_ * turbu_epsilon_prior_[index_i] * k_production_prior_[index_i] / turbu_k_prior_[index_i];
		//epsilon_dissipation = (C_2_ * turbu_epsilon_prior_[index_i] / turbu_k_prior_[index_i]) * turbu_epsilon_i;
		
		depsilon_dt_[index_i] = epsilon_production - epsilon_dissipation + epsilon_lap;
		depsilon_dt_without_disspation_[index_i] = epsilon_production + epsilon_lap;

		//** for test */
		ep_production[index_i] = epsilon_production;
		ep_dissipation_[index_i] = epsilon_dissipation;
		ep_diffusion_[index_i] = epsilon_lap;
	}
	//=================================================================================================//
	void E_TurtbulentModelInner::update(size_t index_i, Real dt)
	{
		//** The near wall epsilon value is updated in wall function part *
		if (is_near_wall_P1_[index_i] != 1)
		{
			turbu_epsilon_[index_i] += depsilon_dt_[index_i] * dt;
			//** If use source term linearisation *
			//Real denominator = 1.0 + C_2_ * turbu_epsilon_[index_i] * dt / turbu_k_[index_i];
			//turbu_epsilon_[index_i] += depsilon_dt_without_disspation_[index_i] * dt;
			//turbu_epsilon_[index_i] /= denominator;
		}
			
	}
//=================================================================================================//
	TKEnergyForce<Inner<>>::TKEnergyForce(BaseInnerRelation& inner_relation)
		: TKEnergyForce<Base, DataDelegateInner>(inner_relation),
		test_k_grad_rslt_(*this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
		B_(*particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
		//tke_acc_inner_(*particles_->registerSharedVariable<Vecd>("TkeAccInner"))
	{
		//this->particles_->registerSharedVariable(tke_acc_inner_, "TkeAccInner");
		//this->particles_->addVariableToWrite<Vecd>("TkeAccInner");
	}
	//=================================================================================================//
    void TKEnergyForce<Inner<>>::interaction(size_t index_i, Real dt)
    {
		Real turbu_k_i = turbu_k_[index_i];
		Vecd force = Vecd::Zero();
		Vecd k_gradient = Vecd::Zero();
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n]* this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
			//** strong form * 
			//k_gradient += -1.0*(turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
			//** weak form * 
			k_gradient += (turbu_k_i + turbu_k_[index_j]) * nablaW_ijV_j;
		    //** If use RKGC *
			//k_gradient += (turbu_k_i * B_[index_j] + turbu_k_[index_j] * B_[index_i]) * nablaW_ijV_j;
			
		}
		force = -1.0 * (2.0 / 3.0) * k_gradient * mass_[index_i];
		force_[index_i] += force;

		//for test
		//tke_acc_inner_[index_i] = force / mass_[index_i];
		test_k_grad_rslt_[index_i] = k_gradient;
    }
	//=================================================================================================//
	TKEnergyForce<Contact<>>::TKEnergyForce(BaseContactRelation& contact_relation)
		: TKEnergyForce<Base, DataDelegateContact>(contact_relation),
		test_k_grad_rslt_(*this->particles_->template getVariableDataByName<Vecd>("TkeGradResult")),
		B_(*particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
	{
		//this->particles_->registerSharedVariable(tke_acc_wall_, "TkeAccWall");
		//this->particles_->addVariableToWrite<Vecd>("TkeAccWall");
	}
	//=================================================================================================//
	void TKEnergyForce<Contact<>>::interaction(size_t index_i, Real dt)
	{
		Real turbu_k_i = turbu_k_[index_i];
		Vecd force = Vecd::Zero();
		Vecd k_gradient = Vecd::Zero();
		for (size_t k = 0; k < DataDelegateContact::contact_configuration_.size(); ++k)
		{
			Neighborhood& contact_neighborhood = (*DataDelegateContact::contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n]* this->Vol_[index_j] * contact_neighborhood.e_ij_[n];
				//** weak form * 
				k_gradient +=  (turbu_k_i + turbu_k_i) * nablaW_ijV_j;
				//** If use RKGC *
				//k_gradient +=  (turbu_k_i + turbu_k_i)* B_[index_i] * nablaW_ijV_j;
			}
		}
		force = -1.0 * (2.0 / 3.0) * k_gradient * mass_[index_i];
		force_[index_i] += force;

		//for test
		//tke_acc_wall_[index_i] = force / mass_[index_i];
		test_k_grad_rslt_[index_i] += k_gradient;
	}
//=================================================================================================//
	TurbuViscousForce<Inner<>>::TurbuViscousForce(BaseInnerRelation& inner_relation)
		: TurbuViscousForce<DataDelegateInner>(inner_relation), 
		ForcePrior(particles_, "ViscousForce"),
		turbu_indicator_(*this->particles_->template getVariableDataByName<int>("TurbulentIndicator")),
		is_extra_viscous_dissipation_(*this->particles_->template getVariableDataByName<int>("TurbulentExtraViscousDissipation")),
		B_(*particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
	{
		//this->particles_->registerSharedVariable(visc_acc_inner_, "ViscousAccInner");
		//this->particles_->addVariableToSort<Vecd>("ViscousAccInner");
		//this->particles_->addVariableToWrite<Vecd>("ViscousAccInner");
	}
//=================================================================================================//
	void TurbuViscousForce<Inner<>>::interaction(size_t index_i, Real dt)
	{
		turbu_indicator_[index_i] = 0;

		Real mu_eff_i = turbu_mu_[index_i] + molecular_viscosity_;
		
		//visc_acc_inner_[index_i] = Vecd::Zero();

		Vecd force = Vecd::Zero();
		Vecd vel_derivative = Vecd::Zero();
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			const Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real mu_eff_j = turbu_mu_[index_j] + molecular_viscosity_;
			Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
			vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
			
			Vecd shear_stress = mu_harmo * vel_derivative;
			Vecd shear_stress_eij = shear_stress.dot(e_ij) * e_ij;
			
			Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
			
			//Real dissipation_judge = 0.0025 * rho_[index_i] * c0_ * smoothing_length_;
			//Real dissipation = 0.1 * rho_[index_i] * c0_ * smoothing_length_;
			
			Real dissipation = rho_[index_i] *  smoothing_length_ * SMIN(3.0 * SMAX(u_jump, Real(0)), c0_);
			Real dissipation_judge = dissipation;


			//** Introduce dissipation *
			Vecd shear_stress_eij_corrected = shear_stress_eij;
			if (mu_harmo < dissipation_judge && is_extra_viscous_dissipation_[index_i] == 1)
			{
				shear_stress_eij_corrected = ((dissipation * vel_derivative).dot(e_ij)) * e_ij;
				turbu_indicator_[index_i]++; //** For test *
			}
			shear_stress = (shear_stress - shear_stress_eij) + shear_stress_eij_corrected;
			
			
			Vecd force_j = 2.0 * mass_[index_i] * shear_stress * inner_neighborhood.dW_ij_[n]* this->Vol_[index_j];
			//Vecd force_j = e_ij.dot((B_[index_i] + B_[index_j]) * e_ij) * mass_[index_i] * shear_stress * inner_neighborhood.dW_ij_[n]* this->Vol_[index_j];
			
			force += force_j;
		}
		viscous_force_[index_i] = force / rho_[index_i];
		//for test
		//visc_acc_inner_[index_i] = force / rho_[index_i]/ mass_[index_i];
	}
//=================================================================================================//
	TurbuViscousForce<Contact<Wall>>::TurbuViscousForce(BaseContactRelation& wall_contact_relation)
		: BaseTurbuViscousForceWithWall(wall_contact_relation),
		wall_particle_spacing_(wall_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
		B_(*particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
	{
		//this->particles_->registerSharedVariable(visc_acc_wall_, "ViscousAccWall");
		//this->particles_->addVariableToSort<Vecd>("ViscousAccWall");
		//this->particles_->addVariableToWrite<Vecd>("ViscousAccWall");
	}
//=================================================================================================//
	void TurbuViscousForce<Contact<Wall>>::interaction(size_t index_i, Real dt)
	{
		//this->visc_acc_wall_[index_i] = Vecd::Zero();
		//int is_near_wall_P2 = this->is_near_wall_P2_[index_i];		
		//** Wall viscous force only affects P2 region fluid particles *
		if (this->is_near_wall_P2_[index_i] != 10)
			return;
		
		Real turbu_k_i = this->turbu_k_[index_i];
		Real turbu_k_i_05 = pow(turbu_k_i, 0.5);
		Real rho_i = this->rho_[index_i];
		const Vecd& vel_i = this->vel_[index_i];

		Real y_p_constant_i = this->y_p_[index_i];

		Vecd force = Vecd::Zero();		
		Vecd e_j_n = Vecd::Zero();
		Vecd e_j_tau = Vecd::Zero();
		Matd WSS_j_tn = Matd::Zero();  //** Wall shear stress of wall particle j on t-n plane *
		Matd WSS_j = Matd::Zero();  //** Wall shear stress of wall particle j on x-y plane *
		Matd Q = Matd::Zero();
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			//StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
			Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
			StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				//Real r_ij = contact_neighborhood.r_ij_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];

				e_j_n = n_k[index_j];
				Q = getTransformationMatrix(e_j_n);

				//** Get tangential unit vector, temporarily only suitable for 2D*
				e_j_tau[0] = e_j_n[1];
				e_j_tau[1] = e_j_n[0] * (-1.0);

				//** Calculate the local friction velocity *
				Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
				
				Real y_p_j = get_distance_from_P_to_wall(y_p_constant_i);
				Real y_star_j = rho_i * C_mu_25_ * turbu_k_i_05 * y_p_j / molecular_viscosity_;
				Real u_star_j = get_dimensionless_velocity(y_star_j);
				Real fric_vel_mag_j = sqrt(C_mu_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

				//** Construct local wall shear stress, if this is on each wall particle j   *
				Real WSS_tn_mag_j = rho_i * fric_vel_mag_j * fric_vel_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));
				
				WSS_j_tn(0, 0) = 0.0;
				WSS_j_tn(0, 1) = WSS_tn_mag_j;
				WSS_j_tn(1, 0) = 0.0;
				WSS_j_tn(1, 1) = 0.0;
				
				//** Transform local wall shear stress to global   *
				WSS_j = Q.transpose() * WSS_j_tn * Q;
				Vecd force_j =  2.0 * mass_[index_i] * WSS_j * e_ij * contact_neighborhood.dW_ij_[n]* this->Vol_[index_j] / rho_i;
				//Vecd force_j = mass_[index_i] * WSS_j * ( (B_[index_i] + B_[index_j])* e_ij ) * contact_neighborhood.dW_ij_[n]* this->Vol_[index_j] / rho_i;

				//force_j = force_j - (force_j.dot(e_j_n)) * e_j_n;

				force += force_j;
			}
		}
		viscous_force_[index_i] += force;
		//** For test *
		//this->visc_acc_wall_[index_i] = force/ mass_[index_i];
	}
//=================================================================================================//
	TurbulentEddyViscosity::
		TurbulentEddyViscosity(SPHBody& sph_body)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body),
		rho_(*particles_->getVariableDataByName<Real>("Density")), 
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation")),
		wall_Y_plus_(*particles_->getVariableDataByName<Real>("WallYplus")),
		wall_Y_star_(*particles_->getVariableDataByName<Real>("WallYstar")),
		mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()) {}
	//=================================================================================================//
	void TurbulentEddyViscosity::update(size_t index_i, Real dt)
	{
		turbu_mu_[index_i] = rho_[index_i] * C_mu_ * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
	}
//=================================================================================================//
	TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL)
		: LocalDynamicsReduce<ReduceMax>(sph_body), DataDelegateSimple(sph_body),
		vel_(*particles_->getVariableDataByName<Vecd>("Velocity")), 
		smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
		speed_ref_turbu_(U_max), advectionCFL_(advectionCFL),
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
	{
		Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
		speed_ref_turbu_ = SMAX(viscous_speed, speed_ref_turbu_);
	}
	//=================================================================================================//
	Real TurbulentAdvectionTimeStepSize::reduce(size_t index_i, Real dt)
	{
		Real turbu_viscous_speed = (fluid_.ReferenceViscosity() + turbu_mu_[index_i])
			/ fluid_.ReferenceDensity() / smoothing_length_min_;
		Real turbu_viscous_speed_squre = turbu_viscous_speed * turbu_viscous_speed;
		Real vel_n_squre = vel_[index_i].squaredNorm();
		Real vel_bigger = SMAX(turbu_viscous_speed_squre, vel_n_squre);

		return vel_bigger;
	}
	//=================================================================================================//
	Real TurbulentAdvectionTimeStepSize::outputResult(Real reduced_value)
	{
		Real speed_max = sqrt(reduced_value);
		return advectionCFL_ * smoothing_length_min_ / (SMAX(speed_max, speed_ref_turbu_) + TinyReal);
	}
//=================================================================================================//
	InflowTurbulentCondition::InflowTurbulentCondition(BodyPartByCell& body_part
		, Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet) :
		BaseFlowBoundaryCondition(body_part), type_turbu_inlet_(type_turbu_inlet),
		relaxation_rate_(relaxation_rate),
		CharacteristicLength_(CharacteristicLength),
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation"))
	{
		TurbulentLength_ = turbulent_length_ratio_for_epsilon_inlet_ * CharacteristicLength_ ;
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
	Real InflowTurbulentCondition::getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k)
	{
		Real u = velocity[0];
		Real temp_in_turbu_k = 1.5 * pow((turbulent_intensity_ * u), 2);
		Real turbu_k_original = turbu_k;
if(type_turbu_inlet_==1)
{		
		Real channel_height = CharacteristicLength_; //** Temporarily treatment *

        //** Impose fully-developed K from PYTHON result */
        //** Calculate the distance to wall, Y. position here is the actual postion in x-y coordinate, no transformation*/
        Real Y = 0.0;
		if(position[1] < channel_height/2.0)
		{
			Y = position[1];
		}
		else if(position[1] > channel_height/2.0)
		{
			Y = channel_height - position[1];
		}

        int polynomial_order = 8 ;
        int num_coefficient = polynomial_order + 1 ;
        //** Coefficient of the polynomia, 8th-order, from py21 dp=0.024 */
        // Real coeff[] = {
        //     1.215679e-02, -6.681989e-02, 5.043783e-01, 
        //     -2.344875e+00,  6.368016e+00, -1.041386e+01, 
        //     1.009652e+01, -5.336236e+00, 1.183368e+00
        // };
		//** Coefficient of the polynomia, 8th-order, from py21 dp=0.1 */
		Real coeff[] = {
            1.159981e-02, -4.662944e-02, 2.837400e-01, 
            -1.193955e+00, 3.034851e+00, -4.766077e+00, 
            4.529136e+00, -2.380854e+00, 5.307586e-01
        };
        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }

        if(Y > channel_height/2.0 || Y < 0.0)
        {
            std::cout<< "position[1]=" <<position[1]<<std::endl;
            std::cout<< "Y=" <<Y<<std::endl;
            std::cout<< "polynomial_value=" <<polynomial_value<<std::endl;
            std::cout<< "Stop" <<std::endl;
            std::cout<< "=================" <<std::endl;
            std::cin.get();
        }
        
        temp_in_turbu_k = polynomial_value;
}
		if (position[0] < 0.0) //** Temporarily treatment *
		{
			turbu_k_original = temp_in_turbu_k;
		}
		return turbu_k_original;
	}
	//=================================================================================================//
	Real InflowTurbulentCondition::getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E)
	{
		//Real temp_in_turbu_E = C_mu_ * pow(turbu_k, 1.5) / (0.1*getTurbulentLength());
		Real temp_in_turbu_E = C_mu_75_ * pow(turbu_k, 1.5) / TurbulentLength_;
		Real turbu_E_original = turbu_E;
if(type_turbu_inlet_==1)
{	
		Real channel_height = CharacteristicLength_; //** Temporarily treatment *

        //** Impose fully-developed K from PYTHON result */
        //** Calculate the distance to wall, Y. position here is the actual postion in x-y coordinate, no transformation*/
        Real Y = 0.0;
		if(position[1] < channel_height/2.0)
		{
			Y = position[1];
		}
		else if(position[1] > channel_height/2.0)
		{
			Y = channel_height - position[1];
		}

        int polynomial_order = 8 ;
        int num_coefficient = polynomial_order + 1 ;
        //** Coefficient of the polynomia, 8th-order, from py21 dp=0.024 */
        // Real coeff[] = {
        //     1.633474e-02,  -2.488756e-01, 1.912092e+00, 
        //     -8.381386e+00,   2.205987e+01, -3.542125e+01, 
        //     3.391904e+01, -1.777442e+01, 3.918818e+00
        // };
		//** Coefficient of the polynomia, 8th-order, from py21 dp=0.1 */
		Real coeff[] = {
            1.428191e-02, -1.766636e-01, 1.153107e+00, 
            -4.515606e+00, 1.103752e+01, -1.694146e+01, 
            1.584534e+01, -8.241577e+00, 1.825421e+00
        };

        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }

        if(Y > channel_height/2.0 || Y < 0.0)
        {
            std::cout<< "position[1]=" <<position[1]<<std::endl;
            std::cout<< "Y=" <<Y<<std::endl;
            std::cout<< "polynomial_value=" <<polynomial_value<<std::endl;
            std::cout<< "Stop" <<std::endl;
            std::cout<< "=================" <<std::endl;
            std::cin.get();
        }
		
		temp_in_turbu_E = polynomial_value;
}
		if (position[0] < 0.0) //** Temporarily treatment *
		{
			turbu_E_original = temp_in_turbu_E;
		}
		return turbu_E_original;
	}
//=================================================================================================//
	JudgeIsNearWall::
		JudgeIsNearWall(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation)
		: LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
		distance_to_dummy_interface_(*particles_->registerSharedVariable<Real>("DistanceToDummyInterface")),
		distance_to_dummy_interface_up_average_(*particles_->registerSharedVariable<Real>("DistanceToDummyInterfaceUpAver")),
		is_near_wall_P1_(*particles_->registerSharedVariable<int>("IsNearWallP1")),
		is_near_wall_P2_(*particles_->registerSharedVariable<int>("IsNearWallP2")),
		index_nearest_(*particles_->registerSharedVariable<int>("NearestIndex")),
		e_nearest_tau_(*particles_->registerSharedVariable<Vecd>("WallNearestTangentialUnitVector")),
		e_nearest_normal_(*particles_->registerSharedVariable<Vecd>("WallNearestNormalUnitVector")),
		pos_(*particles_->getVariableDataByName<Vecd>("Position")),
		dimension_(2),
		fluid_particle_spacing_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
		wall_particle_spacing_(contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
		distance_from_wall_(*particles_->getVariableDataByName<Vecd>("DistanceFromWall"))
	{
		for (size_t k = 0; k != contact_particles_.size(); ++k)
		{
			contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
			contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
		}

		//particles_->registerSharedVariable(distance_to_dummy_interface_up_average_, "DistanceToDummyInterfaceUpAver");
		particles_->addVariableToSort<Real>("DistanceToDummyInterfaceUpAver");
		particles_->addVariableToWrite<Real>("DistanceToDummyInterfaceUpAver");

		//particles_->registerSharedVariable(distance_to_dummy_interface_, "DistanceToDummyInterface");
		particles_->addVariableToSort<Real>("DistanceToDummyInterface");
		particles_->addVariableToWrite<Real>("DistanceToDummyInterface");

		//particles_->registerSharedVariable(index_nearest_, "NearestIndex");
		particles_->addVariableToSort<int>("NearestIndex");
		particles_->addVariableToWrite<int>("NearestIndex");

		//particles_->registerSharedVariable(is_near_wall_P1_, "IsNearWallP1");
		particles_->addVariableToSort<int>("IsNearWallP1");
		particles_->addVariableToWrite<int>("IsNearWallP1");
		
		//particles_->registerSharedVariable(is_near_wall_P2_, "IsNearWallP2");
		particles_->addVariableToSort<int>("IsNearWallP2");
		particles_->addVariableToWrite<int>("IsNearWallP2");

		//particles_->registerSharedVariable(e_nearest_tau_, "WallNearestTangentialUnitVector");
		particles_->addVariableToSort<Vecd>("WallNearestTangentialUnitVector");

		//particles_->registerSharedVariable(e_nearest_normal_, "WallNearestNormalUnitVector");
		particles_->addVariableToSort<Vecd>("WallNearestNormalUnitVector");

		particles_->addVariableToWrite<Vecd>("DistanceFromWall");
	};
	//=================================================================================================//
	void JudgeIsNearWall::interaction(size_t index_i, Real dt)
	{
		//** If not clear the values completely, particles that leave P2 region will still carry the values. *
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

		//** Calculate nearest info. *
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
			StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
			Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
			
			if (contact_neighborhood.current_size_ != 0)
				is_near_contact++;

			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
				Real r_ij = contact_neighborhood.r_ij_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];
				Vecd& n_k_j = n_k[index_j];

				//** The distance to dummy interface is 0.5 dp smaller than the r_ij_normal *  
				r_dummy_normal_j = abs(n_k_j.dot(r_ij * e_ij)) - 0.5 * wall_particle_spacing_;

				/** Get the minimum distance, the distance to wall should not be negative*/
				if (r_ij < r_min && r_dummy_normal_j > 0.0 + TinyReal)
				{
					r_min = r_ij; //** Find the nearest wall particle *
					//**If use level-set,this would not activate.*
					r_dummy_normal = r_dummy_normal_j;
					e_i_nearest_n = n_k[index_j];
					id_nearest_j = index_j;
				}
				//** Only average the bigger value or itself*
				if (r_dummy_normal_j - r_dummy_normal > (-1.0 * TinyReal))
				{
					count_average++;
					//** Sum the projection distances according to the kernel approx. *  
					r_dmy_itfc_n_sum += weight_j * r_dummy_normal_j;
					ttl_weight += weight_j;
				}
			}
		} 
		//** This is a temporary treatment, particles in inlet region is not corrected by wall function *
		//if (is_near_contact > 0 && pos_[index_i][0] > 0.0) 
		if (is_near_contact > 0)
		{
			is_near_wall_P2_[index_i] = 10; //** Particles that have contact are defined as in region P2 *  					
			//** Get the tangential unit vector *
			if (dimension_ == 2)
			{
				e_i_nearest_tau[0] = e_i_nearest_n[1];
				e_i_nearest_tau[1] = e_i_nearest_n[0] * (-1.0);
			}
			//** Check the function *
			if (r_dmy_itfc_n_sum <= 0.0)
			{
				std::cout << "r_dmy_itfc_n_sum is almost zero" << std::endl;
				std::cout << "count=" << count_average << std::endl;
				std::cin.get();
			}
			//** Average the projection distances according to the kernel approx. *  
			distance_to_dummy_interface_up_average_[index_i] = r_dmy_itfc_n_sum / ttl_weight;
			
			//** Store wall-nearest values. *  
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
			//** Choose one kind of the distance to classify *
			Real distance = distance_to_dummy_interface_[index_i];
			 
			//** Classify the wall-nearest paritcles *
			if (distance < 1.0 * fluid_particle_spacing_)
				is_near_wall_P1_[index_i] = 1;
			//** Check the distance. Useless because with contant_y_p this may not make big difference *  
			//if (distance < 0.05 * fluid_particle_spacing_)
			//{
			//	std::cout << "There is a particle too close to wall" << std::endl;
			//	std::cout << "index_i=" << index_i << std::endl;
			//	std::cout << "distance=" << distance << std::endl;
			//	std::cin.get();
			//}
		}
	}
//=================================================================================================//
	StandardWallFunctionCorrection::
		StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation, Real y_p_constant)
		: LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
		y_p_(*particles_->registerSharedVariable<Real>("Y_P")),
		wall_Y_plus_(*particles_->registerSharedVariable<Real>("WallYplus")),
		wall_Y_star_(*particles_->registerSharedVariable<Real>("WallYstar")),
		velo_tan_(*particles_->registerSharedVariable<Real>("TangentialVelocity")),
		velo_friction_(*particles_->registerSharedVariable<Vecd>("FrictionVelocity")),
		vel_(*particles_->getVariableDataByName<Vecd>("Velocity")), rho_(*particles_->getVariableDataByName<Real>("Density")), 
		molecular_viscosity_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()),
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation")),
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
		is_near_wall_P2_(*particles_->getVariableDataByName<int>("IsNearWallP2")),
		velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
		k_production_(*particles_->getVariableDataByName<Real>("K_Production")),
		distance_to_dummy_interface_(*particles_->getVariableDataByName<Real>("DistanceToDummyInterface")),
		distance_to_dummy_interface_up_average_(*particles_->getVariableDataByName<Real>("DistanceToDummyInterfaceUpAver")),
		index_nearest(*particles_->getVariableDataByName<int>("NearestIndex")),
		e_nearest_tau_(*particles_->getVariableDataByName<Vecd>("WallNearestTangentialUnitVector")),
		e_nearest_normal_(*particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector"))
	{
		for (size_t k = 0; k != contact_particles_.size(); ++k)
		{
			contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
			contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
		}

		//particles_->registerSharedVariable(y_p_, "Y_P");
		particles_->addVariableToSort<Real>("Y_P");
		particles_->addVariableToWrite<Real>("Y_P");

		//** Fixed y_p_ as a constant distance *
		std::fill(y_p_.begin(), y_p_.end(), y_p_constant);

		//particles_->registerSharedVariable(wall_Y_plus_, "WallYplus");
		particles_->addVariableToSort<Real>("WallYplus");
		particles_->addVariableToWrite<Real>("WallYplus");

		//** Initial value is important, especially when use log law *
		//particles_->registerSharedVariable(wall_Y_star_, "WallYstar", TinyReal);
		particles_->addVariableToSort<Real>("WallYstar");
		particles_->addVariableToWrite<Real>("WallYstar");

		//particles_->registerSharedVariable(velo_tan_, "TangentialVelocity");
		particles_->addVariableToSort<Real>("TangentialVelocity");
		particles_->addVariableToWrite<Real>("TangentialVelocity");

		//particles_->registerSharedVariable(velo_friction_, "FrictionVelocity");
		particles_->addVariableToSort<Vecd>("FrictionVelocity");
		particles_->addVariableToWrite<Vecd>("FrictionVelocity");
	};
	//=================================================================================================//
	void StandardWallFunctionCorrection::interaction(size_t index_i, Real dt)
	{
		velo_tan_[index_i] = 0.0;
		velo_friction_[index_i] = Vecd::Zero();
		wall_Y_plus_[index_i] = 0.0;
		wall_Y_star_[index_i] = 0.0;

		//** If use level-set to get distance from P to wall, activate this *
		//y_p_[index_i]= distance_to_dummy_interface_levelset_[index_i];

		if (is_near_wall_P2_[index_i] == 10)
		{
			Real y_p_constant_i = y_p_[index_i];

			Real turbu_k_i_05 = pow(turbu_k_[index_i], 0.5);
			Real turbu_k_i_15 = pow(turbu_k_[index_i], 1.5);

			//** Choose one kind of the distance to calculate the wall-nearest values *
			//Real r_dummy_normal = distance_to_dummy_interface_up_average_[index_i];
			//Real r_dummy_normal = distance_to_dummy_interface_[index_i];
			//Real r_dummy_normal = distance_to_dummy_interface_levelset_[index_i];

			//if (r_dummy_normal <= TinyReal)
			//{
				//std::cout << "r_dummy_normal <= TinyReal" << std::endl;
				//std::cin.get();
			//}
			Vecd e_i_nearest_tau = e_nearest_tau_[index_i];
			Vecd e_i_nearest_n = e_nearest_normal_[index_i];
			const Vecd& vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];
			Real nu_i = molecular_viscosity_ / rho_i;

			//** Calcualte Y_star, note the current code is based on Y_star *
			wall_Y_star_[index_i] = y_p_constant_i * C_mu_25_ * turbu_k_i_05  / nu_i;
			
			//** Calculate friction velocity, including P2 region. *  
			Real velo_fric_mag = 0.0;
			Real velo_tan_mag = 0.0; //** tangitial velo magnitude for fluid particle i *

			velo_tan_mag = abs(e_i_nearest_tau.dot(vel_i));
			velo_tan_[index_i] = velo_tan_mag;

			if (wall_Y_star_[index_i] != static_cast<Real>(wall_Y_star_[index_i]))
			{
				std::cout << "y* is not a real value, please check" << std::endl;
				std::cin.get();
			}
				

			Real u_star = get_dimensionless_velocity(wall_Y_star_[index_i]);
			velo_fric_mag = sqrt(C_mu_25_ * turbu_k_i_05 * velo_tan_mag / u_star);

			if (velo_fric_mag != static_cast<Real>(velo_fric_mag))
			{
				std::cout << "friction velocity is not a real, please check" << std::endl;
				std::cout << "velo_fric=" << velo_fric_mag << std::endl << "velo_tan_mag=" << velo_tan_mag << std::endl;
				std::cout << "turbu_k_=" << pow(turbu_k_[index_i], 0.5) << std::endl;
				std::cout << "sum=" << (Karman_ * velo_tan_mag * C_mu_25_ * pow(turbu_k_[index_i], 0.5) /
					log(turbu_const_E_ * C_mu_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_)) << std::endl;
				std::cout << "numerator=" << Karman_ * velo_tan_mag * C_mu_25_ * pow(turbu_k_[index_i], 0.5) << std::endl;
				std::cout << "denominator=" << log(turbu_const_E_ * C_mu_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_) << std::endl;
				Real temp = C_mu_25_* pow(turbu_k_[index_i], 0.5)* velo_tan_mag / u_star;

				std::cout << "temp =" <<temp<< std::endl;

				std::cout << "pow(turbu_k_[index_i], 0.5) =" << pow(turbu_k_[index_i], 0.5) << std::endl;
				std::cout << "velo_tan_mag / u_star =" << velo_tan_mag / u_star << std::endl;
				std::cout << "velo_tan_mag =" << velo_tan_mag << std::endl;
				std::cout << " u_star =" << u_star << std::endl;
				std::cin.get();
			}

			//** friction velocity have the same direction of vel_i, if not, change its direction *
			velo_friction_[index_i] = velo_fric_mag * e_i_nearest_tau;
			if (vel_i.dot(velo_friction_[index_i]) < 0.0)
				velo_friction_[index_i] = -1.0 * velo_friction_[index_i];

			//** Calcualte Y_plus  *
			wall_Y_plus_[index_i] = y_p_constant_i * velo_fric_mag  / nu_i;
			
			// ** Correct the near wall values, only for P1 region *
			if (is_near_wall_P1_[index_i] == 1)
			{
				Matd vel_grad_i_tn = Matd::Zero();  //** velocity gradient of wall-nearest fluid particle i on t-n plane *
				Matd Q = Matd::Zero();
				Real total_weight = 0.0;
				
				Real epsilon_p_weighted_sum = 0.0;
				Real dudn_p_weighted_sum = 0.0;
				Real G_k_p_weighted_sum = 0.0;
				
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
					Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];

					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						Real epsilon_p_j = 0.0;
						Real dudn_p_j = 0.0;
						Real G_k_p_j = 0.0;
						
						Real y_p_j = 0.0;

						Vecd e_j_tau = Vecd::Zero();

						size_t index_j = contact_neighborhood.j_[n];
						//Real r_ij = contact_neighborhood.r_ij_[n];
						//Vecd& e_ij = contact_neighborhood.e_ij_[n];
						Vecd e_j_n  = n_k[index_j];
						
						//** Get tangential unit vector, temporarily only suitable for 2D*
						e_j_tau[0] = e_j_n[1];
						e_j_tau[1] = e_j_n[0] * (-1.0);

						y_p_j = get_distance_from_P_to_wall(y_p_constant_i);

						Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
						total_weight += weight_j;

						Real denominator_log_law_j = C_mu_25_ * turbu_k_i_05 * Karman_ * y_p_j;

						Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
						Real y_star_j =  C_mu_25_ * turbu_k_i_05 * y_p_j / nu_i;
						Real u_star_j = get_dimensionless_velocity(y_star_j);
						Real fric_vel_mag_j = sqrt(C_mu_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

						Real dudn_p_mag_j = get_near_wall_velocity_gradient_magnitude(y_star_j, fric_vel_mag_j, denominator_log_law_j, nu_i);
						dudn_p_j =  dudn_p_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));
						dudn_p_weighted_sum += weight_j * dudn_p_j;

						if (y_star_j < y_star_threshold_laminar_ && GlobalStaticVariables::physical_time_ > start_time_laminar_)
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

				//** Correct normal velocity at particle P *
				//vel_[index_i] = vel_i - (vel_i.dot(e_i_nearest_n)) * e_i_nearest_n;
			}
			//** Correct normal velocity at particle P *
			//vel_[index_i] = vel_i - (vel_i.dot(e_i_nearest_n)) * e_i_nearest_n;
		}
	}
//=================================================================================================//
	ConstrainNormalVelocityInRegionP::
		ConstrainNormalVelocityInRegionP(SPHBody& sph_body)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body), 
		vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
		is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
		e_nearest_normal_(*particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")){}
	//=================================================================================================//
	void ConstrainNormalVelocityInRegionP::update(size_t index_i, Real dt)
	{
		if (is_near_wall_P1_[index_i] == 1)
		{
			vel_[index_i] = vel_[index_i] - (vel_[index_i].dot(e_nearest_normal_[index_i])) * e_nearest_normal_[index_i];
		}
	}
//=================================================================================================//

//=================================================================================================//
	ConstrainVelocityAt_Y_Direction::
		ConstrainVelocityAt_Y_Direction(SPHBody& sph_body, Real Length_channel)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body), 
		vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
		pos_(*particles_->getVariableDataByName<Vecd>("Position")),
		length_channel_(Length_channel){}
	//=================================================================================================//
	void ConstrainVelocityAt_Y_Direction::update(size_t index_i, Real dt)
	{
		if(pos_[index_i][0]> 0.5 * length_channel_) //** Very temporary treatment *
		{
			vel_[index_i][1] = 0.0;
		}
		
	}
//=================================================================================================//
	UpdateTurbulentPlugFlowIndicator::
		UpdateTurbulentPlugFlowIndicator(SPHBody& sph_body, Real DH)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body), 
		turbu_plug_flow_indicator_(*particles_->registerSharedVariable<int>("TurbulentPlugFlowIndicator")),
		pos_(*particles_->getVariableDataByName<Vecd>("Position")), channel_width_(DH)
		{
			//particles_->registerSharedVariable(turbu_plug_flow_indicator_, "TurbulentPlugFlowIndicator");
			particles_->addVariableToSort<int>("TurbulentPlugFlowIndicator");
			particles_->addVariableToWrite<int>("TurbulentPlugFlowIndicator");
		}
	//=================================================================================================//
	void UpdateTurbulentPlugFlowIndicator::update(size_t index_i, Real dt)
	{
		turbu_plug_flow_indicator_[index_i] = 0 ;
		if(pos_[index_i][0]> 0.0 ) //** Buffer region is still applied tran.vel. Very temporary treatment *
		{
			if(pos_[index_i][1]> 0.25 * channel_width_ && pos_[index_i][1] < 0.75 * channel_width_ ) //** Very temporary treatment * 
			{
				turbu_plug_flow_indicator_[index_i] = 1 ;
			}
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
	GetLimiterOfTransportVelocityCorrection::
		GetLimiterOfTransportVelocityCorrection(SPHBody& sph_body, Real slope)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body),
		h_ref_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
		zero_gradient_residue_(*particles_->getVariableDataByName<Vecd>("ZeroGradientResidue")),
		slope_(slope),
		limiter_tvc_(*particles_->registerSharedVariable<Real>("LimiterOfTVC"))
		{
			particles_->addVariableToWrite<Real>("LimiterOfTVC");
		}
	//=================================================================================================//
	void GetLimiterOfTransportVelocityCorrection::update(size_t index_i, Real dt)
	{
		Real squared_norm = zero_gradient_residue_[index_i].squaredNorm();
		limiter_tvc_[index_i] = SMIN(slope_ * squared_norm * h_ref_ * h_ref_, Real(1) ) ;
	}
//=================================================================================================//
//*********************TESTING MODULES*********************
//=================================================================================================//

//=================================================================================================//
	BaseGetTimeAverageData::BaseGetTimeAverageData(BaseInnerRelation& inner_relation, int num_observer_points)
		: BaseTurtbulentModel<Base, DataDelegateInner>(inner_relation),  plt_engine_(),
		pos_(*particles_->getVariableDataByName<Vecd>("Position")), 
		turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
		turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
		turbu_epsilon_(*particles_->getVariableDataByName<Real>("TurbulentDissipation")),
		num_cell(num_observer_points)
	{
		num_data = 5;
		file_name_.push_back("vel_x_sto_");
		file_name_.push_back("turbu_k_sto_");
		file_name_.push_back("turbu_epsilon_sto_");
		file_name_.push_back("turbu_mu_sto_");
		file_name_.push_back("vel_sto_");

		num_in_cell_.resize(num_cell);
		data_time_aver_sto_.resize(num_cell); //Rows

		data_sto_.resize(num_cell); //Rows
		for (int i = 0; i != num_cell; ++i)
		{
			data_sto_[i].resize(num_data); //Cols
		}

		for (int j = 0; j != num_data; ++j)
		{
			file_path_output_ = "../bin/output/" + file_name_[j] + ".dat";
			std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
			out_file << "run_time" << "   ";
			for (int i = 0; i != num_cell; ++i)
			{
				std::string quantity_name_i = file_name_[j] + "[" + std::to_string(i) + "]";
				plt_engine_.writeAQuantityHeader(out_file, data_sto_[i][j], quantity_name_i);
			}
			out_file << "\n";
			out_file.close();
		}
	}
	//=================================================================================================//
	void BaseGetTimeAverageData::output_time_history_data(Real cutoff_time)
	{
		/** Output for .dat file. */
		for (int j = 0; j != num_data; ++j)
		{
			file_path_output_ = "../bin/output/" + file_name_[j] + ".dat";
			std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
			out_file << GlobalStaticVariables::physical_time_ << "   ";
			for (int i = 0; i != num_cell; ++i)
			{
				//if (num_in_cell_[i] == 0 && GlobalStaticVariables::physical_time_ > cutoff_time)
				//{
				//	std::cout << "There is a empaty monitoring cell, cell number=" << i << std::endl;
				//	std::cin.get();
				//}
				num_in_cell_[i] == 0 ? plt_engine_.writeAQuantity(out_file, 0.0) :
					plt_engine_.writeAQuantity(out_file, data_sto_[i][j] / num_in_cell_[i]);
			}
			out_file << "\n";
			out_file.close();
		}
		//** Clear data *
		for (int i = 0; i < num_cell; i++)
		{
			num_in_cell_[i] = 0;
			for (int j = 0; j != num_data; ++j)
			{
				data_sto_[i][j] = 0.0;
			}
		}
	}
	//=================================================================================================//
	void BaseGetTimeAverageData::get_time_average_data(Real cutoff_time)
	{
		for (int j = 0; j != num_data; ++j)
		{
			data_loaded_.clear();
			int num_line_data = 0;
			//** Load data *
			file_path_input_ = "../bin/output/" + file_name_[j] + ".dat";
			std::ifstream in_file(file_path_input_.c_str());
			bool skipFirstLine = true;
			std::string line;
			while (std::getline(in_file, line))
			{
				if (skipFirstLine)
				{
					skipFirstLine = false;
					continue;
				}
				num_line_data++;
				std::vector<Real> data_point;
				std::istringstream iss(line);
				Real value;
				while (iss >> value)
				{
					data_point.push_back(value);
				}
				data_loaded_.push_back(data_point);
			}

			in_file.close();
			//** Deal with data *
			for (int k = 0; k != num_cell; ++k)
			{
				Real sum = 0.0;
				int count = 0;
				for (int i = 0; i != num_line_data; ++i)
				{
					if (data_loaded_[i][0] > cutoff_time)
					{
						count++;
						Real delta_t = data_loaded_[i][0] - data_loaded_[i - 1][0];
						sum += data_loaded_[i][k + 1] * delta_t;
						//sum += data_loaded_[i][k + 1]; //**the first col is time*
					}
				}
				//data_time_aver_sto_[k] = sum / count;
				data_time_aver_sto_[k] = sum / (data_loaded_[num_line_data - 1][0] - cutoff_time);
			}
			//** Output data *
			file_path_output_ = "../bin/output/TimeAverageData.dat";
			std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
			out_file << file_name_[j] << "\n";
			for (int k = 0; k != num_cell; ++k)
			{
				plt_engine_.writeAQuantity(out_file, data_time_aver_sto_[k]);
			}
			out_file << "\n";
			out_file.close();
		}
		std::cout << "The cutoff_time is " << cutoff_time << std::endl;
	}
	//=================================================================================================//
	GetTimeAverageCrossSectionData::GetTimeAverageCrossSectionData(BaseInnerRelation& inner_relation, int num_observer_points, const StdVec<Real>& bound_x, Real offset_dist_y)
		: BaseGetTimeAverageData(inner_relation, num_observer_points)
	{
		x_min_ = bound_x[0];
		x_max_ = bound_x[1];
		offset_dist_y_ = offset_dist_y;
		//** Get the center coordinate of the monitoring cell *
		for (int i = 0; i < num_cell; i++)
		{
			Real upper_bound = ((i + 1) * particle_spacing_min_ + offset_dist_y_);
			Real lower_bound = (i * particle_spacing_min_ + offset_dist_y_);
			monitor_cellcenter_y.push_back((lower_bound + upper_bound) / 2.0);
		}
		file_path_output_ = "../bin/output/monitor_cell_center_y.dat";
		std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
		for (int i = 0; i != num_cell; ++i)
		{
			plt_engine_.writeAQuantity(out_file, monitor_cellcenter_y[i]);
			out_file << "\n";
		}
		out_file << "\n";
		out_file.close();
	}
	//=================================================================================================//
	void GetTimeAverageCrossSectionData::update(size_t index_i, Real dt)
	{
		//** Get data *
		if (pos_[index_i][0] > x_min_ && pos_[index_i][0] <= x_max_)
		{
			for (int i = 0; i < num_cell; i++)
			{
				if (pos_[index_i][1] > (i * particle_spacing_min_ + offset_dist_y_) &&
					pos_[index_i][1] <= ((i + 1) * particle_spacing_min_ + offset_dist_y_))
				{
					num_in_cell_[i] += 1;
					data_sto_[i][0] += vel_[index_i][0];
					data_sto_[i][1] += turbu_k_[index_i];
					data_sto_[i][2] += turbu_epsilon_[index_i];
					data_sto_[i][3] += turbu_mu_[index_i];
					data_sto_[i][4] += vel_[index_i].norm();
				}
			}
		}
	}
	//=================================================================================================//
	GetTimeAverageCrossSectionData_Y::GetTimeAverageCrossSectionData_Y(BaseInnerRelation& inner_relation,
		int num_observer_points, Real observe_x_ratio, const StdVec<Real>& bound_y, const StdVec<Real>& bound_x_f, const StdVec<Real>& bound_x_b)
		: BaseGetTimeAverageData(inner_relation, num_observer_points), observe_x_ratio_(observe_x_ratio),
		bound_x_f_(bound_x_f), bound_x_b_(bound_x_b), bound_y_(bound_y)
	{
		observe_x_spacing_ = particle_spacing_min_ * observe_x_ratio_;
	}
	//=================================================================================================//
	void GetTimeAverageCrossSectionData_Y::update(size_t index_i, Real dt)
	{
		//** Get data *
		if (pos_[index_i][1] > bound_y_[0] && pos_[index_i][1] <= bound_y_[1])
		{
			for (int i = 0; i < num_cell; i++)
			{
				if (i < static_cast<int>(bound_x_f_.size()) - 1) //* Front of cylinder
				{
					if (pos_[index_i][0] > bound_x_f_[i] && pos_[index_i][0] <= bound_x_f_[i + 1])
					{
						num_in_cell_[i] += 1;
						data_sto_[i][0] += vel_[index_i][0];
						data_sto_[i][1] += turbu_k_[index_i];
						data_sto_[i][2] += turbu_epsilon_[index_i];
						data_sto_[i][3] += turbu_mu_[index_i];
						data_sto_[i][4] += vel_[index_i].norm();
					}
				}
				else if (i >= static_cast<int>(bound_x_f_.size()) - 1) //* behind of cylinder
				{
					int j = i - (static_cast<int>(bound_x_f_.size()) - 1);
					if (pos_[index_i][0] > bound_x_b_[j] && pos_[index_i][0] <= bound_x_b_[j + 1])
					{
						num_in_cell_[i] += 1;
						data_sto_[i][0] += vel_[index_i][0];
						data_sto_[i][1] += turbu_k_[index_i];
						data_sto_[i][2] += turbu_epsilon_[index_i];
						data_sto_[i][3] += turbu_mu_[index_i];
						data_sto_[i][4] += vel_[index_i].norm();
					}
				}
			}

		}
	}
	//=================================================================================================//
	void GetTimeAverageCrossSectionData_Y::output_monitor_x_coordinate()
	{
		StdVec<Real> monitor_cellcenter_x;
		if (static_cast<int>(bound_x_f_.size()) != 0)
		{
			for (int i = 0; i < static_cast<int>(bound_x_f_.size()) - 1; i++)
			{
				monitor_cellcenter_x.push_back((bound_x_f_[i] + bound_x_f_[i + 1]) / 2.0);
			}
		}
		if (static_cast<int>(bound_x_b_.size()) != 0)
		{
			for (int i = 0; i < static_cast<int>(bound_x_b_.size()) - 1; i++)
			{
				monitor_cellcenter_x.push_back((bound_x_b_[i] + bound_x_b_[i + 1]) / 2.0);
			}
		}
		

		file_path_output_ = "../bin/output/monitor_cell_center_x.dat";
		std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
		for (int i = 0; i != num_cell; ++i)
		{
			plt_engine_.writeAQuantity(out_file, monitor_cellcenter_x[i]);
			out_file << "\n";
		}
		out_file << "\n";
		out_file.close();
	}
//=================================================================================================//
	ClearYPositionForTest::
		ClearYPositionForTest(SPHBody& sph_body)
		: LocalDynamics(sph_body), DataDelegateSimple(sph_body),
		pos_(*particles_->getVariableDataByName<Vecd>("Position")), vel_(*particles_->getVariableDataByName<Vecd>("Velocity")) {}
	//=================================================================================================//
	void ClearYPositionForTest::update(size_t index_i, Real dt)
	{
		vel_[index_i][1] = 0.0;
	}
//=================================================================================================//
}
//=================================================================================================//
}
//=================================================================================================//