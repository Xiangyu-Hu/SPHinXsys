#pragma once
#include "k-epsilon_turbulent_model.hpp"
namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
	BaseTurbuClosureCoeff::BaseTurbuClosureCoeff()
		: Karman(0.4187), C_mu(0.09), TurbulentIntensity(5.0e-2), sigma_k(1.0),
		C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793) {}
//=================================================================================================//
    void GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
    {
		//** The near wall velo grad is updated in wall function part *
		if (is_near_wall_P1_[index_i] == 1)
		{
			return;
		}
		Vecd vel_i = vel_[index_i];
		velocity_gradient_[index_i] = Matd::Zero();
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			//** Strong form *
			velocity_gradient_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
			//** Weak form *
			//velocity_gradient_[index_i] += (vel_i + vel_[index_j]) * nablaW_ijV_j.transpose();
		}
    }
//=================================================================================================//
	TKEnergyAcc<Inner<>>::TKEnergyAcc(BaseInnerRelation& inner_relation)
		: TKEnergyAcc<Base, FluidDataInner>(inner_relation),
		test_k_grad_rslt_(*this->particles_->template getVariableByName<Vecd>("TkeGradResult"))
	{
		this->particles_->registerVariable(tke_acc_inner_, "TkeAccInner");
		this->particles_->addVariableToWrite<Vecd>("TkeAccInner");
	}

//=================================================================================================//
    void TKEnergyAcc<Inner<>>::interaction(size_t index_i, Real dt)
    {
		Real turbu_k_i = turbu_k_[index_i];
		Vecd acceleration = Vecd::Zero();
		Vecd k_gradient = Vecd::Zero();
		//Vecd nablaW_ijV_j_test= Vecd::Zero();
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			//** strong form * 
			//k_gradient += -1.0*(turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
			//** weak form * 
			k_gradient += -1.0 * (-1.0) * (turbu_k_i + turbu_k_[index_j]) * nablaW_ijV_j;

			//nablaW_ijV_j_test = nablaW_ijV_j;

			//** For test *
			//if (GlobalStaticVariables::physical_time_ > 1.0 && pos_[index_i][1] > 0.95)
			//{
			//	std::cout << "index_i=" << index_i << std::endl;
			//	std::cout << "index_j=" << index_j << std::endl;
			//	std::cout << "nablaW_ijV_j=" << nablaW_ijV_j_test << std::endl;
			//}
		}
		acceleration = -1.0 * (2.0 / 3.0) * k_gradient;


		//** For test *
		//acceleration[0] = 0.0;
		
		acc_prior_[index_i] += acceleration;

		//for test
		tke_acc_inner_[index_i] = acceleration;
		test_k_grad_rslt_[index_i] = k_gradient;
    }
//=================================================================================================//
	TKEnergyAcc<Contact<>>::TKEnergyAcc(BaseContactRelation& contact_relation)
		: TKEnergyAcc<Base, FluidContactData>(contact_relation),
		test_k_grad_rslt_(*this->particles_->template getVariableByName<Vecd>("TkeGradResult"))
	{
		this->particles_->registerVariable(tke_acc_wall_, "TkeAccWall");
		this->particles_->addVariableToWrite<Vecd>("TkeAccWall");
	}
//=================================================================================================//
	void TKEnergyAcc<Contact<>>::interaction(size_t index_i, Real dt)
	{
		Real turbu_k_i = turbu_k_[index_i];
		Vecd acceleration = Vecd::Zero();
		Vecd k_gradient = Vecd::Zero();
		for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
		{
			Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
				//** weak form * 
				k_gradient += -1.0 * (-1.0) * (turbu_k_i + turbu_k_i) * nablaW_ijV_j;
			}
		}
		acceleration = -1.0 * (2.0 / 3.0) * k_gradient;

		//** For test *
		//acceleration[0] = 0.0;

		acc_prior_[index_i] += acceleration;

		//for test
		tke_acc_wall_[index_i] = acceleration;
		test_k_grad_rslt_[index_i] += k_gradient;
	}
//=================================================================================================//
	TurbuViscousAcceleration<Inner<>>::TurbuViscousAcceleration(BaseInnerRelation& inner_relation)
		: TurbuViscousAcceleration<FluidDataInner>(inner_relation)
	{
		this->particles_->registerVariable(visc_acc_inner_, "ViscousAccInner");
		this->particles_->addVariableToWrite<Vecd>("ViscousAccInner");
		//this->particles_->registerVariable(visc_acc_wall_, "ViscousAccWall");
		//this->particles_->addVariableToWrite<Vecd>("ViscousAccWall");
	}
//=================================================================================================//
	void TurbuViscousAcceleration<Inner<>>::interaction(size_t index_i, Real dt)
	{
		Real mu_eff_i = turbu_mu_[index_i] + mu_;
		Vecd acceleration = Vecd::Zero();
		Vecd vel_derivative = Vecd::Zero();
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			const Vecd& e_ij = inner_neighborhood.e_ij_[n];

			Real mu_eff_j = turbu_mu_[index_j] + mu_;
			Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
			vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);

			Vecd acc_j = 2.0 * mu_harmo * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
			acceleration += acc_j;
		}

		//** For test *
		//acceleration[1] = 0.0;

		acc_prior_[index_i] += acceleration / rho_[index_i];
		//for test
		visc_acc_inner_[index_i] = acceleration / rho_[index_i];
	}
//=================================================================================================//
	TurbuViscousAcceleration<ContactWall<>>::TurbuViscousAcceleration(BaseContactRelation& wall_contact_relation)
		: BaseTurbuViscousAccelerationWithWall(wall_contact_relation)
	{
		this->particles_->registerVariable(visc_acc_wall_, "ViscousAccWall");
		this->particles_->addVariableToWrite<Vecd>("ViscousAccWall");
		//this->particles_->registerVariable(visc_direction_matrix_, "ViscDirectionMatrix");
		//this->particles_->addVariableToWrite<Vecd>("ViscDirectionMatrix");
	}
//=================================================================================================//

	void TurbuViscousAcceleration<ContactWall<>>::interaction(size_t index_i, Real dt)
	{
		Real turbu_mu_i = this->turbu_mu_[index_i];
		Real rho_i = this->rho_[index_i];
		const Vecd& vel_i = this->vel_[index_i];
		const Vecd& vel_fric_i = this->velo_friction_[index_i];

		Vecd e_x = { 1.0, 0.0 };
		Vecd e_y = { 0.0, 1.0 };
		Real vel_fric_i_x = vel_fric_i.dot(e_x);
		Real vel_fric_i_y = vel_fric_i.dot(e_y);

		Vecd e_tau = vel_fric_i.normalized();
		Real y_plus_i = this->wall_Y_plus_[index_i];
		Real y_p = this->y_p_[index_i];

		Real u_plus_i = 0.0;
		Real mu_w = 0.0;
		Real mu_p = 0.0;
		Real theta = 0.0;

		Vecd acceleration = Vecd::Zero();
		Vecd vel_derivative = Vecd::Zero();
		Matd direc_matrix = Matd::Zero();
		
		Vecd e_j_n = Vecd::Zero();
		Vecd e_j_tau = Vecd::Zero();
		Matd WSS_j_tn = Matd::Zero();  //** Wall shear stress of wall particle j on t-n plane *
		Matd WSS_j = Matd::Zero();  //** Wall shear stress of wall particle j on x-y plane *
		Matd Q = Matd::Zero();
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
			Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
			StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);

			
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real r_ij = contact_neighborhood.r_ij_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];

				//Real e_ij_x = e_ij.dot(e_x);
				//Real e_ij_y = e_ij.dot(e_y);
				
				e_j_n = n_k[index_j];
				Q = getTransformationMatrix(e_j_n);

				//** Get tangential unit vector, temporarily only suitable for 2D*
				e_j_tau[0] = e_j_n[1];
				e_j_tau[1] = e_j_n[0] * (-1.0);
				//if (vel_i.dot(e_j_tau) < 0.0)
					//e_j_tau = -1.0 * e_j_tau;  //** Assume the tangential unit vector has the same direction of velocity *

				//** Construct local wall shear stress, if this is on each wall particle j   *
				WSS_j_tn(0, 0) = 0.0;
				WSS_j_tn(0, 1) = rho_i* vel_fric_i.dot(vel_fric_i)* vel_i.dot(e_j_tau)/ (vel_i.norm() + 0.01 * smoothing_length_);
				WSS_j_tn(1, 0) = 0.0;
				WSS_j_tn(1, 1) = rho_i * vel_fric_i.dot(vel_fric_i) * vel_i.dot(e_j_n) / (vel_i.norm() + 0.01 * smoothing_length_);;
				
				//** Transform local wall shear stress to global   *
				WSS_j = Q.transpose() * WSS_j_tn * Q;
				Vecd acc_j = -1.0 * -1.0 * 2.0 * WSS_j * e_ij * contact_neighborhood.dW_ijV_j_[n]/rho_i;
				//std::cout << "acc_j=" << acc_j << std::endl;
				//std::cout << "vel_i=" << vel_i << std::endl;
				//std::cout << "vel_fric_i=" << vel_fric_i << std::endl;
				
				//** Calculate the direction matrix of wall shear stress *
				//direc_matrix = e_tau * e_j_n.transpose() + (e_tau * e_j_n.transpose()).transpose();
				//Vecd acc_j = -1.0 * -1.0 * 2.0 * vel_fric_i.dot(vel_fric_i) * direc_matrix * e_ij * contact_neighborhood.dW_ijV_j_[n];
				
				//Vecd acc_j = Vecd::Zero();
				//acc_j[0] = -1.0 * -1.0 * 2.0 * vel_fric_i_x * vel_fric_i_x  * e_ij_y * contact_neighborhood.dW_ijV_j_[n];
				//acc_j[1] = -1.0 * -1.0 * 2.0 * vel_fric_i_y * vel_fric_i_y * e_ij_x * contact_neighborhood.dW_ijV_j_[n];;

				acceleration += acc_j;
			}
		}
		//** For test *
		//acceleration[1] = 0.0;

		this->acc_prior_[index_i] += acceleration;
		//** For test *
		this->visc_acc_wall_[index_i] = acceleration;
		//this->visc_direction_matrix_[index_i] = direc_matrix;
	}


//=================================================================================================//
}
//=================================================================================================//
}
//=================================================================================================//