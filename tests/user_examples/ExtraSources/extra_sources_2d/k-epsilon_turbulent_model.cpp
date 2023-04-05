#include "k-epsilon_turbulent_model.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbulentClosureCoefficient::BaseTurbulentClosureCoefficient()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(1.0e-2), sigma_k(1.0),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793) {}
		//=================================================================================================//
		BaseTurtbulentData::BaseTurtbulentData(BaseInnerRelation& inner_relation)
			: FluidDataInner(inner_relation), BaseTurbulentClosureCoefficient(),
			particle_spacing_min_(inner_relation.real_body_->sph_adaptation_->MinimumSpacing()),
			Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			vel_(particles_->vel_),acc_prior_(particles_->acc_prior_),
			mu_(particles_->fluid_.ReferenceViscosity()),dimension_(Vecd(0).size()) {}
		//=================================================================================================//
		BaseTurtbulentModel::BaseTurtbulentModel(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()),
			BaseTurtbulentData(inner_relation),
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
			grad_calculated_(*particles_->getVariableByName<Matd>("GradCalculated")) {}
		//=================================================================================================//
		K_TurtbulentModelInner::K_TurtbulentModelInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModel(inner_relation)
		{
			particles_->registerVariable(dk_dt_, "ChangeRateOfTKE");
			particles_->registerSortableVariable<Real>("ChangeRateOfTKE");

			//particles_->registerAVariable(turbu_k_, "TurbulenceKineticEnergy",6.43e-9);
			//particles_->registerAVariable(turbu_k_, "TurbulenceKineticEnergy", 0.00015);
			particles_->registerVariable(turbu_k_, "TurbulenceKineticEnergy", 0.000180001);
			particles_->registerSortableVariable<Real>("TurbulenceKineticEnergy");
			particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

			particles_->registerVariable(turbu_mu_, "TurbulentViscosity", 1.0e-9);
			particles_->registerSortableVariable<Real>("TurbulentViscosity");
			particles_->addVariableToWrite<Real>("TurbulentViscosity");


			//particles_->registerAVariable(turbu_epsilon_, "TurbulentDissipation", 3.72e-9);
			//particles_->registerAVariable(turbu_epsilon_, "TurbulentDissipation", 2.156208e-5);
			particles_->registerVariable(turbu_epsilon_, "TurbulentDissipation", 3.326679e-5);
			particles_->registerSortableVariable<Real>("TurbulentDissipation");
			particles_->addVariableToWrite<Real>("TurbulentDissipation");

			particles_->registerVariable(grad_vel_ij, "VelocityGradient");
			particles_->registerSortableVariable<Matd>("VelocityGradient");

			particles_->registerVariable(Rij, "RenoyldShearStress");
			particles_->registerSortableVariable<Matd>("RenoyldShearStress");

			particles_->registerVariable(transpose_grad_vel_ij, "TransposeVelocityGradient");
			particles_->registerSortableVariable<Matd>("TransposeVelocityGradient");


			particles_->registerVariable(production_k_, "Production_K");
			particles_->registerSortableVariable<Real>("Production_K");
			particles_->addVariableToWrite<Real>("Production_K");

			particles_->registerVariable(vel_x_n_, "Velocity_X");
			particles_->registerSortableVariable<Real>("Velocity_X");

			particles_->registerVariable(lap_k_, "Lap_K");
			particles_->registerSortableVariable<Real>("Lap_K");
			particles_->addVariableToWrite<Real>("Lap_K");

			particles_->registerVariable(lap_k_term_, "Lap_K_term");
			particles_->registerSortableVariable<Real>("Lap_K_term");
			particles_->addVariableToWrite<Real>("Lap_K_term");

			particles_->registerVariable(temp_dW_, "dW_term");
			particles_->registerSortableVariable<Real>("dW_term");
			particles_->addVariableToWrite<Real>("dW_term");

			particles_->registerVariable(temp_rij_, "rij_term");
			particles_->registerSortableVariable<Real>("rij_term");
			particles_->addVariableToWrite<Real>("rij_term");

			particles_->registerSortableVariable<int>("SurfaceIndicator");
		}
		//=================================================================================================//
		void K_TurtbulentModelInner::interaction(size_t index_i, Real dt)
		{
			Vecd vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];
			Real turbu_mu_i = turbu_mu_[index_i];
			Real turbu_k_i = turbu_k_[index_i];

			dk_dt_[index_i] = 0.0;
			Real k_production(0.0);
			Real k_derivative(0.0);
			Real k_lap(0.0);
			Matd strain_rate = Matd::Zero();
			Matd Re_stress = Matd::Zero();
			Matd velocity_gradient = Matd::Zero();
			
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				velocity_gradient += (vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
				
				k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				k_lap += 2.0 * (mu_+turbu_mu_i/sigma_k)*k_derivative * inner_neighborhood.dW_ijV_j_[n]/ rho_i;
			}
			strain_rate = 0.5 * (velocity_gradient.transpose() + velocity_gradient);
			Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
			Matd k_production_matrix = Re_stress * velocity_gradient.transpose();
			k_production = k_production_matrix.sum();

			dk_dt_[index_i] = k_production +turbu_epsilon_[index_i] + k_lap;
		}
		//=================================================================================================//
		void K_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_k_[index_i] += dk_dt_[index_i] * dt;
		}

	}
	//=================================================================================================//
}
//=================================================================================================//