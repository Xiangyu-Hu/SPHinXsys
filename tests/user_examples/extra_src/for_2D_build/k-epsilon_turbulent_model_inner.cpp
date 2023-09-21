#include "k-epsilon_turbulent_model_inner.hpp"
#include "k-epsilon_turbulent_model_inner.h"
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbuClosureCoeff::BaseTurbuClosureCoeff()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(5.0e-2), sigma_k(0.8),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793){}
		//=================================================================================================//
		BaseTurtbulentModelInner::BaseTurtbulentModelInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			particle_spacing_min_(inner_relation.real_body_->sph_adaptation_->MinimumSpacing()),
			 rho_(particles_->rho_), vel_(particles_->vel_), 
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), dimension_(Vecd(0).size()),
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		K_TurtbulentModelInner::K_TurtbulentModelInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation)
		{
				particles_->registerVariable(dk_dt_, "ChangeRateOfTKE");
				particles_->registerSortableVariable<Real>("ChangeRateOfTKE");

				particles_->registerVariable(turbu_k_, "TurbulenceKineticEnergy", 0.000180001);
				particles_->registerSortableVariable<Real>("TurbulenceKineticEnergy");
				particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

				particles_->registerVariable(turbu_mu_, "TurbulentViscosity", 1.0e-9);
				particles_->registerSortableVariable<Real>("TurbulentViscosity");
				particles_->addVariableToWrite<Real>("TurbulentViscosity");

				particles_->registerVariable(turbu_epsilon_, "TurbulentDissipation", 3.326679e-5);
				particles_->registerSortableVariable<Real>("TurbulentDissipation");
				particles_->addVariableToWrite<Real>("TurbulentDissipation");

				particles_->registerVariable(k_production_, "K_Production");
				particles_->registerSortableVariable<Real>("K_Production");
				particles_->addVariableToWrite<Real>("K_Production");
				
				particles_->registerVariable(B_, "CorrectionMatrix");
				particles_->registerSortableVariable<Matd>("CorrectionMatrix");
				particles_->addVariableToWrite<Matd>("CorrectionMatrix");

				particles_->registerVariable(is_near_wall_P1_, "IsNearWallP1");
				particles_->registerSortableVariable<int>("IsNearWallP1");
				particles_->addVariableToWrite<int>("IsNearWallP1");


				//** for test */
				particles_->registerVariable(k_diffusion_, "K_Diffusion");
				particles_->registerSortableVariable<Real>("K_Diffusion");
				particles_->addVariableToWrite<Real>("K_Diffusion");

				particles_->addVariableToWrite<Real>("ChangeRateOfTKE");

				particles_->registerVariable(velocity_gradient, "Velocity_Gradient");
				particles_->registerSortableVariable<Matd>("Velocity_Gradient");
				particles_->addVariableToWrite<Matd>("Velocity_Gradient");
				
				//for test
				//particles_->registerVariable(velocity_gradient_wall, "Velocity_Gradient_Wall");
				//particles_->registerSortableVariable<Matd>("Velocity_Gradient_Wall");
				//particles_->addVariableToWrite<Matd>("Velocity_Gradient_Wall");

				particles_->registerVariable(vel_x_, "Velocity_X");
				particles_->registerSortableVariable<Real>("Velocity_X");
				
		}
		//=================================================================================================//
		GetVelocityGradientInner::GetVelocityGradientInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			vel_(particles_->vel_),velocity_gradient(*particles_->getVariableByName<Matd>("Velocity_Gradient")),
			is_near_wall_P1_(*particles_->getVariableByName<int>("IsNearWallP1"))
		{
			//for test
			particles_->registerVariable(velocity_gradient_wall, "Velocity_Gradient_Wall");
			particles_->registerSortableVariable<Matd>("Velocity_Gradient_Wall");
			particles_->addVariableToWrite<Matd>("Velocity_Gradient_Wall");
		}
		//=================================================================================================//
		void K_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_k_[index_i] += dk_dt_[index_i] * dt;
		}
		//=================================================================================================//
		E_TurtbulentModelInner::E_TurtbulentModelInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation), 
			k_production_(*particles_->getVariableByName<Real>("K_Production")),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation"))
		{
			particles_->registerVariable(dE_dt_, "ChangeRateOfTDR");
			particles_->registerSortableVariable<Real>("ChangeRateOfTDR");
			particles_->addVariableToWrite<Real>("ChangeRateOfTDR");

			particles_->registerVariable(ep_production, "Ep_Production");
			particles_->registerSortableVariable<Real>("Ep_Production");
			particles_->addVariableToWrite<Real>("Ep_Production");
			particles_->registerVariable(ep_dissipation_, "Ep_Dissipation_");
			particles_->registerSortableVariable<Real>("Ep_Dissipation_");
			particles_->addVariableToWrite<Real>("Ep_Dissipation_");
			particles_->registerVariable(ep_diffusion_, "Ep_Diffusion_");
			particles_->registerSortableVariable<Real>("Ep_Diffusion_");
			particles_->addVariableToWrite<Real>("Ep_Diffusion_");

		}
		//=================================================================================================//
		void E_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_epsilon_[index_i] += dE_dt_[index_i] * dt;
		}
		//=================================================================================================//
		TKEnergyAccInner::
			TKEnergyAccInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation), acc_prior_(particles_->acc_prior_),
			indicator_(particles_->indicator_), pos_(particles_->pos_),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			B_(*particles_->getVariableByName<Matd>("CorrectionMatrix"))//,
			//is_near_wall_P1_(*particles_->getVariableByName<int>("IsNearWallP1")), //for test
			//is_near_wall_P2_(*particles_->getVariableByName<int>("IsNearWallP2"))  //for test

		{
			particles_->registerVariable(tke_acc_inner_, "TkeAccInner");
			particles_->addVariableToWrite<Vecd>("TkeAccInner");
			particles_->registerVariable(tke_acc_wall_, "TkeAccWall");
			particles_->addVariableToWrite<Vecd>("TkeAccWall");

			particles_->registerVariable(test_k_grad_rslt_, "TkeGradResult");
			particles_->addVariableToWrite<Vecd>("TkeGradResult");
		}
		//=================================================================================================//
		TurbuViscousAccInner::TurbuViscousAccInner(BaseInnerRelation& inner_relation)
			: BaseViscousAccelerationInner(inner_relation),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			wall_Y_plus_(*particles_->getVariableByName<Real>("WallYplus")),
			velo_friction_(*particles_->getVariableByName<Vecd>("FrictionVelocity")),
			distance_to_wall_(*particles_->getVariableByName<Real>("DistanceToWall"))
		{
			particles_->registerVariable(visc_acc_inner_, "ViscousAccInner");
			particles_->addVariableToWrite<Vecd>("ViscousAccInner");
			particles_->registerVariable(visc_acc_wall_, "ViscousAccWall");
			particles_->addVariableToWrite<Vecd>("ViscousAccWall");

			particles_->registerVariable(shear_stress_, "ShearStress");
			//particles_->addVariableToWrite<Matd>("ShearStress");
			particles_->registerVariable(shear_stress_wall_, "ShearStressWall");
			//particles_->addVariableToWrite<Matd>("ShearStressWall");
		}
		//=================================================================================================//
		TurbulentEddyViscosity::
			TurbulentEddyViscosity(SPHBody& sph_body)
			: LocalDynamics(sph_body), FluidDataSimple(sph_body),
			rho_(particles_->rho_), wall_Y_star_(*particles_->getVariableByName<Real>("WallYstar")),
			wall_Y_plus_(*particles_->getVariableByName<Real>("WallYplus")),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")){}
		//=================================================================================================//
		void TurbulentEddyViscosity::update(size_t index_i, Real dt)
		{
			turbu_mu_[index_i] = rho_[index_i] * C_mu * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
		}
		
		//=================================================================================================//
		TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max* U_max),FluidDataSimple(sph_body), 
			vel_(particles_->vel_), advectionCFL_(advectionCFL),
			smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity"))
		{
			Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
			reference_ = SMAX(viscous_speed * viscous_speed, reference_);
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
			return advectionCFL_ * smoothing_length_min_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		InflowTurbulentCondition::InflowTurbulentCondition(BodyPartByCell& body_part
			, Real CharacteristicLength, Real relaxation_rate): 
			BaseFlowBoundaryCondition(body_part),
			relaxation_rate_(relaxation_rate), 
			CharacteristicLength_(CharacteristicLength), 
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation"))
		{
			TurbulentLength_ = 0.07 * CharacteristicLength_ / pow(C_mu, 0.75);
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
		Real InflowTurbulentCondition:: getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k)
		{
			Real u = velocity[0];
			Real temp_in_turbu_k = 1.5 * pow((TurbulentIntensity * u), 2);
			Real turbu_k_original = turbu_k;
			if (position[0] < 0.0)
			{
				turbu_k_original = temp_in_turbu_k;
				//std::cout << "temp_in_turbu_k="<< temp_in_turbu_k << std::endl;
			}
			return turbu_k_original;
		}
		//=================================================================================================//
		Real InflowTurbulentCondition::getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E)
		{
			//Real temp_in_turbu_E = C_mu * pow(turbu_k, 1.5) / (0.1*getTurbulentLength());
			Real temp_in_turbu_E = pow(turbu_k, 1.5) / TurbulentLength_;
			Real turbu_E_original = turbu_E;
			if (position[0] < 0.0)
			{
				turbu_E_original = temp_in_turbu_E;
			}
			return turbu_E_original;
		}
	}
	//=================================================================================================//
}
//=================================================================================================//