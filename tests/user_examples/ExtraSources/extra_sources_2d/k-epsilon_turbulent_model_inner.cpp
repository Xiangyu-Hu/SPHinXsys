#include "k-epsilon_turbulent_model_inner.hpp"
#include "k-epsilon_turbulent_model_inner.h"
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbuClosureCoeff::BaseTurbuClosureCoeff()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(1.0e-2), sigma_k(1.0),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793){}
		//=================================================================================================//
		BaseTurtbulentModelInner::BaseTurtbulentModelInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			particle_spacing_min_(inner_relation.real_body_->sph_adaptation_->MinimumSpacing()),
			 rho_(particles_->rho_), vel_(particles_->vel_), 
			mu_(particles_->fluid_.ReferenceViscosity()), dimension_(Vecd(0).size()), 
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()){}
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
				
				//** for test */
				particles_->registerVariable(lap_k_, "Lap_K");
				particles_->registerSortableVariable<Real>("Lap_K");
				particles_->addVariableToWrite<Real>("Lap_K");

				particles_->registerVariable(lap_k_term_, "Lap_K_term");
				particles_->registerSortableVariable<Real>("Lap_K_term");
				particles_->addVariableToWrite<Real>("Lap_K_term");

				particles_->addVariableToWrite<Real>("ChangeRateOfTKE");

				particles_->registerVariable(velocity_gradient, "Velocity_Gradient");
				particles_->registerSortableVariable<Matd>("Velocity_Gradient");
				particles_->addVariableToWrite<Matd>("Velocity_Gradient");
				
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
			surface_indicator_(particles_->surface_indicator_), pos_(particles_->pos_),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")){}
		//=================================================================================================//
		TurbuViscousAccInner::TurbuViscousAccInner(BaseInnerRelation& inner_relation)
			: BaseViscousAccelerationInner(inner_relation),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")) {}
		//=================================================================================================//
		TurbulentEddyViscosity::
			TurbulentEddyViscosity(SPHBody& sph_body)
			: LocalDynamics(sph_body), FluidDataSimple(sph_body),
			rho_(particles_->rho_), wall_Y_star_(*particles_->getVariableByName<Real>("WallYstar")),
			wall_Y_plus_(*particles_->getVariableByName<Real>("WallYplus")),
			is_near_wall_P1_(*particles_->getVariableByName<int>("IsNearWallP1")),
			mu_(particles_->fluid_.ReferenceViscosity()),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")){}
		//=================================================================================================//
		void TurbulentEddyViscosity::update(size_t index_i, Real dt)
		{
			if (is_near_wall_P1_[index_i] == 0)
			{
				turbu_mu_[index_i] = rho_[index_i] * C_mu * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
			}
			else //for the near wall particles, wall function effects
			{
				turbu_mu_[index_i] = wall_Y_star_[index_i] * mu_ * Karman / log(turbu_const_E * wall_Y_star_[index_i]);
			}
		}
		
		//=================================================================================================//
		TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max* U_max),FluidDataSimple(sph_body), 
			vel_(particles_->vel_),
			smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			advectionCFL_(advectionCFL), fluid_(particles_->fluid_),
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