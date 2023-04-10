#include "k-epsilon_turbulent_model_inner.hpp"
#include "k-epsilon_turbulent_model_inner.h"
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbulentClosureCoefficient::BaseTurbulentClosureCoefficient()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(1.0e-2), sigma_k(1.0),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793), TurbulentLength(0.0){}
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
		TurbulentKineticEnergyAccelerationInner::
			TurbulentKineticEnergyAccelerationInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation), acc_prior_(particles_->acc_prior_),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")){}
		//=================================================================================================//




	}
	//=================================================================================================//
}
//=================================================================================================//