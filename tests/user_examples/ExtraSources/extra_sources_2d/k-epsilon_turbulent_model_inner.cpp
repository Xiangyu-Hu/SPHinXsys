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
			 rho_(particles_->rho_), 
			vel_(particles_->vel_), 
			mu_(particles_->fluid_.ReferenceViscosity()), dimension_(Vecd(0).size()), 
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()){}
		//=================================================================================================//
		K_TurtbulentModelInner::K_TurtbulentModelInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation) 
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
		}
		//=================================================================================================//
		void K_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_k_[index_i] += dk_dt_[index_i] * dt;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//