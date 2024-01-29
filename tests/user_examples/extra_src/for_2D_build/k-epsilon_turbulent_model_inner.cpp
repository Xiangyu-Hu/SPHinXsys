#include "k-epsilon_turbulent_model_inner.hpp"
#include "k-epsilon_turbulent_model_inner.h"
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbuClosureCoeffInner::BaseTurbuClosureCoeffInner()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(5.0e-2), sigma_k(1.0),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793) {}
		//=================================================================================================//
		BaseTurtbulentModelInner::BaseTurtbulentModelInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			particle_spacing_min_(inner_relation.real_body_->sph_adaptation_->MinimumSpacing()),
			rho_(particles_->rho_), vel_(particles_->vel_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), dimension_(Vecd(0).size()),
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//

		//=================================================================================================//

		//=================================================================================================//

		//=================================================================================================//

	}
	//=================================================================================================//
}
//=================================================================================================//