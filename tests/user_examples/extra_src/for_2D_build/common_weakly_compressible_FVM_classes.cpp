
#include "common_weakly_compressible_FVM_classes.h"

namespace SPH
{
	//=================================================================================================//
	WCAcousticTimeStepSizeInFVM::WCAcousticTimeStepSizeInFVM(SPHBody &sph_body) 
		: AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(particles_->p_), vel_(particles_->vel_), fluid_(particles_->fluid_){};
	//=================================================================================================//
	Real WCAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
	{
        // I chose a time-step size according to Eulerian method
        return 0.1 / (reduced_value + TinyReal);
	}
	//=================================================================================================//
	BaseRelaxationInFVM::BaseRelaxationInFVM(BaseInnerRelationInFVM &inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<FluidParticles>(inner_relation),
		fluid_(particles_->fluid_), p_(particles_->p_), rho_(particles_->rho_), drho_dt_(particles_->drho_dt_), 
		vel_(particles_->vel_), mom_(*particles_->getVariableByName<Vecd>("Momentum")),
		dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
		dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
		pos_(particles_->pos_), mu_(particles_->fluid_.ReferenceViscosity()){};
	//=================================================================================================//
	BaseForceFromFluidInFVM::BaseForceFromFluidInFVM(BaseInnerRelationInFVM& inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<FluidParticles>(inner_relation), Vol_(particles_->Vol_) {};
	//=================================================================================================//
	ViscousForceFromFluidInFVM::ViscousForceFromFluidInFVM(BaseInnerRelationInFVM& inner_relation)
		: BaseForceFromFluidInFVM(inner_relation), fluid_(particles_->fluid_), 
		vel_(particles_->vel_), mu_(particles_->fluid_.ReferenceViscosity())
	{
		particles_->registerVariable(force_from_fluid_, "ViscousForceFromFluid");
	};
	//=================================================================================================//
	void ViscousForceFromFluidInFVM::interaction(size_t index_i, Real dt)
	{
		Real Vol_i = Vol_[index_i];
		const Vecd &vel_i = vel_[index_i];
		Vecd force = Vecd::Zero();
		const NeighborhoodInFVM &inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				Vecd vel_j = -vel_i;
				Vecd vel_derivative = (vel_j - vel_i) / (inner_neighborhood.r_ij_[n] + TinyReal);
				force += 2.0 * mu_ * vel_derivative * Vol_i * inner_neighborhood.dW_ijV_j_[n];
			}
		}
		force_from_fluid_[index_i] = force;
	}
	//=================================================================================================//
}
//=================================================================================================//