
#include "common_compressible_FVM_classes.h"

namespace SPH
{
	//=================================================================================================//
	CompressibleAcousticTimeStepSizeInFVM::CompressibleAcousticTimeStepSizeInFVM(SPHBody& sph_body)
		: AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_), 
		compressible_fluid_(CompressibleFluid(1.0, 1.4)) {};
	//=================================================================================================//
	Real CompressibleAcousticTimeStepSizeInFVM::reduce(size_t index_i, Real dt)
	{
		return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
	}
	//=================================================================================================//
	Real CompressibleAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
	{
		return 0.0005 / (reduced_value + TinyReal);
	}
	//=================================================================================================//
	BaseIntegrationInCompressibleFVM::BaseIntegrationInCompressibleFVM(BaseInnerRelationInFVM& inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<BaseParticles>(inner_relation), 
		compressible_fluid_(CompressibleFluid(1.0, 1.4)), E_(*particles_->getVariableByName<Real>("TotalEnergy")),
		dE_dt_(*particles_->getVariableByName<Real>("TotalEnergyChangeRate")),
		dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")), 
		rho_(particles_->rho_), drho_dt_(*particles_->registerSharedVariable<Real>("DensityChangeRate")), p_(*particles_->getVariableByName<Real>("Pressure")),
		mom_(*particles_->getVariableByName<Vecd>("Momentum")),
		dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
		dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")), vel_(particles_->vel_), pos_(particles_->pos_) {};
	//=================================================================================================//
}
//=================================================================================================//