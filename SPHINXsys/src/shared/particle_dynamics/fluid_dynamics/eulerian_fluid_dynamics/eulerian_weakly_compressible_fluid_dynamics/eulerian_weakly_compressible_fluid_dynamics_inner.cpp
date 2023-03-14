#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.hpp"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{	//=================================================================================================//
		WeaklyCompressibleFluidInitialCondition::
			WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body), EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			  pos_(particles_->pos_), vel_(particles_->vel_), mom_(particles_->mom_),
			  rho_(particles_->rho_), p_(particles_->p_) {}
		//=================================================================================================//
		EulerianFlowTimeStepInitialization::
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
			: BaseTimeStepInitialization(sph_body, gravity_ptr),
			  EulerianWeaklyCompressibleFluidDataSimple(sph_body), rho_(particles_->rho_),
			  pos_(particles_->pos_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		void EulerianFlowTimeStepInitialization::update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()),
			  EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			  vel_(particles_->vel_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  mu_(particles_->fluid_.ReferenceViscosity()),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			  fluid_(particles_->fluid_), rho_(particles_->rho_),
			  p_(particles_->p_), vel_(particles_->vel_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		BaseIntegration::BaseIntegration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()),
			  EulerianWeaklyCompressibleFluidDataInner(inner_relation), fluid_(particles_->fluid_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_), rho_(particles_->rho_),
			  p_(particles_->p_), drho_dt_(particles_->drho_dt_), vel_(particles_->vel_), mom_(particles_->mom_),
			  dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//