#include "eulerian_compressible_fluid_dynamics_inner.h"

using namespace std;
//=========================================================================================================//
namespace SPH
{
	//=====================================================================================================//
	namespace eulerian_compressible_fluid_dynamics
	{
		//=================================================================================================//
		CompressibleFlowTimeStepInitialization::
			CompressibleFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
			: BaseTimeStepInitialization(sph_body, gravity_ptr), CompressibleFluidDataSimple(sph_body),
			  rho_(particles_->rho_), dE_dt_prior_(particles_->dE_dt_prior_), pos_(particles_->pos_),
			  vel_(particles_->vel_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		void CompressibleFlowTimeStepInitialization::update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
			dE_dt_prior_[index_i] = rho_[index_i] * (gravity_->InducedAcceleration(pos_[index_i])).dot(vel_[index_i]);
		}
		//=================================================================================================//
		CompressibleFluidInitialCondition::
			CompressibleFluidInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body), CompressibleFluidDataSimple(sph_body),
			  pos_(particles_->pos_), vel_(particles_->vel_), mom_(particles_->mom_),
			  rho_(particles_->rho_), E_(particles_->E_), p_(particles_->p_),
			  gamma_(particles_->compressible_fluid_.HeatCapacityRatio()) {}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), CompressibleFluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			  mass_(particles_->mass_), dE_dt_prior_(particles_->dE_dt_prior_),
			  vel_(particles_->vel_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
			  mu_(particles_->compressible_fluid_.ReferenceViscosity()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd &vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				// viscous force
				vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative  * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}
			dmom_dt_prior_[index_i] += rho_[index_i] * acceleration;
			dE_dt_prior_[index_i] += rho_[index_i] * acceleration.dot(vel_[index_i]);
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  CompressibleFluidDataSimple(sph_body),
			  compressible_fluid_(particles_->compressible_fluid_),
			  rho_(particles_->rho_), p_(particles_->p_), vel_(particles_->vel_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		BaseIntegration::BaseIntegration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), CompressibleFluidDataInner(inner_relation),
			  compressible_fluid_(particles_->compressible_fluid_),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			  drho_dt_(particles_->drho_dt_), E_(particles_->E_), dE_dt_(particles_->dE_dt_),
			  dE_dt_prior_(particles_->dE_dt_prior_),
			  vel_(particles_->vel_), mom_(particles_->mom_),
			  dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//