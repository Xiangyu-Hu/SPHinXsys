/**
 * @file 	eulerian_compressible_fluid_dynamics.cpp
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
 */

#include "eulerian_compressible_fluid_dynamics_inner.h"

//=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_compressible_fluid_dynamics
	{
		//=================================================================================================//
		CompressibleFlowTimeStepInitialization::CompressibleFlowTimeStepInitialization(SPHBody &sph_body)
			: ParticleDynamicsSimple(sph_body), CompressibleFluidDataSimple(sph_body),
			  rho_n_(particles_->rho_n_), dE_dt_prior_(particles_->dE_dt_prior_), pos_n_(particles_->pos_n_),
			  vel_n_(particles_->vel_n_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  gravity_(gravity_ptr_keeper_.createPtr<Gravity>(Vecd(0))) {}
		//=================================================================================================//
		CompressibleFlowTimeStepInitialization::
			CompressibleFlowTimeStepInitialization(SPHBody &sph_body, Gravity &gravity)
			: ParticleDynamicsSimple(sph_body), CompressibleFluidDataSimple(sph_body),
			  rho_n_(particles_->rho_n_), dE_dt_prior_(particles_->dE_dt_prior_), pos_n_(particles_->pos_n_),
			  vel_n_(particles_->vel_n_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  gravity_(&gravity) {}
		//=================================================================================================//
		void CompressibleFlowTimeStepInitialization::setupDynamics(Real dt)
		{
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void CompressibleFlowTimeStepInitialization::Update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_n_[index_i] * gravity_->InducedAcceleration(pos_n_[index_i]);
			dE_dt_prior_[index_i] = rho_n_[index_i] * SimTK::dot(gravity_->InducedAcceleration(pos_n_[index_i]), vel_n_[index_i]);
		}
		//=================================================================================================//
		CompressibleFluidInitialCondition::
			CompressibleFluidInitialCondition(EulerianFluidBody &body)
			: ParticleDynamicsSimple(body), CompressibleFluidDataSimple(body),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), mom_(particles_->mom_),
			  rho_n_(particles_->rho_n_), E_(particles_->E_), p_(particles_->p_),
			  gamma_(material_->HeatCapacityRatio()) {}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  CompressibleFluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
			  mass_(particles_->mass_), dE_dt_prior_(particles_->dE_dt_prior_),
			  vel_n_(particles_->vel_n_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()),
			  mu_(material_->ReferenceViscosity()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			const Vecd &vel_i = vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//viscous force
				vel_derivative = (vel_i - vel_n_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}
			dmom_dt_prior_[index_i] += rho_n_[index_i] * acceleration;
			dE_dt_prior_[index_i] += rho_n_[index_i] * dot(acceleration, vel_n_[index_i]);
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(EulerianFluidBody &body)
			: ParticleDynamicsReduce<Real, ReduceMax>(body),
			  CompressibleFluidDataSimple(body), rho_n_(particles_->rho_n_),
			  p_(particles_->p_), vel_n_(particles_->vel_n_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength())
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return material_->getSoundSpeed(p_[index_i], rho_n_[index_i]) + vel_n_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			particles_->signal_speed_max_ = reduced_value;
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		BaseRelaxation::BaseRelaxation(BaseBodyRelationInner &inner_relation)
			: ParticleDynamics1Level(*inner_relation.sph_body_),
			  CompressibleFluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
			  drho_dt_(particles_->drho_dt_), E_(particles_->E_), dE_dt_(particles_->dE_dt_),
			  dE_dt_prior_(particles_->dE_dt_prior_),
			  vel_n_(particles_->vel_n_), mom_(particles_->mom_),
			  dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BasePressureRelaxation::Initialization(size_t index_i, Real dt)
		{
			E_[index_i] += dE_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Real rho_e = E_[index_i] - 0.5 * mom_[index_i].normSqr() / rho_n_[index_i];
			p_[index_i] = material_->getPressure(rho_n_[index_i], rho_e);
		}
		//=================================================================================================//
		void BasePressureRelaxation::Update(size_t index_i, Real dt)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_n_[index_i] = mom_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		BaseDensityAndEnergyRelaxation::
			BaseDensityAndEnergyRelaxation(BaseBodyRelationInner &inner_relation)
			: BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BaseDensityAndEnergyRelaxation::Update(size_t index_i, Real dt)
		{
			E_[index_i] += dE_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//