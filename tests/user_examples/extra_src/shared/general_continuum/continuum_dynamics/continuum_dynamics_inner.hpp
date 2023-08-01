#pragma once
#include "continuum_dynamics_inner.h"
namespace SPH
{
	namespace continuum_dynamics
	{
        template <class FluidDynamicsType>
        BaseIntegration1stHalf<FluidDynamicsType>::BaseIntegration1stHalf(BaseInnerRelation& inner_relation)
            : FluidDynamicsType(inner_relation), 
            acc_shear_(*this->particles_->template getVariableByName<Vecd>("AccelerationByShear")) {}

        template <class FluidDynamicsType>
        void BaseIntegration1stHalf<FluidDynamicsType>::initialization(size_t index_i, Real dt)
        {
            FluidDynamicsType::initialization(index_i, dt);
        }
        template <class FluidDynamicsType>
        void BaseIntegration1stHalf<FluidDynamicsType>::interaction(size_t index_i, Real dt)
        {
            FluidDynamicsType::interaction(index_i, dt);
        }
        template <class FluidDynamicsType>
        void BaseIntegration1stHalf<FluidDynamicsType>::update(size_t index_i, Real dt)
        {
            this->vel_[index_i] += (this->acc_prior_[index_i] + this->acc_[index_i] + this->acc_shear_[index_i]) * dt;
        }
	}
} // namespace SPH