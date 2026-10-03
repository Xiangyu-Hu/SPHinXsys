#pragma once

#include "udf_k-omega_turbulent_model.h"

namespace SPH
{

namespace fluid_dynamics
{

namespace udf
{

template <class DataDelegationType>
template <class BaseRelationType>
kOmega_BaseTurbulentModel<Base, DataDelegationType>::kOmega_BaseTurbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      viscosity_(this->sph_body_-> template getMaterialProperty<Viscosity>()),
      mu_(viscosity_.ReferenceViscosity()),
      smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().getSPHAdaptation().MinimumSpacing()),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}

template <class DataDelegationType>
template <class BaseRelationType>
kOmega_GetVelocityGradient<DataDelegationType>::
    kOmega_GetVelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(this->particles_->template registerStateVariableData<Matd>("TurbulentVelocityGradient")),
      velocity_gradient_wall(this->particles_->template registerStateVariableData<Matd>("Velocity_Gradient_Wall")) {}

}

}

}
