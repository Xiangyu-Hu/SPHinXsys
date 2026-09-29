

#pragma once

#include "udf_k-epsilon_turbulent_model.h"

namespace SPH
{

namespace fluid_dynamics
{
namespace udf
{

template <class DataDelegationType>
template <class BaseRelationType>
kEpsilon_BaseTurbulentModel<Base, DataDelegationType>::kEpsilon_BaseTurbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      turbu_strain_rate_(this->particles_->template registerStateVariableData<Matd>("TurbulentStrainRate")),
      viscosity_(this->sph_body_-> template getMaterialProperty<Viscosity>()),
      mu_(viscosity_.ReferenceViscosity()),
      smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().getSPHAdaptation().MinimumSpacing()),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}

}
}

}
