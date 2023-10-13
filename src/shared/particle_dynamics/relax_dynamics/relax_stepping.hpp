#ifndef RELAX_STEPPING_HPP
#define RELAX_STEPPING_HPP

#include "relax_stepping.h"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
ParticleRelaxation<Base, DataDelegationType>::ParticleRelaxation(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      sph_adaptation_(sph_body.sph_adaptation_),
      pos_(this->particles_->rho_) {}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_STEPPING_HPP
