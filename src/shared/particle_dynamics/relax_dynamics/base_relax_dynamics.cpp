#include "base_relax_dynamics.h"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
RandomizeParticlePosition::RandomizeParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      randomize_scale_(sph_body.getSPHAdaptation().MinimumSpacing()) {}
//=================================================================================================//
void RandomizeParticlePosition::update(size_t index_i, Real dt)
{
    Vecd &pos_n_i = pos_[index_i];
    for (int k = 0; k < pos_n_i.size(); ++k)
    {
        pos_n_i[k] += dt * rand_uniform(-1.0, 1.0) * randomize_scale_;
    }
}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
