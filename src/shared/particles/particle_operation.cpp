#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
SpawnRealParticle::SpawnRealParticle(BaseParticles *particles)
    : variables_to_sort_(particles->VariablesToSort()),
      copyable_states_(),
      dv_original_id_(particles->getVariableByName<UnsignedInt>("OriginalID")),
      sv_total_real_particles_(particles->svTotalRealParticles()),
      real_particles_bound_(particles->RealParticlesBound()) {}
//=================================================================================================//
} // namespace SPH
