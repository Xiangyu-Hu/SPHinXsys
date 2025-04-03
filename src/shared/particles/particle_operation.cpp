#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
SpawnRealParticle::SpawnRealParticle(BaseParticles *particles)
    : evolving_variables_(particles->EvolvingVariables()),
      copyable_states_(),
      dv_original_id_(particles->getVariableByName<UnsignedInt>("OriginalID")),
      sv_total_real_particles_(particles->svTotalRealParticles()),
      particles_bound_(particles->ParticlesBound()) {}
//=================================================================================================//
RemoveRealParticle::RemoveRealParticle(BaseParticles *particles)
    : evolving_variables_(particles->EvolvingVariables()),
      copyable_states_(),
      dv_original_id_(particles->getVariableByName<UnsignedInt>("OriginalID")),
      sv_total_real_particles_(particles->svTotalRealParticles()){}
} // namespace SPH
