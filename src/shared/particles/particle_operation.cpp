#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
SpawnRealParticle::SpawnRealParticle(BaseParticles *particles)
    : evolving_variables_(particles->EvolvingVariables()),
      dv_original_id_(particles->getVariableByName<UnsignedInt>("OriginalID")),
      group_manager_(particles->getParticleGroupManager()),
      sv_total_real_particles_(particles->svTotalRealParticles()),
      particles_bound_(particles->ParticlesBound()) {}
//=================================================================================================//
VariableArrayAssemble &SpawnRealParticle::getCopyableStates()
{

    if (!copyable_states_keeper_.getPtr())
    {
        VariableArrayAssemble &copyable_states =
            *copyable_states_keeper_.createPtr<VariableArrayAssemble>();
        OperationBetweenDataAssembles<
            DiscreteVariables, VariableArrayAssemble, VariableArrayAssembleInitialization>
            initialize_discrete_variable_array;
        initialize_discrete_variable_array(evolving_variables_, copyable_states);
    }

    return *copyable_states_keeper_.getPtr();
}
//=================================================================================================//
RemoveRealParticle::RemoveRealParticle(BaseParticles *particles)
    : evolving_variables_(particles->EvolvingVariables()),
      copyable_states_(),
      dv_original_id_(particles->getVariableByName<UnsignedInt>("OriginalID")),
      group_manager_(particles->getParticleGroupManager()),
      life_status_(group_manager_.getGroupMask("LifeStatus")),
      sv_total_real_particles_(particles->svTotalRealParticles()) {}
//=================================================================================================//
} // namespace SPH
