#include "particle_sort_ck.hpp"

namespace SPH
{
//=================================================================================================//
UpdateSortableVariables::UpdateSortableVariables(BaseParticles *particles)
    : initialize_temp_variables_()
{
    initialize_temp_variables_(temp_variables_, particles->ParticlesBound());
}
//=================================================================================================//
} // namespace SPH