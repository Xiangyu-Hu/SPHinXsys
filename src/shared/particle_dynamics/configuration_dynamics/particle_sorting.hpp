#ifndef PARTICLE_SORTING_HPP
#define PARTICLE_SORTING_HPP

#include "particle_sorting.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
ParticleSorting<ExecutionPolicy>::ParticleSorting(RealBody &real_body)
    : BaseDynamics<void>(),
      particle_sequence_(real_body), particle_data_sort_(real_body),
      update_sorted_id_(real_body) {}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleSorting<ExecutionPolicy>::exec(Real dt)
{
    particle_sequence_.exec();
    particle_data_sort_.exec();
    update_sorted_id_.exec();
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORTING_HPP
