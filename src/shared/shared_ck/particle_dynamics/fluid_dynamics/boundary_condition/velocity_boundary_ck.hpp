#ifndef VELOCITY_BOUNDARY_CK_HPP
#define VELOCITY_BOUNDARY_CK_HPP

#include "velocity_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <typename ParallelExecutionPolicy, class KernelCorrectionType, class VelocityConditionFunction>
VelocityBidirectionalConditionCK<
    ParallelExecutionPolicy, KernelCorrectionType, VelocityConditionFunction>::
    VelocityBidirectionalConditionCK(AlignedBoxPartByCell &emitter_by_cell,
                                     ParticleBuffer<Base> &inlet_buffer)
    : tag_buffer_particles_(emitter_by_cell),
      velocity_condition_(emitter_by_cell),
      pressure_condition_(emitter_by_cell),
      emitter_injection_(emitter_by_cell, inlet_buffer),
      disposer_outflow_deletion_(emitter_by_cell) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // VELOCITY_BOUNDARY_CK_HPP
