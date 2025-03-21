/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	velocity_boundary_ck.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef VELOCITY_BOUNDARY_CK_H
#define VELOCITY_BOUNDARY_CK_H

#include "fluid_boundary_ck.h"
#include "fluid_boundary_state.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename ParallelExecutionPolicy, typename SequencedExecutionPolicy, class KernelCorrectionType, class VelocityConditionFunction>
class VelocityBidirectionalConditionCK
{
  public:
    StateDynamics<ParallelExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;

    StateDynamics<ParallelExecutionPolicy, fluid_dynamics::InflowConditionCK<AlignedBoxPartByCell, VelocityConditionFunction>> velocity_condition_;

    StateDynamics<ParallelExecutionPolicy, PressureConditionCK<AlignedBoxPartByCell, KernelCorrectionType, DummyPressure>>
        pressure_condition_;

    StateDynamics<ParallelExecutionPolicy, BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, NonPrescribedPressure>> emitter_injection_;

    StateDynamics<ParallelExecutionPolicy, DisposerOutflowDeletionCK> disposer_outflow_deletion_;

    VelocityBidirectionalConditionCK(AlignedBoxPartByCell &emitter_by_cell, ParticleBuffer<Base> &inlet_buffer);

    /// Tag (or flag) particles in the buffer.
    void tagBufferParticles() { tag_buffer_particles_.exec(); }

    /// Apply the pressure condition (note that this usually takes a time-step dt).
    void applyPressureCondition(Real dt) { pressure_condition_.exec(dt); }

    /// Apply the velocity condition.
    void applyVelocityCondition() { velocity_condition_.exec(); }

    /// Perform the injection step (for inflow).
    void injectParticles() { emitter_injection_.exec(); }

    /// Perform the deletion step (for outflow).
    void deleteParticles() { disposer_outflow_deletion_.exec(); }
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // VELOCITY_BOUNDARY_CK_H