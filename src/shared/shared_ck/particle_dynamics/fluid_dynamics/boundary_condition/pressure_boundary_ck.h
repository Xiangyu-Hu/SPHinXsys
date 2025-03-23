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
 * @file 	pressure_boundary_ck.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef PRESSURE_BOUNDARY_CK_H
#define PRESSURE_BOUNDARY_CK_H

#include "fluid_boundary_ck.h"
#include "fluid_boundary_state.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... T>
class PressureConditionCK;
template <class AlignedBoxPartType, class KernelCorrectionType, class ConditionFunction>
class PressureConditionCK<AlignedBoxPartType, KernelCorrectionType, ConditionFunction>
    : public BaseLocalDynamics<AlignedBoxPartType>
{
    using ConditionKernel = typename ConditionFunction::ComputingKernel;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    PressureConditionCK(AlignedBoxPartType &aligned_box_part);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
        int *buffer_particle_indicator_;
        Vecd *zero_gradient_residue_;
        Vecd *pos_;
        Real *physical_time_;
        Vecd *vel_;
        CorrectionKernel correction_;
        Real *rho_;
        Transform *transform_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    ConditionFunction condition_function_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_;
    DiscreteVariable<Vecd> *dv_pos_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_vel_;
    KernelCorrectionType kernel_correction_;
    DiscreteVariable<Real> *dv_rho_;
};

template <typename ParallelExecutionPolicy, typename SequencedExecutionPolicy, class KernelCorrectionType, class PressureConditionFunction>
class PressureBidirectionalConditionCK
{
  public:
    StateDynamics<ParallelExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;

    StateDynamics<ParallelExecutionPolicy,
                  PressureConditionCK<AlignedBoxPartByCell, KernelCorrectionType, PressureConditionFunction>>
        pressure_condition_;

    StateDynamics<SequencedExecutionPolicy, BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, PressureConditionFunction>> emitter_injection_;

    StateDynamics<SequencedExecutionPolicy, DisposerOutflowDeletionCK> disposer_outflow_deletion_;

    PressureBidirectionalConditionCK(AlignedBoxPartByCell &emitter_by_cell,
                                     ParticleBuffer<Base> &inlet_buffer)
        : tag_buffer_particles_(emitter_by_cell),
          pressure_condition_(emitter_by_cell),
          emitter_injection_(emitter_by_cell, inlet_buffer),
          disposer_outflow_deletion_(emitter_by_cell)
    {
    }

    /// Tag (or flag) particles in the buffer.
    void tagBufferParticles() { tag_buffer_particles_.exec(); }

    /// Apply the pressure condition (note that this usually takes a time-step dt).
    void applyPressureCondition(Real dt) { pressure_condition_.exec(dt); }

    /// Perform the injection step (for inflow).
    void injectParticles() { emitter_injection_.exec(); }

    /// Perform the deletion step (for outflow).
    void deleteParticles() { disposer_outflow_deletion_.exec(); }
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // PRESSURE_BOUNDARY_CK_H