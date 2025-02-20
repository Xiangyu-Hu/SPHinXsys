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
 * @file 	fluid_boundary_ck.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef FLUID_BOUNDARY_CK_H
#define FLUID_BOUNDARY_CK_H

#include "base_body_part.h"
#include "base_data_package.h"
#include "base_fluid_dynamics.h"
#include "fluid_boundary_state.hpp"
#include "particle_operation.hpp"
#include "particle_reserve.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... T>
class InflowConditionCK;
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
class InflowConditionCK<AlignedBoxPartType, ConditionFunction>
    : public BaseLocalDynamics<AlignedBoxPartType>
{
    using ConditionKernel = typename ConditionFunction::ComputingKernel;

  public:
    InflowConditionCK(AlignedBoxPartType &aligned_box_part);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    ConditionFunction condition_function_;
};
//=================================================================================================//
template <typename AlignedBoxPartType>
class EmitterInflowInjectionCK : public BaseLocalDynamics<AlignedBoxPartType>
{
    using CreateRealParticleKernel = typename SpawnRealParticle::ComputingKernel;

  public:
    EmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer);
    virtual ~EmitterInflowInjectionCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        CreateRealParticleKernel create_real_particle_;
        Real rho0_;
        Vecd *pos_;
        Real *rho_, *p_;
    };

    class FinishDynamics
    {
        BaseParticles *particles_;
        ParticleBuffer<Base> &buffer_;

      public:
        FinishDynamics(EmitterInflowInjectionCK &encloser);
        void operator()();
    };

  protected:
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SpawnRealParticle create_real_particle_method_;
    Real rho0_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
};
//=================================================================================================//
template <typename AlignedBoxPartType>
class BufferEmitterInflowInjectionCK : public BaseLocalDynamics<AlignedBoxPartType>
{
    using CreateRealParticleKernel = typename SpawnRealParticle::ComputingKernel;

  public:
    BufferEmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer);
    virtual ~BufferEmitterInflowInjectionCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        CreateRealParticleKernel create_real_particle_;
        Real rho0_;
        Vecd *pos_;
        Real *rho_, *p_;
        int *buffer_particle_indicator_;
    };

    class FinishDynamics
    {
        BaseParticles *particles_;
        ParticleBuffer<Base> &buffer_;

      public:
        FinishDynamics(BufferEmitterInflowInjectionCK &encloser);
        void operator()();
    };

  protected:
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SpawnRealParticle create_real_particle_method_;
    Real rho0_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
};
//=================================================================================================//
class DisposerOutflowDeletionCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{
    using RemoveRealParticleKernel = typename DespawnRealParticle::ComputingKernel;

  public:
    DisposerOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part);
    virtual ~DisposerOutflowDeletionCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        RemoveRealParticleKernel remove_real_particle_;
        Real rho0_;
        Vecd *pos_;
        Real *rho_, *p_;
    };

    class FinishDynamics
    {
        BaseParticles *particles_;

      public:
        FinishDynamics(DisposerOutflowDeletionCK &encloser);
        void operator()();
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    DespawnRealParticle remove_real_particle_method_;
    Real rho0_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
};
} // namespace fluid_dynamics
} // namespace SPH
namespace SPH
{
namespace fluid_dynamics
{
struct NonPrescribedPressure
{
    template <class BoundaryConditionType>
    NonPrescribedPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

class TagBufferParticlesCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{

  public:
    TagBufferParticlesCK(AlignedBoxPartByCell &aligned_box_part);
    virtual ~TagBufferParticlesCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        int *part_id_;
        Vecd *pos_;
        int *buffer_particle_indicator_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
};

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

template <typename ExecutionPolicy, class KernelCorrectionType, class ConditionFunction>
class BidirectionalBufferCK
{
  public:
    StateDynamics<ExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;

    StateDynamics<ExecutionPolicy,
                  PressureConditionCK<AlignedBoxPartByCell, KernelCorrectionType, ConditionFunction>>
        pressure_condition_;

    StateDynamics<execution::SequencedPolicy, EmitterInflowInjectionCK<AlignedBoxPartByCell>> emitter_injection_;

    StateDynamics<execution::SequencedPolicy, DisposerOutflowDeletionCK> disposer_outflow_deletion_;

    BidirectionalBufferCK(AlignedBoxPartByCell &emitter_by_cell,
                          ParticleBuffer<Base> &inlet_buffer)
        : tag_buffer_particles_(emitter_by_cell),
          pressure_condition_(emitter_by_cell),
          emitter_injection_(emitter_by_cell, inlet_buffer),
          disposer_outflow_deletion_(emitter_by_cell)
    {
    }

    // Member functions to call the internal dynamics

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
#endif // FLUID_BOUNDARY_CK_H