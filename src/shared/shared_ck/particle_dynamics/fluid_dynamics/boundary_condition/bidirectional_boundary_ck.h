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
 * @file 	bidirectional_boundary_ck.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef BIDIRECTIONAL_BOUNDARY_CK_H
#define BIDIRECTIONAL_BOUNDARY_CK_H

#include "base_body_part.h"
#include "base_fluid_dynamics.h"
#include "fluid_boundary_state.hpp"
#include "particle_operation.hpp"
#include "particle_reserve.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
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
        void update(size_t index_i, Real dt = 0.0);

      protected:
        int part_id_;
        AlignedBox *aligned_box_;
        Vecd *pos_;
        int *buffer_particle_indicator_;
    };

  protected:
    int part_id_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
};

template <class TargetPressure>
class BufferInflowInjectionCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{
    using CreateRealParticleKernel = typename SpawnRealParticle::ComputingKernel;
    using PressureKernel = typename TargetPressure::ComputingKernel;

  public:
    BufferInflowInjectionCK(AlignedBoxPartByCell &aligned_box_part, ParticleBuffer<Base> &buffer);
    virtual ~BufferInflowInjectionCK() {};

    class FinishDynamics
    {
      public:
        FinishDynamics(BufferInflowInjectionCK &encloser)
            : particles_(encloser.particles_), buffer_(encloser.buffer_) {}
        void operator()() { buffer_.checkEnoughBuffer(*particles_); }

      private:
        BaseParticles *particles_;
        ParticleBuffer<Base> &buffer_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      private:
        int part_id_;
        AlignedBox *aligned_box_;
        UnsignedInt *total_real_particles_;
        CreateRealParticleKernel create_real_particle_;
        Vecd *pos_;
        int *buffer_particle_indicator_;
        PressureKernel pressure_kernel_;
        Real *physical_time_;
        Real upper_bound_fringe_;
    };

  protected:
    int part_id_;
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    SpawnRealParticle create_real_particle_method_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
    TargetPressure target_pressure_method_;
    SingularVariable<Real> *sv_physical_time_;
    Real upper_bound_fringe_;
};

class BufferOutflowDeletionCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{
    using RemoveRealParticleKernel = typename DeleteRealParticle::ComputingKernel;

  public:
    BufferOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part);
    virtual ~BufferOutflowDeletionCK() {};

    class UpdateKernel
    {
        struct IsDeletable
        {
            int part_id_;
            AlignedBox *aligned_box_;
            Vecd *pos_;
            int *buffer_particle_indicator_;
            IsDeletable(int part_id, AlignedBox *aligned_box,
                        Vecd *pos, int *buffer_particle_indicator)
                : part_id_(part_id), aligned_box_(aligned_box), pos_(pos),
                  buffer_particle_indicator_(buffer_particle_indicator) {}
            bool operator()(size_t index_i) const
            {
                return buffer_particle_indicator_[index_i] == part_id_ &&
                       aligned_box_->checkLowerBound(pos_[index_i]);
            }
        };

      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        Vecd *pos_;
        IsDeletable is_deltable_;
        UnsignedInt *total_real_particles_;
        RemoveRealParticleKernel remove_real_particle_;
    };

  protected:
    int part_id_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    DeleteRealParticle remove_real_particle_method_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
    DiscreteVariable<Vecd> *dv_pos_;
};

template <class TargetVelocity>
class VelocityConditionCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{
    using VelocityKernel = typename TargetVelocity::ComputingKernel;

  public:
    VelocityConditionCK(AlignedBoxPartByCell &aligned_box_part);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        Transform *transform_;
        VelocityKernel velocity_kernel_;
        Vecd *dv_pos_;
        Vecd *dv_vel_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    TargetVelocity target_velocity_method_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Vecd> *dv_vel_;
};

template <class KernelCorrectionType, class TargetPressure>
class PressureConditionCK : public BaseLocalDynamics<AlignedBoxPartByCell>,
                            public ForcePriorCK
{
    using PressureKernel = typename TargetPressure::ComputingKernel;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    PressureConditionCK(AlignedBoxPartByCell &aligned_box_part);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        CorrectionKernel correction_;
        PressureKernel pressure_kernel_;
        Real *physical_time_;
        Vecd *zero_gradient_residue_;
        Vecd *pos_;
        Real *mass_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    KernelCorrectionType kernel_correction_;
    TargetPressure target_pressure_method_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_mass_;
};

template <typename ExecutionPolicy, class KernelCorrectionType, class TargetPressure, class TargetVelocity>
class BidirectionalBoundaryCK
{
    StateDynamics<ExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;
    StateDynamics<ExecutionPolicy, VelocityConditionCK<TargetVelocity>> velocity_condition_;
    StateDynamics<ExecutionPolicy, PressureConditionCK<KernelCorrectionType, TargetPressure>> pressure_condition_;
    StateDynamics<ExecutionPolicy, BufferInflowInjectionCK<TargetPressure>> inflow_injection_;
    StateDynamics<ExecutionPolicy, BufferOutflowDeletionCK> outflow_deletion_;

  public:
    BidirectionalBoundaryCK(AlignedBoxPartByCell &emitter_by_cell, ParticleBuffer<Base> &particle_buffer);
    void tagBufferParticles() { tag_buffer_particles_.exec(); }
    void applyVelocityCondition() { velocity_condition_.exec(); }
    void applyPressureCondition() { pressure_condition_.exec(); }
    void injectParticles() { inflow_injection_.exec(); }
    void deleteParticles() { outflow_deletion_.exec(); }
};

template <typename ExecutionPolicy, class KernelCorrectionType, class TargetPressure>
using BidirectionalPressureBoundaryCK =
    BidirectionalBoundaryCK<ExecutionPolicy, KernelCorrectionType, TargetPressure, AlignedVelocity>;

template <typename ExecutionPolicy, class KernelCorrectionType, class TargetVelocity>
using BidirectionalVelocityBoundaryCK =
    BidirectionalBoundaryCK<ExecutionPolicy, KernelCorrectionType, ZeroGradientPressure, TargetVelocity>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BOUNDARY_CK_H