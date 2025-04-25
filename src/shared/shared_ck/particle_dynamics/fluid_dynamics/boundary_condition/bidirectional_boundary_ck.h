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
class BufferIndicationCK : public BaseLocalDynamics<AlignedBoxByCell>
{

  public:
    BufferIndicationCK(AlignedBoxByCell &aligned_box_part);
    virtual ~BufferIndicationCK() {};

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
        int *buffer_indicator_;
    };

  protected:
    int part_id_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_buffer_indicator_;
};

template <class ConditionType>
class BufferInflowInjectionCK : public BaseLocalDynamics<AlignedBoxByCell>
{
    using SpawnRealParticleKernel = typename SpawnRealParticle::ComputingKernel;
    using FluidType = typename ConditionType::Fluid;
    using EosKernel = typename FluidType::EosKernel;

  public:
    template <typename... Args>
    BufferInflowInjectionCK(AlignedBoxByCell &aligned_box_part,
                            ParticleBuffer<Base> &buffer, Args &&...args);
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
        EosKernel eos_;
        ConditionType condition_;
        AlignedBox *aligned_box_;
        UnsignedInt *total_real_particles_;
        SpawnRealParticleKernel spawn_real_particle_;
        Vecd *pos_;
        int *buffer_indicator_;
        Real *physical_time_;
        Real *p_, *rho_;
        Real upper_bound_fringe_;
    };

  protected:
    int part_id_;
    ParticleBuffer<Base> &buffer_;
    FluidType &fluid_;
    ConditionType condition_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    SpawnRealParticle spawn_real_particle_method_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_buffer_indicator_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Real> *dv_p_, *dv_rho_;
    Real upper_bound_fringe_;
};

class BufferOutflowDeletionCK : public BaseLocalDynamics<AlignedBoxByCell>
{
    using RemoveRealParticleKernel = typename RemoveRealParticle::ComputingKernel;

  public:
    BufferOutflowDeletionCK(AlignedBoxByCell &aligned_box_part);
    virtual ~BufferOutflowDeletionCK() {};

    class UpdateKernel
    {
        struct IsDeletable
        {
            int part_id_;
            AlignedBox *aligned_box_;
            Vecd *pos_;
            int *buffer_indicator_;
            IsDeletable(int part_id, AlignedBox *aligned_box, Vecd *pos, int *buffer_particle_indicator);

            bool operator()(size_t index_i) const
            {
                return buffer_indicator_[index_i] == part_id_ &&
                       aligned_box_->checkLowerBound(pos_[index_i]);
            };
        };

      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

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
    RemoveRealParticle remove_real_particle_method_;
    DiscreteVariable<int> *dv_buffer_indicator_;
    DiscreteVariable<Vecd> *dv_pos_;
};

template <class KernelCorrectionType, typename ConditionType>
class PressureVelocityCondition : public BaseLocalDynamics<AlignedBoxByCell>,
                                  public BaseStateCondition
{
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <typename... Args>
    PressureVelocityCondition(AlignedBoxByCell &aligned_box_part, Args &&...args);

    class UpdateKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        CorrectionKernel correction_kernel_;
        ConditionType condition_;
        Real *physical_time_;
        Vecd *zero_gradient_residue_;
        int axis_;
        Transform *transform_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    KernelCorrectionType kernel_correction_method_;
    ConditionType condition_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_;
};

template <typename ExecutionPolicy, class KernelCorrectionType, class ConditionType>
class BidirectionalBoundaryCK
{
    StateDynamics<ExecutionPolicy, BufferIndicationCK> tag_buffer_particles_;
    StateDynamics<ExecutionPolicy, PressureVelocityCondition<KernelCorrectionType, ConditionType>> boundary_condition_;
    StateDynamics<ExecutionPolicy, BufferInflowInjectionCK<ConditionType>> inflow_injection_;
    StateDynamics<ExecutionPolicy, BufferOutflowDeletionCK> outflow_deletion_;

  public:
    template <typename... Args>
    BidirectionalBoundaryCK(AlignedBoxByCell &aligned_box_part,
                            ParticleBuffer<Base> &particle_buffer, Args &&...args);
    void tagBufferParticles() { tag_buffer_particles_.exec(); }
    void applyBoundaryCondition(Real dt) { boundary_condition_.exec(dt); }
    void injectParticles() { inflow_injection_.exec(); }
    void deleteParticles() { outflow_deletion_.exec(); }
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BOUNDARY_CK_H
