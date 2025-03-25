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
        Real *physical_time_;
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
    };

  protected:
    SingularVariable<Real> *sv_physical_time_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    ConditionFunction condition_function_;
};

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
        void update(size_t index_i, Real dt = 0.0);

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
        Real rho0_;
        Real sound_speed_;
        Vecd *pos_;
        Real *rho_, *p_;
        int *buffer_particle_indicator_;
        ConditionKernel condition_;
        int *previous_surface_indicator_;
        Real *physical_time_;
        Real upper_bound_fringe_;
    };

  protected:
    int part_id_;
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    SpawnRealParticle create_real_particle_method_;
    Real rho0_;
    Real sound_speed_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
    TargetPressure target_pressure_method_;
    DiscreteVariable<int> *dv_previous_surface_indicator_;
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
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_H