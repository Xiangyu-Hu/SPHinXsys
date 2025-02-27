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
        Real *physical_time_;
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
    };

  protected:
    SingularVariable<Real> *sv_physical_time_;
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
template <typename AlignedBoxPartType, class ConditionFunction>
class BufferEmitterInflowInjectionCK : public BaseLocalDynamics<AlignedBoxPartType>
{
    using CreateRealParticleKernel = typename SpawnRealParticle::ComputingKernel;
    using ConditionKernel = typename ConditionFunction::ComputingKernel;

  public:
    // 1) Move constructor body here, inline:
    BufferEmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer)
        : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
          buffer_(buffer),
          sv_aligned_box_(aligned_box_part.svAlignedBox()),
          create_real_particle_method_(this->particles_),
          rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
          sound_speed_(0.0),
          dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
          dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
          dv_p_(this->particles_->template getVariableByName<Real>("Pressure")),
          dv_buffer_particle_indicator_(this->particles_->template getVariableByName<int>("BufferParticleIndicator")),
          condition_function_(this->particles_),
          dv_previous_surface_indicator_(this->particles_->template getVariableByName<int>("PreviousSurfaceIndicator")),
          sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
          upper_bound_fringe_(0.5 * this->sph_body_.getSPHBodyResolutionRef())

    {
        buffer_.checkParticlesReserved();
        Fluid &fluid_ = DynamicCast<Fluid>(this, this->particles_->getBaseMaterial());
        sound_speed_ = fluid_.getSoundSpeed(rho0_);
    }

    // 2) Inline definition for FinishDynamics constructor:
    class FinishDynamics
    {
      public:
        FinishDynamics(BufferEmitterInflowInjectionCK &encloser)
            : particles_(encloser.particles_), buffer_(encloser.buffer_) {}

        void operator()()
        {
            buffer_.checkEnoughBuffer(*particles_);
        }

      private:
        BaseParticles *particles_;
        ParticleBuffer<Base> &buffer_;
    };

    // 3) Inline definition for UpdateKernel constructor & update(...) if needed:
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : aligned_box_(nullptr),
              create_real_particle_(ex_policy, encloser.create_real_particle_method_),
              rho0_(encloser.rho0_),
              sound_speed_(encloser.sound_speed_),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
              p_(encloser.dv_p_->DelegatedData(ex_policy)),
              buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)),
              condition_(ex_policy, encloser.condition_function_),
              previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
              physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
              upper_bound_fringe_(encloser.upper_bound_fringe_)
        {
            aligned_box_ = encloser.sv_aligned_box_->DelegatedData(ex_policy);
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            if (!aligned_box_->checkInBounds(pos_[index_i]))
            {
                if (aligned_box_->checkUpperBound(pos_[index_i]))
                {
                    if (buffer_particle_indicator_[index_i] == 1)
                    {
                        // if (index_i < this->particles_->TotalRealParticles())
                        {
                            size_t new_particle_index = create_real_particle_(index_i);
                            buffer_particle_indicator_[new_particle_index] = 0;
                            pos_[index_i] = aligned_box_->getUpperPeriodic(pos_[index_i]); // Periodic bounding.
                            p_[index_i] = condition_(index_i, *physical_time_);
                            rho_[index_i] = p_[index_i] / pow(sound_speed_, 2.0) + rho0_;
                            previous_surface_indicator_[index_i] = 1;
                        }
                    }
                }
            }
        }

      private:
        AlignedBox *aligned_box_;
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

    virtual ~BufferEmitterInflowInjectionCK() {};

  protected:
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    SpawnRealParticle create_real_particle_method_;
    Real rho0_;
    Real sound_speed_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
    DiscreteVariable<int> *dv_buffer_particle_indicator_;
    ConditionFunction condition_function_;
    DiscreteVariable<int> *dv_previous_surface_indicator_;
    SingularVariable<Real> *sv_physical_time_;
    Real upper_bound_fringe_;
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

template <typename ExecutionPolicy, class KernelCorrectionType, class PressureConditionFunction>
class PressureBidirectionalConditionCK
{
  public:
    StateDynamics<ExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;

    StateDynamics<ExecutionPolicy,
                  PressureConditionCK<AlignedBoxPartByCell, KernelCorrectionType, PressureConditionFunction>>
        pressure_condition_;

    StateDynamics<execution::SequencedDevicePolicy, BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, PressureConditionFunction>> emitter_injection_;

    StateDynamics<execution::SequencedDevicePolicy, DisposerOutflowDeletionCK> disposer_outflow_deletion_;

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
class NonPrescribedPressure : public BaseStateCondition
{
  public:
    NonPrescribedPressure(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    // The "ComputingKernel" is what your SPHinXsys code will instantiate
    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser) {}

        // Return a fixed or trivial pressure value
        Real operator()(size_t index_i, Real time)
        {
            return p_[index_i]; // Do nothing!
        }
    };
};
template <typename ExecutionPolicy, class KernelCorrectionType, class VelocityConditionFunction>
class VelocityBidirectionalConditionCK
{
  public:
    StateDynamics<ExecutionPolicy, TagBufferParticlesCK> tag_buffer_particles_;

    StateDynamics<ExecutionPolicy, fluid_dynamics::InflowConditionCK<AlignedBoxPartByCell, VelocityConditionFunction>> velocity_condition_;

    StateDynamics<execution::SequencedDevicePolicy, BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, NonPrescribedPressure>> emitter_injection_; // Using NonPrescribedPressure as we do not change pressure

    StateDynamics<execution::SequencedDevicePolicy, DisposerOutflowDeletionCK> disposer_outflow_deletion_;

    VelocityBidirectionalConditionCK(AlignedBoxPartByCell &emitter_by_cell,
                                     ParticleBuffer<Base> &inlet_buffer)
        : tag_buffer_particles_(emitter_by_cell),
          velocity_condition_(emitter_by_cell),
          emitter_injection_(emitter_by_cell, inlet_buffer),
          disposer_outflow_deletion_(emitter_by_cell)
    {
    }

    /// Tag (or flag) particles in the buffer.
    void tagBufferParticles() { tag_buffer_particles_.exec(); }

    /// Apply the velocity condition.
    void applyVelocityCondition() { velocity_condition_.exec(); }

    /// Perform the injection step (for inflow).
    void injectParticles() { emitter_injection_.exec(); }

    /// Perform the deletion step (for outflow).
    void deleteParticles() { disposer_outflow_deletion_.exec(); }
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_H