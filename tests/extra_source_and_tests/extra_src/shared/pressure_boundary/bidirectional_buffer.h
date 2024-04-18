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
 * @file 	bidirectional_buffer.h
 * @brief 	Here, we define the algorithm classes for bidirectiontal buffer.
 * @details The buffer particle index is periodically updated at each time step.
            The bidirectional buffer can serve for unidirectional, bidirectional and mixed flows.
 * @author	Shuoguo Zhang and Xiangyu Hu
 */

#ifndef BIDIRECTIONAL_BUFFER_H
#define BIDIRECTIONAL_BUFFER_H

#include "sphinxsys.h"

namespace SPH
{
namespace fluid_dynamics
{
struct NonPrescribedPressure
{
    template <class BoundaryConditionType>
    NonPrescribedPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

template <typename TargetPressure, class ExecutionPolicy = ParallelPolicy>
class BidirectionalBuffer
{
  protected:
    TargetPressure target_pressure_;

    class TagBufferParticles : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
    {
      public:
        TagBufferParticles(BodyAlignedBoxByCell &aligned_box_part, int axis)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part), FluidDataSimple(sph_body_),
              pos_(*base_particles_.getVariableByName<Vecd>("Position")), aligned_box_(aligned_box_part.aligned_box_), axis_(axis),
              buffer_particle_indicator_(*particles_->registerSharedVariable<int>("BufferParticleIndicator"))
        {
            particles_->addVariableToSort<int>("BufferParticleIndicator");
        };
        virtual ~TagBufferParticles(){};

        virtual void update(size_t index_i, Real dt = 0.0)
        {
            buffer_particle_indicator_[index_i] = aligned_box_.checkInBounds(axis_, pos_[index_i]) ? 1 : 0;
        };

      protected:
        StdLargeVec<Vecd> &pos_;
        AlignedBoxShape &aligned_box_;
        const int axis_;
        StdLargeVec<int> &buffer_particle_indicator_;
    };

    class Injection : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
    {
      public:
        Injection(BodyAlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer,
                  int axis, TargetPressure &target_pressure)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part), FluidDataSimple(sph_body_),
              axis_(axis), particle_buffer_(particle_buffer), aligned_box_(aligned_box_part.aligned_box_),
              fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
              pos_n_(*particles_->getVariableByName<Vecd>("Position")),
              rho_n_(*particles_->getVariableByName<Real>("Density")),
              p_(*particles_->getVariableByName<Real>("Pressure")),
              previous_surface_indicator_(*particles_->getVariableByName<int>("PreviousSurfaceIndicator")),
              buffer_particle_indicator_(*particles_->getVariableByName<int>("BufferParticleIndicator")),
              target_pressure_(target_pressure)
        {
            particle_buffer_.checkParticlesReserved();
        };
        virtual ~Injection(){};

        void update(size_t index_i, Real dt = 0.0)
        {
            if (aligned_box_.checkUpperBound(axis_, pos_n_[index_i]) && buffer_particle_indicator_[index_i] == 1)
            {
                mutex_switch_to_real_.lock();
                particle_buffer_.checkEnoughBuffer(*particles_);
                /** Buffer Particle state copied from real particle. */
                particles_->copyFromAnotherParticle(particles_->total_real_particles_, index_i);
                /** Realize the buffer particle by increasing the number of real particle in the body.  */
                particles_->total_real_particles_ += 1;
                mutex_switch_to_real_.unlock();

                /** Periodic bounding. */
                pos_n_[index_i] = aligned_box_.getUpperPeriodic(axis_, pos_n_[index_i]);
                Real sound_speed = fluid_.getSoundSpeed(rho_n_[index_i]);
                p_[index_i] = target_pressure_(p_[index_i]);
                rho_n_[index_i] = p_[index_i] / pow(sound_speed, 2.0) + fluid_.ReferenceDensity();
                previous_surface_indicator_[index_i] = 1;
            }
        }

      protected:
        std::mutex mutex_switch_to_real_;
        const int axis_;
        ParticleBuffer<Base> &particle_buffer_;
        AlignedBoxShape &aligned_box_;
        Fluid &fluid_;
        StdLargeVec<Vecd> &pos_n_;
        StdLargeVec<Real> &rho_n_, &p_;
        StdLargeVec<int> &previous_surface_indicator_, &buffer_particle_indicator_;

      private:
        TargetPressure &target_pressure_;
    };

  public:
    BidirectionalBuffer(BodyAlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer, int axis)
        : target_pressure_(*this), tag_buffer_particles(aligned_box_part, axis),
          injection(aligned_box_part, particle_buffer, axis, target_pressure_){};
    virtual ~BidirectionalBuffer(){};

    SimpleDynamics<TagBufferParticles, ExecutionPolicy> tag_buffer_particles;
    SimpleDynamics<Injection, ExecutionPolicy> injection;
};

using NonPrescribedPressureBidirectionalBuffer = BidirectionalBuffer<NonPrescribedPressure>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BUFFER_H