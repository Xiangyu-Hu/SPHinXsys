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

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

template <typename TargetPressure, class ExecutionPolicy = ParallelPolicy>
class BidirectionalBuffer
{
  protected:
    TargetPressure target_pressure_;

    class TagBufferParticles : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        TagBufferParticles(BodyAlignedBoxByCell &aligned_box_part)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              part_id_(aligned_box_part.getPartID()),
              pos_(particles_->getVariableDataByName<Vecd>("Position")),
              aligned_box_(aligned_box_part.getAlignedBoxShape()),
              buffer_particle_indicator_(particles_->registerStateVariable<int>("BufferParticleIndicator"))
        {
            particles_->addVariableToSort<int>("BufferParticleIndicator");
        };
        virtual ~TagBufferParticles(){};

        virtual void update(size_t index_i, Real dt = 0.0)
        {
            buffer_particle_indicator_[index_i] = aligned_box_.checkInBounds(pos_[index_i]) ? part_id_ : 0;
        };

      protected:
        int part_id_;
        Vecd *pos_;
        AlignedBoxShape &aligned_box_;
        int *buffer_particle_indicator_;
    };

    class InjectionAndDeletion : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        InjectionAndDeletion(BodyAlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer,
                             TargetPressure &target_pressure)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              part_id_(aligned_box_part.getPartID()),
              particle_buffer_(particle_buffer),
              aligned_box_(aligned_box_part.getAlignedBoxShape()),
              fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
              pos_(particles_->getVariableDataByName<Vecd>("Position")),
              rho_(particles_->getVariableDataByName<Real>("Density")),
              p_(particles_->getVariableDataByName<Real>("Pressure")),
              previous_surface_indicator_(particles_->getVariableDataByName<int>("PreviousSurfaceIndicator")),
              buffer_particle_indicator_(particles_->getVariableDataByName<int>("BufferParticleIndicator")),
              physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")),
              target_pressure_(target_pressure)
        {
            particle_buffer_.checkParticlesReserved();
        };
        virtual ~InjectionAndDeletion(){};

        void update(size_t index_i, Real dt = 0.0)
        {
            if (aligned_box_.checkUpperBound(pos_[index_i]) &&
                buffer_particle_indicator_[index_i] == part_id_ &&
                index_i < particles_->TotalRealParticles())
            {
                mutex_switch.lock();
                particle_buffer_.checkEnoughBuffer(*particles_);
                size_t new_particle_index = particles_->createRealParticleFrom(index_i);
                buffer_particle_indicator_[new_particle_index] = 0;
                mutex_switch.unlock();

                /** Periodic bounding. */
                pos_[index_i] = aligned_box_.getUpperPeriodic(pos_[index_i]);
                Real sound_speed = fluid_.getSoundSpeed(rho_[index_i]);
                p_[index_i] = target_pressure_(p_[index_i], *physical_time_);
                rho_[index_i] = p_[index_i] / pow(sound_speed, 2.0) + fluid_.ReferenceDensity();
                previous_surface_indicator_[index_i] = 1;
            }

            while (aligned_box_.checkLowerBound(pos_[index_i]) &&
                   buffer_particle_indicator_[index_i] == part_id_ &&
                   index_i < particles_->TotalRealParticles())
            {
                mutex_switch.lock();
                particles_->switchToBufferParticle(index_i);
                mutex_switch.unlock();
            }
        }

      protected:
        int part_id_;
        std::mutex mutex_switch;
        ParticleBuffer<Base> &particle_buffer_;
        AlignedBoxShape &aligned_box_;
        Fluid &fluid_;
        Vecd *pos_;
        Real *rho_, *p_;
        int *previous_surface_indicator_, *buffer_particle_indicator_;
        Real *physical_time_;

      private:
        TargetPressure &target_pressure_;
    };

  public:
    BidirectionalBuffer(BodyAlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer)
        : target_pressure_(*this), tag_buffer_particles(aligned_box_part),
          injection_deletion(aligned_box_part, particle_buffer, target_pressure_){};
    virtual ~BidirectionalBuffer(){};

    SimpleDynamics<TagBufferParticles, ExecutionPolicy> tag_buffer_particles;
    SimpleDynamics<InjectionAndDeletion, ExecutionPolicy> injection_deletion;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BUFFER_H