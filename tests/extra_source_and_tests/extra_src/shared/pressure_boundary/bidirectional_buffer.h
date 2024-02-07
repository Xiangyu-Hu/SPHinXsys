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
 * @file 	fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BIDIRECTIONAL_BUFFER_H
#define BIDIRECTIONAL_BUFFER_H

#include "sphinxsys.h"

namespace SPH
{
namespace fluid_dynamics
{
class BidirectionalBuffer
{
  protected:
    ConcurrentIndexVector buffer_particle_list;
    virtual Real getTargetPressure(Real dt) = 0;

    template <class ExecutionPolicy>
    class TagBufferParticles : public BaseDynamics<void>, public LocalDynamics, public GeneralDataDelegateSimple 
    {
      private:
        SharedPtrKeeper<Shape> shape_ptr_keeper_;

      public:
        Shape &body_part_shape_;

        TagBufferParticles(ConcurrentIndexVector &buffer_particle_list, RealBody &real_body, SharedPtr<Shape> shape_ptr)
            : BaseDynamics<void>(real_body), LocalDynamics(real_body), GeneralDataDelegateSimple(real_body),
              body_part_shape_(shape_ptr_keeper_.assignRef(shape_ptr)), buffer_particle_list_(buffer_particle_list),
              base_particles_(real_body.getBaseParticles()), 
              buffer_particle_indicator_(*particles_->registerSharedVariable<int>("BufferParticleIndicator"))
        {
            particles_->registerSortableVariable<int>("BufferParticleIndicator");    
        };
        virtual ~TagBufferParticles(){};
    
        virtual void update(size_t index_i, Real dt = 0.0) 
        {
            if (body_part_shape_.checkContain(base_particles_.pos_[index_i]))          
            {
                buffer_particle_list_.push_back(index_i);
                buffer_particle_indicator_[index_i] = 1;               
            }
        };
        
        virtual void exec(Real dt = 0.0) override
        {
            setupDynamics(dt);
            
            particle_for(ExecutionPolicy(), base_particles_.total_real_particles_,
                         [&](size_t index_i)
                         { update(index_i, dt); });
        };

      protected:
        ConcurrentIndexVector &buffer_particle_list_;
        BaseParticles &base_particles_;
        StdLargeVec<int> &buffer_particle_indicator_;  
    };

    template <class ExecutionPolicy>
    class Injection : public BaseDynamics<void>, public LocalDynamics, public GeneralDataDelegateSimple 
    {
      public:
        Injection(ConcurrentIndexVector &buffer_particle_list, RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                  size_t body_buffer_width, int axis_direction, BidirectionalBuffer& buffer) : 
            BaseDynamics<void>(real_body), LocalDynamics(real_body),GeneralDataDelegateSimple(real_body),
            buffer_particle_list_(buffer_particle_list), aligned_box_(*shape_ptr.get()),
            fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), 
            pos_n_(particles_->pos_), rho_n_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), 
            axis_(axis_direction), body_buffer_width_(body_buffer_width),                                                                                                          
            previous_surface_indicator_(*particles_->getVariableByName<int>("PreviousSurfaceIndicator")),                                                                                                                 
            buffer_particle_indicator_(*particles_->getVariableByName<int>("BufferParticleIndicator")),
            buffer_(buffer)           
        {
            size_t total_body_buffer_particles = 1000.0 * body_buffer_width_;
            particles_->addBufferParticles(total_body_buffer_particles);
            real_body.allocateConfigurationMemoriesForBufferParticles();
        };
        virtual ~Injection(){};

        /** This class is only implemented in sequential due to memory conflicts. */
        virtual void update(size_t index_i, Real dt = 0.0);

        virtual void exec(Real dt = 0.0) override
        {
            setupDynamics(dt);

            particle_for(ExecutionPolicy(), buffer_particle_list_.size(),
                         [&](size_t index_i)
                         { update(index_i, dt); });

            buffer_particle_list_.clear();
        };

      protected:
        std::mutex mutex_switch_to_real_;
        ConcurrentIndexVector &buffer_particle_list_;
        AlignedBoxShape &aligned_box_;
        Fluid &fluid_;
        StdLargeVec<Vecd> &pos_n_;
        StdLargeVec<Real> &rho_n_, &p_;
        const int axis_; 
        size_t body_buffer_width_;
        StdLargeVec<int> &previous_surface_indicator_, &buffer_particle_indicator_;

      private:
        BidirectionalBuffer &buffer_;
    };

  public:
    BidirectionalBuffer(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                        size_t body_buffer_width, int axis_direction) : 
        tag_buffer_particles(this->buffer_particle_list, real_body, shape_ptr),                                                                                                                     
        injection(this->buffer_particle_list, real_body, shape_ptr, body_buffer_width, axis_direction, *this)
    {
        buffer_particle_list.clear();
    };
    virtual ~BidirectionalBuffer(){};
    
    TagBufferParticles<execution::ParallelPolicy> tag_buffer_particles;
    Injection<execution::ParallelPolicy> injection;   
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BUFFER_H