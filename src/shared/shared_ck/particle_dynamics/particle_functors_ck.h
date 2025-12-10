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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */

/**
 * @file    particle_functors.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_FUNCTORS_CK_H
#define PARTICLE_FUNCTORS_CK_H

namespace SPH
{
/*
 * ParticleScopeTypeCK determines the scope of particles for which transport
 * velocity correction is applied. It is a template class specialized for
 * different particle types, defining rules for inclusion in the computation.
 */
template <typename ScopeType>
class ParticleScopeTypeCK;
//-------------------------------------------------------------------------------------------------
// 1) Specialization for AllParticles
//-------------------------------------------------------------------------------------------------
template <>
class ParticleScopeTypeCK<AllParticles> : public WithinScope
{
  public:
    // Constructor
    explicit ParticleScopeTypeCK(BaseParticles *particles)
        : WithinScope()
    {
    }

    // Nested functor-like class:
    class ComputingKernel
    {
      public:
        // The signature typically follows the style of other SPH "ComputingKernel" constructors
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<AllParticles> &encloser,
                        ComputingKernelType &computing_kernel){};

        constexpr bool operator()(size_t /*index_i*/) const
        {
            return true; // Always in scope
        };
    };
};

//-------------------------------------------------------------------------------------------------
// 2) Specialization for BulkParticles (which is typically IndicatedParticles<0> in SPH code)
//    i.e. we only want particles that have indicator_[i] == 0
//-------------------------------------------------------------------------------------------------
template <>
class ParticleScopeTypeCK<BulkParticles> : public WithinScope
{
  public:
    explicit ParticleScopeTypeCK(BaseParticles *particles)
        : WithinScope(),
          dv_indicator_(particles->getVariableByName<int>("Indicator"))
    {
    }

    // Nested class implementing the boolean check
    class ComputingKernel
    {
      public:
        // Typically, we pass the "encloser" object to get the pointer from dv_indicator_
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<BulkParticles> &encloser,
                        ComputingKernelType &computing_kernel)
            : indicator_(encloser.dv_indicator_->DelegatedData(ex_policy))
        {
        }

        bool operator()(size_t index_i) const
        {
            return (indicator_[index_i] == 0);
        }

      protected:
        int *indicator_;
    };

  protected:
    DiscreteVariable<int> *dv_indicator_;
};

//-------------------------------------------------------------------------------------------------
// 3) Specialization for NotIndicatedParticles (which is typically "indicator != 0")
//-------------------------------------------------------------------------------------------------
template <int INDICATOR>
class ParticleScopeTypeCK<NotIndicatedParticles<INDICATOR>> : public WithinScope
{
  public:
    explicit ParticleScopeTypeCK(BaseParticles *particles)
        : WithinScope(),
          dv_indicator_(particles->getVariableByName<int>("Indicator"))
    {
    }

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<NotIndicatedParticles<INDICATOR>> &encloser,
                        ComputingKernelType &computing_kernel)
            : indicator_(encloser.dv_indicator_->DelegatedData(ex_policy))
        {
        }

        bool operator()(size_t index_i) const
        {
            return (indicator_[index_i] != INDICATOR);
        }

      protected:
        int *indicator_;
    };

  protected:
    DiscreteVariable<int> *dv_indicator_;
};

class ExcludeBufferParticles;
template <>
class ParticleScopeTypeCK<ExcludeBufferParticles> : public WithinScope
{
  public:
    explicit ParticleScopeTypeCK(BaseParticles *particles)
        : WithinScope(),
          dv_buffer_particles_indicator_(particles->registerStateVariable<int>("BufferIndicator"))
    {
    }

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<ExcludeBufferParticles> &encloser,
                        ComputingKernelType &computing_kernel)
            : buffer_particles_indicator_(encloser.dv_buffer_particles_indicator_->DelegatedData(ex_policy))
        {
        }

        bool operator()(size_t index_i) const
        {
            return (buffer_particles_indicator_[index_i] == 0);
        }

      protected:
        int *buffer_particles_indicator_;
    };

  protected:
    DiscreteVariable<int> *dv_buffer_particles_indicator_;
};

} // namespace SPH
#endif // PARTICLE_FUNCTORS_CK_H
