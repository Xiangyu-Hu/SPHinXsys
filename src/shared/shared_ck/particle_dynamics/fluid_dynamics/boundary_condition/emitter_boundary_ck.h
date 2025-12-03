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
 * @file 	emitter_boundary_ck.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef EMITTER_BOUNDARY_CK_H
#define EMITTER_BOUNDARY_CK_H

#include "base_body_part.h"
#include "base_data_type_package.h"
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
class EmitterInflowConditionCK;

template <class AlignedBoxPartType, class ConditionFunction>
class EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>
    : public BaseLocalDynamics<AlignedBoxPartType>
{
    using ConditionKernel = typename ConditionFunction::ComputingKernel;

  public:
    EmitterInflowConditionCK(AlignedBoxPartType &aligned_box_part);

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
    using SpawnRealParticleKernel = typename SpawnRealParticle::ComputingKernel;

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
        SpawnRealParticleKernel spawn_real_particle_;
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
    SpawnRealParticle spawn_real_particle_method_;
    Real rho0_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EMITTER_BOUNDARY_CK_H