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

#include "base_data_package.h"
#include "base_fluid_dynamics.h"
#include "particle_reserve.h"
#include "sphinxsys_variable_array.h"

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
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    ConditionFunction condition_function_;
};

/**
 * @class EmitterInflowInjectionCK
 * @brief Inject particles into the computational domain.
 * Note that the axis is at the local coordinate and upper bound direction is
 * the local positive direction.
 */
class EmitterInflowInjectionCK : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    EmitterInflowInjectionCK(BodyAlignedBoxByParticle &aligned_box_part, ParticleBuffer<Base> &buffer);
    virtual ~EmitterInflowInjectionCK() {};

    void update(size_t original_index_i, Real dt = 0.0);

  protected:
    std::mutex mutex_switch_to_real_; /**< mutex exclusion for memory conflict */
    Fluid &fluid_;
    UnsignedInt *original_id_;
    UnsignedInt *sorted_id_;
    Vecd *pos_;
    Real *rho_, *p_;
    ParticleBuffer<Base> &buffer_;
    AlignedBoxShape &aligned_box_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_H