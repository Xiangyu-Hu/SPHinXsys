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
#include "particle_reserve.h"
#include "sphinxsys_variable_array.h"
#include "fluid_boundary_state.hpp"

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
        AlignedBox *aligned_box_;
        ConditionKernel condition_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    ConditionFunction condition_function_;
};

class EmitterInflowInjectionCK : public BaseLocalDynamics<AlignedBoxPartByParticle>
{
  public:
    EmitterInflowInjectionCK(AlignedBoxPartByParticle &aligned_box_part, ParticleBuffer<Base> &buffer);
    virtual ~EmitterInflowInjectionCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        Real rho0_;
        VariableDataArrays copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK> copy_particle_state_;
        Vecd *pos_;
        Real *rho_, *p_;
        UnsignedInt *total_real_particles_;
    };

  protected:
    ParticleBuffer<Base> &buffer_;
    SingularVariable<AlignedBox> *sv_aligned_box_;
    Real rho0_;
    DiscreteVariableArrays copyable_states_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;

};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_H