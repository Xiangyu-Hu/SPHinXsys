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
 * @file 	fluid_time_step_ck.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_TIME_STEP_CK_H
#define FLUID_TIME_STEP_CK_H

#include "base_fluid_dynamics.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class FluidType = WeaklyCompressibleFluid>
class AcousticTimeStepCK : public LocalDynamicsReduce<ReduceMax>
{
    using EosKernel = typename FluidType::EosKernel;

  public:
    explicit AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepCK() {};

    class FinishDynamics
    {
        Real h_min_, acousticCFL_;

      public:
        using OutputType = Real;
        FinishDynamics(AcousticTimeStepCK<FluidType> &encloser);
        Real Result(Real reduced_value);
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy>
        ReduceKernel(const ExecutionPolicy &ex_policy, AcousticTimeStepCK<FluidType> &encloser);

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            return eos_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
        };

      protected:
        EosKernel eos_;
        Real *rho_, *p_;
        Vecd *vel_;
        Real h_min_;
    };

  protected:
    FluidType &fluid_;
    DiscreteVariable<Real> *dv_rho_, *dv_p_;
    DiscreteVariable<Vecd> *dv_vel_;
    Real h_min_;
    Real acousticCFL_;
};

class AdvectionTimeStepCK : public LocalDynamicsReduce<ReduceMax>
{
  public:
    AdvectionTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepCK() {};

    class FinishDynamics
    {
        Real h_min_;
        Real speed_ref_, advectionCFL_;

      public:
        using OutputType = Real;
        FinishDynamics(AdvectionTimeStepCK &encloser);
        Real Result(Real reduced_value);
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy>
        ReduceKernel(const ExecutionPolicy &ex_policy, AdvectionTimeStepCK &encloser)
            : h_min_(encloser.h_min_),
              mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              force_(encloser.dv_force_->DelegatedData(ex_policy)),
              force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)){};

        Real reduce(size_t index_i, Real dt)
        {
            Real acceleration_scale =
                4.0 * h_min_ * (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
            return SMAX(vel_[index_i].squaredNorm(), acceleration_scale);
        };

      protected:
        Real h_min_;
        Real *mass_;
        Vecd *vel_, *force_, *force_prior_;
    };

  protected:
    Real h_min_;
    Real speed_ref_, advectionCFL_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_force_, *dv_force_prior_;
};

/**
 * @class AdvectionViscousTimeStepCK
 * @brief Computing the advection time step size
 */
class AdvectionViscousTimeStepCK : public AdvectionTimeStepCK
{
  public:
    AdvectionViscousTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionViscousTimeStepCK() {};
};

class AdvectionStepSetup : public LocalDynamics
{
  public:
    explicit AdvectionStepSetup(SPHBody &sph_body);
    virtual ~AdvectionStepSetup() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser);

        void update(size_t index_i, Real dt = 0.0)
        {
            Vol_[index_i] = mass_[index_i] / rho_[index_i];
            dpos_[index_i] = Vecd::Zero();
        };

      protected:
        Real *Vol_, *mass_, *rho_;
        Vecd *dpos_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_, *dv_mass_, *dv_rho_;
    DiscreteVariable<Vecd> *dv_dpos_;
};

class UpdateParticlePosition : public LocalDynamics
{
  public:
    explicit UpdateParticlePosition(SPHBody &sph_body);
    virtual ~UpdateParticlePosition() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, UpdateParticlePosition &encloser);

        void update(size_t index_i, Real dt = 0.0)
        {
            pos_[index_i] += dpos_[index_i];
        };

      protected:
        Vecd *pos_, *dpos_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_dpos_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_CK_H
