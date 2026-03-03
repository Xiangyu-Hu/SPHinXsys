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
 * @file 	fluid_boundary_state.h
 * @brief 	tbd
 * @author	Xiangyu Hu
 */

#ifndef FLUID_BOUNDARY_STATE_H
#define FLUID_BOUNDARY_STATE_H

#include "base_data_type_package.h"
#include "base_particles.hpp"
#include "weakly_compressible_fluid.h"
namespace SPH
{

class BaseStateCondition
{
  public:
    BaseStateCondition(BaseParticles *particles);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

      protected:
        Vecd *pos_, *vel_;
        Real *p_, *rho_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
    DiscreteVariable<Real> *dv_p_, *dv_rho_;
};

template <class FluidType = WeaklyCompressibleFluid>
struct PressurePrescribed
{
    typedef FluidType Fluid;
    Real target_pressure_;
    PressurePrescribed(Real target_pressure) : target_pressure_(target_pressure) {};
    Real getPressure(const Real &input_pressure, Real time) { return target_pressure_; };
    Real getAxisVelocity(const Vecd &input_position, const Real &input_axis_velocity, Real time)
    {
        return input_axis_velocity;
    };
};

template <class FluidType = WeaklyCompressibleFluid>
struct VelocityPrescribed
{
    typedef FluidType Fluid;
    Real getPressure(const Real &input_pressure, Real time) { return input_pressure; };
    // Real operator()(const Vecd &input_position, const Real &input_axis_velocity, Real time)
    // to be implemented in derived class
};

} // namespace SPH
#endif // FLUID_BOUNDARY_STATE_H