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
 * @file 	fluid_dynamics_pressure.h
 * @brief 	Here, we define different ways to handle particle_pressure variable in fluid dynamics.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_PRESSURE_H
#define FLUID_DYNAMICS_PRESSURE_H

#include "base_particles.hpp"
#include "weakly_compressible_fluid.h"

namespace SPH
{
class BasePressure
{
  public:
    BasePressure(Fluid &fluid, BaseParticles *base_particles)
        : fluid_(fluid),
          p_(*base_particles->getVariableByName<Real>("Pressure")),
          rho_(*base_particles->getVariableByName<Real>("Density")){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &p_;
    StdLargeVec<Real> &rho_;
};

class PlainPressure : public BasePressure
{
  public:
    PlainPressure(Fluid &fluid, BaseParticles *base_particles)
        : BasePressure(fluid, base_particles){};
    void update(size_t index_i)
    {
        p_[index_i] = fluid_.getPressure(rho_[index_i]);
    };
    Real atInterface(size_t index_i, size_t index_j)
    {
        return 0.5 * (p_[index_i] + p_[index_j]);
    };
};

class BaseCorrectedPressure : public BasePressure
{
  public:
    BaseCorrectedPressure(Fluid &fluid, BaseParticles *base_particles)
        : BasePressure(fluid, base_particles),
          B_(*base_particles->getVariableByName<Matd>("KernelCorrectionMatrix")){};

  protected:
    StdLargeVec<Matd> &B_;
};

class PlainCorrectedPressure : public BaseCorrectedPressure
{
  public:
    PlainCorrectedPressure(Fluid &fluid, BaseParticles *base_particles)
        : BaseCorrectedPressure(fluid, base_particles)
    {
        base_particles->registerVariable(p_B_, "CorrectedPressure");
    };
    void update(size_t index_i)
    {
        p_[index_i] = fluid_.getPressure(rho_[index_i]);
        p_B_[index_i] = B_[index_i] * p_[index_i];
    };
    Matd atInterface(size_t index_i, size_t index_j)
    {
        return 0.5 * (p_B_[index_i] + p_B_[index_j]);
    };

  protected:
    StdLargeVec<Matd> p_B_;
};

class StaggeredCorrectedPressure : public BaseCorrectedPressure
{
  public:
    StaggeredCorrectedPressure(Fluid &fluid, BaseParticles *base_particles)
        : BaseCorrectedPressure(fluid, base_particles){};
    void update(size_t index_i)
    {
        p_[index_i] = fluid_.getPressure(rho_[index_i]);
    };
    Matd atInterface(size_t index_i, size_t index_j)
    {
        return 0.5 * (B_[index_j] * p_[index_i] + B_[index_i] * p_[index_j]);
    };
};
} // namespace SPH
#endif // FLUID_DYNAMICS_PRESSURE_H