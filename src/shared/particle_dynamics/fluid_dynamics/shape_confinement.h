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
 * @file 	shape_confinement.h
 * @brief 	Here, we define the boundary condition classes for fluid dynamics.
 * @details The boundary conditions very often based on different types of buffers.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SHAPE_CONFINEMENT_H
#define SHAPE_CONFINEMENT_H

#include "base_fluid_dynamics.h"
#include "general_constraint.h"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class StaticConfinementDensity
 * @brief static confinement condition for density summation
 */
class StaticConfinementDensity : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    StaticConfinementDensity(NearShapeSurface &near_surface);
    virtual ~StaticConfinementDensity(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real rho0_, inv_sigma0_;
    Real *mass_, *rho_sum_;
    Vecd *pos_;
    LevelSetShape *level_set_shape_;
};

/**
 * @class StaticConfinementIntegration1stHalf
 * @brief static confinement condition for pressure relaxation
 */
class StaticConfinementIntegration1stHalf : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    StaticConfinementIntegration1stHalf(NearShapeSurface &near_surface);
    virtual ~StaticConfinementIntegration1stHalf(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real *rho_, *p_, *mass_;
    Vecd *pos_, *vel_, *force_;
    LevelSetShape *level_set_shape_;
    AcousticRiemannSolver riemann_solver_;
};

/**
 * @class StaticConfinementIntegration2ndHalf
 * @brief static confinement condition for density relaxation
 */
class StaticConfinementIntegration2ndHalf : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    StaticConfinementIntegration2ndHalf(NearShapeSurface &near_surface);
    virtual ~StaticConfinementIntegration2ndHalf(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real *rho_, *p_, *drho_dt_;
    Vecd *pos_, *vel_;
    LevelSetShape *level_set_shape_;
    AcousticRiemannSolver riemann_solver_;
};

/**
 * @class StaticConfinement
 * @brief Static confined boundary condition for complex structures.
 */
class StaticConfinement
{
  public:
    SimpleDynamics<StaticConfinementDensity> density_summation_;
    SimpleDynamics<StaticConfinementIntegration1stHalf> pressure_relaxation_;
    SimpleDynamics<StaticConfinementIntegration2ndHalf> density_relaxation_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;

    StaticConfinement(NearShapeSurface &near_surface);
    virtual ~StaticConfinement(){};
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // SHAPE_CONFINEMENT_H