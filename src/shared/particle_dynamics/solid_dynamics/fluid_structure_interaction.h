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
 * @file 	fluid_structure_interaction.h
 * @brief 	Here, we define the algorithm classes for fluid structure interaction.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_STRUCTURE_INTERACTION_H
#define FLUID_STRUCTURE_INTERACTION_H

#include "all_particle_dynamics.h"
#include "base_material.h"
#include "elastic_dynamics.h"
#include "force_prior.hpp"
#include "riemann_solver.h"

namespace SPH
{
namespace solid_dynamics
{
/**
 * @class BaseForceFromFluid
 * @brief Base class for computing the forces from the fluid
 */
class BaseForceFromFluid : public ForcePrior, public DataDelegateContact
{
  public:
    explicit BaseForceFromFluid(BaseContactRelation &contact_relation, const std::string &force_name);
    virtual ~BaseForceFromFluid() {};
    Vecd *getForceFromFluid() { return force_from_fluid_; };

  protected:
    Solid &solid_;
    Real *Vol_;
    StdVec<Fluid *> contact_fluids_;
    Vecd *force_from_fluid_;
};

/**
 * @class ViscousForceFromFluid
 * @brief Computing the viscous force from the fluid
 */
class ViscousForceFromFluid : public BaseForceFromFluid
{
  public:
    explicit ViscousForceFromFluid(BaseContactRelation &contact_relation);
    virtual ~ViscousForceFromFluid() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vel_ave_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real> mu_;
    StdVec<Real> smoothing_length_;
};

/**
 * @class PressureForceFromFluid
 * @brief Template class fro computing the pressure force from the fluid with different Riemann solvers.
 * The pressure force is added on the viscous force of the latter is computed.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
template <class FluidIntegration2ndHalfType>
class PressureForceFromFluid : public BaseForceFromFluid
{
    using RiemannSolverType = typename FluidIntegration2ndHalfType::RiemannSolver;

  public:
    explicit PressureForceFromFluid(BaseContactRelation &contact_relation);
    virtual ~PressureForceFromFluid() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vel_ave_, *acc_ave_, *n_;
    StdVec<Real *> contact_rho_, contact_mass_, contact_p_, contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_force_prior_;
    StdVec<RiemannSolverType> riemann_solvers_;
};

/**
 * @class InitializeDisplacement
 * @brief initialize the displacement for computing average velocity.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
class InitializeDisplacement : public LocalDynamics
{
  protected:
    Vecd *pos_, *pos_temp_;

  public:
    explicit InitializeDisplacement(SPHBody &sph_body);
    virtual ~InitializeDisplacement() {};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class UpdateAverageVelocityAndAcceleration
 * @brief Computing average velocity.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
class UpdateAverageVelocityAndAcceleration : public LocalDynamics
{
  protected:
    Vecd *pos_, *pos_temp_, *vel_ave_, *acc_ave_;

  public:
    explicit UpdateAverageVelocityAndAcceleration(SPHBody &sph_body);
    virtual ~UpdateAverageVelocityAndAcceleration() {};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class AverageVelocityAndAcceleration
 * @brief Impose force matching between fluid and solid dynamics.
 * Note that the fluid time step should be larger than that of solid time step.
 * Otherwise numerical instability may occur.
 */
class AverageVelocityAndAcceleration
{
  public:
    SimpleDynamics<InitializeDisplacement> initialize_displacement_;
    SimpleDynamics<UpdateAverageVelocityAndAcceleration> update_averages_;

    explicit AverageVelocityAndAcceleration(SPHBody &sph_body);
    ~AverageVelocityAndAcceleration() {};
};
} // namespace solid_dynamics
} // namespace SPH
#endif // FLUID_STRUCTURE_INTERACTION_H