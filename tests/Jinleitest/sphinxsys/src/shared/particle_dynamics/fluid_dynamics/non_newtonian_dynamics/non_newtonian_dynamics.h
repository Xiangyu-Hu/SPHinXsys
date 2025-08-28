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
 * @file non_newtonian_dynamics.h
 * @brief Here, we define the time integration algorithm classes for non-newtonian fluids.
 * @details We consider here weakly compressible fluids.
 * @author	Theodor Hennings and Xiangyu Hu
 */

#ifndef NON_NEWTONIAN_DYNAMICS_H
#define NON_NEWTONIAN_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "fluid_integration.hpp"
#include "viscosity.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class Oldroyd_BIntegration1stHalf;

template <>
class Oldroyd_BIntegration1stHalf<Inner<>> : public Integration1stHalfInnerRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Matd *tau_, *dtau_dt_;
};

using Integration1stHalfContactWallRiemann =
    Integration1stHalf<Contact<Wall>, AcousticRiemannSolver, NoKernelCorrection>;

template <>
class Oldroyd_BIntegration1stHalf<Contact<Wall>> : public Integration1stHalfContactWallRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Matd *tau_;
};

template <typename... InteractionTypes>
class Oldroyd_BIntegration2ndHalf;

template <>
class Oldroyd_BIntegration2ndHalf<Inner<>> : public Integration2ndHalfInnerRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *vel_grad_, *tau_, *dtau_dt_;
    Real mu_p_, lambda_;
};

using Integration2ndHalfContactWallRiemann =
    Integration2ndHalf<Contact<Wall>, AcousticRiemannSolver>;

template <>
class Oldroyd_BIntegration2ndHalf<Contact<Wall>> : public Integration2ndHalfContactWallRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseContactRelation &wall_contact_relation)
        : Integration2ndHalfContactWallRiemann(wall_contact_relation){};
    virtual ~Oldroyd_BIntegration2ndHalf(){};
};

using Oldroyd_BIntegration1stHalfWithWall = ComplexInteraction<Oldroyd_BIntegration1stHalf<Inner<>, Contact<Wall>>>;
using Oldroyd_BIntegration2ndHalfWithWall = ComplexInteraction<Oldroyd_BIntegration2ndHalf<Inner<>, Contact<Wall>>>;

/**
 * @class SRDViscousTimeStepSize
 * @brief Computing the viscous time step size using the SRD viscosity
 */
class SRDViscousTimeStepSize : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit SRDViscousTimeStepSize(SPHBody &sph_body, Real diffusionCFL = 0.125);
    virtual ~SRDViscousTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Real smoothing_length_;
    Real *rho_;
    Real *mu_srd_;
    Real diffusionCFL;
    Real max_viscosity = 1e-12;
};

class ShearRateDependentViscosity : public LocalDynamics
{
  public:
    explicit ShearRateDependentViscosity(SPHBody &sph_body);
    virtual ~ShearRateDependentViscosity(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *vel_grad_;
    GeneralizedNewtonianViscosity &generalized_viscosity_;
    Real *mu_srd_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // NON_NEWTONIAN_DYNAMICS_H