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
 * @file surface_tension.h
 * @brief A momentum-conservative formulation for surface tension is used here
 * to reach a long-term stable simulation.
 * @details The zero-surface-energy modes in surface tension are identified and 
 * sovlved by introducting a penalty force. The method can be used to simulate 
 * surface tension in multiphase flows with extra high Reynolds numbers and Weber numbers.
 * @author	Shuaihao Zhang and Xiangyu Hu
 */

#ifndef SURFACE_TENSION_H
#define SURFACE_TENSION_H

#include "base_fluid_dynamics.h"
#include "force_prior.h"

namespace SPH
{
namespace fluid_dynamics
{
class SurfaceTensionStress : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit SurfaceTensionStress(BaseContactRelation &contact_relation, Real surface_tension_coeff);
    virtual ~SurfaceTensionStress(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *color_gradient_, *norm_direction_;
    Matd *surface_tension_stress_;
    StdVec<Real> contact_fraction_;
    StdVec<Real *> contact_Vol_;
    Real &surface_tension_coeff_;
};

template <typename... T>
class SurfaceStressForce;

template <class DataDelegationType>
class SurfaceStressForce<DataDelegationType>
    : public ForcePrior, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit SurfaceStressForce(BaseRelationType &base_relation);
    virtual ~SurfaceStressForce(){};

  protected:
    Real *rho_, *mass_, *Vol_;
    Vecd *color_gradient_, *surface_tension_force_, *norm_direction_;
    Matd *surface_tension_stress_;
    Real &surface_tension_coeff_;
};

template <>
class SurfaceStressForce<Inner<>> : public SurfaceStressForce<DataDelegateInner>
{
  public:
    explicit SurfaceStressForce(BaseInnerRelation &inner_relation, Real hourglass_control_coeff = 4.5);
    template <typename BodyRelationType, typename FirstArg>
    explicit SurfaceStressForce(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : SurfaceStressForce(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~SurfaceStressForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real hourglass_control_coeff_;
};

template <>
class SurfaceStressForce<Contact<>> : public SurfaceStressForce<DataDelegateContact>
{
  public:
    explicit SurfaceStressForce(BaseContactRelation &contact_relation, Real hourglass_control_coeff = 4.5);
    template <typename BodyRelationType, typename FirstArg>
    explicit SurfaceStressForce(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : SurfaceStressForce(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~SurfaceStressForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_color_gradient_, contact_norm_direction_;
    StdVec<Matd *> contact_surface_tension_stress_;
    StdVec<Real> contact_surface_tension_, contact_fraction_;

  protected:
    Real hourglass_control_coeff_;
};

using SurfaceStressForceComplex = ComplexInteraction<SurfaceStressForce<Inner<>, Contact<>>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // SURFACE_TENSION_H
