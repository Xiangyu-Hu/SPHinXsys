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
 * @file    interaction_ck.h
 * @brief 	This is for the base classes of erosion particle dynamics, which describe the
 * 			erosion dynamics of a particle and it neighbors.
 * @author	Shuang Li and Xiangyu Hu
 */

#ifndef EROSION_INTERACTION_CK_H
#define EROSION_INTERACTION_CK_H

#include "base_local_dynamics.h"
#include "neighborhood_ck.h"
#include "relation_ck.hpp"
#include "interaction_ck.hpp"

namespace SPH
{
template <>
class Interaction<Soil>
{
  public:
    template <class SoilContactRelationType>
    Interaction(SoilContactRelationType &soil_contact_relation);
    virtual ~Interaction() {};

  protected:
    StdVec<DiscreteVariable<Vecd> *> dv_soil_n_, dv_soil_vel_;
    StdVec<DiscreteVariable<Real> *> dv_soil_Vol_, dv_soil_p_;
};

template <>
class Interaction<Fluid>
{
  public:
    template <class FluidContactRelationType>
    Interaction(FluidContactRelationType &fluid_contact_relation);
    virtual ~Interaction() {};

  protected:
    StdVec<DiscreteVariable<Vecd> *> dv_fluid_vel_;
    StdVec<DiscreteVariable<Real> *> dv_fluid_Vol_, dv_fluid_p_;
};
} // namespace SPH
#endif // EROSION_INTERACTION_CK_H
