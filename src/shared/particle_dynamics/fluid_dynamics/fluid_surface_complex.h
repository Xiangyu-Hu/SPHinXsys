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
 * @file 	fluid_surface_complex.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_SURFACE_COMPLEX_H
#define FLUID_SURFACE_COMPLEX_H

#include "fluid_dynamics_complex.hpp"
#include "fluid_surface_inner.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FreeSurfaceIndicationContact
 * @brief indicate the particles near the free fluid surface.
 */
class FreeSurfaceIndicationContact : public FreeSurfaceIndication<FluidContactData>
{
  public:
    FreeSurfaceIndicationContact(BaseContactRelation &contact_relation)
        : FreeSurfaceIndication<FluidContactData>(contact_relation){};
    virtual ~FreeSurfaceIndicationContact(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ColorFunctionGradientContact
 * @brief indicate the particles near the free fluid surface.
 */
class ColorFunctionGradientContact : public ColorFunctionGradient<FluidContactData>
{
  public:
    ColorFunctionGradientContact(BaseContactRelation &contact_relation)
        : ColorFunctionGradient<FluidContactData>(contact_relation){};
    virtual ~ColorFunctionGradientContact(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class 	SurfaceNormWithWall
 * @brief  Modify surface norm when contact with wall
 */
class SurfaceNormWithWall : public LocalDynamics, public FSIContactData
{
  public:
    SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle);
    virtual ~SurfaceNormWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real contact_angle_;
    Real smoothing_length_;
    Real particle_spacing_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> &surface_norm_;
    StdLargeVec<Real> &pos_div_;
    StdVec<StdLargeVec<Vecd> *> wall_n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_SURFACE_COMPLEX_H
