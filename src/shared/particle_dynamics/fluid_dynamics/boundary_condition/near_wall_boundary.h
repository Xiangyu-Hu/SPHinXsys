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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file near_wall_boundary.h
 * @brief To evaluate the distance from the fluid particle to wall surface.
 * @details The vector is pointing from the wall surface point to the particle.
 * Note that if the fluid particle is too close to the wall surface
 * (less then 1/4 of particle spacing), its position will be corrected to
 * half of the particle spacing.
 * @author Xiangyu Hu
 */

#ifndef NEAR_WALL_BOUNDARY_H
#define NEAR_WALL_BOUNDARY_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
class NearWallDistance : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit NearWallDistance(BaseContactRelation &wall_contact_relation);
    virtual ~NearWallDistance() {};

  protected:
    Real spacing_ref_, distance_default_;
    Vecd *pos_;
    StdVec<Vecd *> wall_pos_, wall_n_;
    StdVec<Real *> wall_phi_;

    void evaluateDistanceAndNormal(size_t index_i, Vecd &distance, Vecd &normal);
};

class DistanceFromWall : public NearWallDistance
{
  public:
    explicit DistanceFromWall(BaseContactRelation &wall_contact_relation);
    virtual ~DistanceFromWall() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *distance_from_wall_;
};

class BoundingFromWall : public NearWallDistance
{
  public:
    explicit BoundingFromWall(BaseContactRelation &wall_contact_relation);
    virtual ~BoundingFromWall() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real distance_min_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // NEAR_WALL_BOUNDARY_H
