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
 * @file    general_geometric.h
 * @brief   This is the particle dynamics applicable for all type bodies
 * @author	Xiangyu Hu
 */

#ifndef GENERAL_GEOMETRIC_H
#define GENERAL_GEOMETRIC_H

#include "base_general_dynamics.h"

namespace SPH
{
/**
 * @class NormalDirectionFromBodyShape
 * @brief normal direction at particles
 */
class NormalDirectionFromBodyShape : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit NormalDirectionFromBodyShape(SPHBody &sph_body);
    virtual ~NormalDirectionFromBodyShape(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Shape &initial_shape_;
    StdLargeVec<Vecd> &pos_, &n_, &n0_;
    StdLargeVec<Real> &phi_, &phi0_;
};

/**
 * @class NormalDirectionFromSubShapeAndOp
 * @brief normal direction at particles
 */
class NormalDirectionFromSubShapeAndOp : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit NormalDirectionFromSubShapeAndOp(SPHBody &sph_body, const std::string &shape_name);
    virtual ~NormalDirectionFromSubShapeAndOp(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    SubShapeAndOp *shape_and_op_;
    Shape *shape_;
    const Real switch_sign_;
    StdLargeVec<Vecd> &pos_, &n_, &n0_;
    StdLargeVec<Real> &phi_, &phi0_;
};

class NormalDirectionFromParticles : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NormalDirectionFromParticles(BaseInnerRelation &inner_relation);
    virtual ~NormalDirectionFromParticles(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Shape &initial_shape_;
    StdLargeVec<Vecd> &pos_, &n_, &n0_;
    StdLargeVec<Real> &phi_, &phi0_, &Vol_;
};
} // namespace SPH
#endif // GENERAL_GEOMETRIC_H
