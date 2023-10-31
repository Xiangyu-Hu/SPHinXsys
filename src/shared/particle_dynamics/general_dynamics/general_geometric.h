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

#include "general_dynamics.h"

namespace SPH
{
/**
 * @class NormalDirectionFromBodyShape
 * @brief normal direction at particles
 */
class NormalDirectionFromBodyShape : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    explicit NormalDirectionFromBodyShape(SPHBody &sph_body);
    virtual ~NormalDirectionFromBodyShape(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Shape &body_shape_;
    StdLargeVec<Vecd> &pos_, &n_, &n0_;
};

/**
 * @class NormalDirectionFromShapeAndOp
 * @brief normal direction at particles
 */
class NormalDirectionFromShapeAndOp : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    explicit NormalDirectionFromShapeAndOp(SPHBody &sph_body, const std::string &shape_name);
    virtual ~NormalDirectionFromShapeAndOp(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    ShapeAndOp *shape_and_op_;
    Shape *shape_;
    const Real switch_sign_;
    StdLargeVec<Vecd> &pos_, &n_, &n0_;
};
} // namespace SPH
#endif // GENERAL_GEOMETRIC_H
