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
 * @file sdf_primitive.h
 * @brief All signed distance function (sdf) primitives use local coordinates.
 * All rotation is around the x-axis.
 * @details Here, we only give popular primitives,
 * For more signed distance function, please check the website:
 * https://iquilezles.org/articles/.
 * @author	Xiangyu Hu
 */

#ifndef SDF_PRIMITIVE_H
#define SDF_PRIMITIVE_H

#include "base_data_type.h"
#include "geometric_primitive.h"
#include "scalar_functions.h"

namespace SPH
{
class SDFBall
{
    Real radius_;

  public:
    explicit SDFBall(Real radius) : radius_(radius) {}
    virtual ~SDFBall() {}
    void setParameters(Real radius) { radius_ = radius; }
    Real operator()(const Vec3d &point) const { return point.norm() - radius_; }
    BoundingBox3d findBounds() const { return BoundingBox3d(Vec3d::Constant(radius_)); }
};

class SDFBox
{
    Vec3d halfsize_;

  public:
    explicit SDFBox(const Vec3d &halfsize) : halfsize_(halfsize) {}
    void setParameters(const Vec3d &halfsize) { halfsize_ = halfsize; }
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const { return BoundingBox3d(halfsize_); }
};

class SDFCylinder
{
    Real halflength_, radius_;

  public:
    explicit SDFCylinder(Real halflength, Real radius) : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCapsule
{
    Real halflength_, radius_;

  public:
    explicit SDFCapsule(Real halflength, Real radius) : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCone
{
    Real halfheight_, radius_;

  public:
    explicit SDFCone(Real halfheight, Real radius) : halfheight_(halfheight), radius_(radius) {}
    void setParameters(Real halfheight, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCappedCone
{
    Real halfheight_, radius1_, radius2_;

  public:
    explicit SDFCappedCone(Real halfheight, Real radius1, Real radius2)
        : halfheight_(halfheight), radius1_(radius1), radius2_(radius2) {}
    void setParameters(Real halfheight, Real radius1, Real radius2);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};
} // namespace SPH

#endif // SDF_PRIMITIVE_H
