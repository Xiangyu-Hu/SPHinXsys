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
 * @file 	geometric_element_and_shape_3d.h
 * @brief tbd.
 * @author Xiangyu Hu
 */

#ifndef GEOMETRIC_ELEMENT_AND_SHAPE_3D_H
#define GEOMETRIC_ELEMENT_AND_SHAPE_3D_H

#include "data_type.h"
#include "transform_geometry.h"

namespace SPH
{
class GeometricCylinder
{
  public:
    explicit GeometricCylinder(Real radius, Real halflength);
    ~GeometricCylinder(){};

    bool checkContain(const Vec3d &probe_point)
    {
        if (ABS(probe_point[0]) > halflength_)
            return false;
        return probe_point.tail(Dimensions - 1).norm() <= radius_;
    };

    Vec3d findClosestPoint(const Vec3d &probe_point);
    BoundingBox3d findBounds();

  protected:
    Real radius_;
    Real halflength_;
};

using TransformGeometryCylinder = TransformGeometry<GeometricCylinder>;

class GeometricShapeCylinder : public TransformShape<GeometricCylinder>
{
  public:
    GeometricShapeCylinder(const Transform &transform, Real radius, Real halflength,
                           const std::string &name = "GeometricShapeCylinder");
    virtual ~GeometricShapeCylinder(){};
};
} // namespace SPH
#endif // GEOMETRIC_ELEMENT_AND_SHAPE_3D_H
