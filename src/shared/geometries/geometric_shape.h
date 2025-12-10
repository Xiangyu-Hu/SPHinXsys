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
 * @file 	geometric_shape.h
 * @brief 	Here, we define shapes represented directly by geometric elements.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GEOMETRIC_SHAPE_H
#define GEOMETRIC_SHAPE_H

#include "base_geometry.h"
#include "geometric_element.h"
#include "transform_geometry.h"

namespace SPH
{
class GeometricShapeBox : public TransformShape<GeometricBox>
{
  public:
    GeometricShapeBox(const Transform &transform, const Vecd &halfsize,
                      const std::string &name = "GeometricShapeBox");
    explicit GeometricShapeBox(const BoundingBoxd &bounding_box,
                               const std::string &name = "GeometricShapeBox");
    GeometricShapeBox(const Vecd &lower_bound, const Vecd &upper_bound,
                      const std::string &name = "GeometricShapeBox");
    virtual ~GeometricShapeBox() {};
};

class GeometricShapeBall : public GeometricBall, public Shape
{
    Vecd center_;

  public:
    explicit GeometricShapeBall(const Vecd &center, Real radius,
                                const std::string &name = "GeometricShapeBall");
    virtual ~GeometricShapeBall() {};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;

  protected:
    virtual BoundingBoxd findBounds() override;
};
} // namespace SPH

#endif // GEOMETRIC_SHAPE_H
