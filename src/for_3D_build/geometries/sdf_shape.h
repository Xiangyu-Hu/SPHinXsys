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
 * @file 	  sdf_shape.h
 * @brief   Here, we define the 3D geometries based on the analytical signed distance function.
 * @details The idea is to define complex geometry by using sdf primitives and binary operations.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SDF_SHAPE_H
#define SDF_SHAPE_H

#include "base_geometry.h"
#include "sdf_primitive.h"
#include "sdf_extension.hpp"
#include "sdf_operation.hpp"

namespace SPH
{

class SDFBase : public Entity
{
  public:
    explicit SDFBase(const std::string &name) : Entity(name) {};
    virtual ~SDFBase() {};
    virtual Real operator()(const Vec3d &point) const = 0;
    virtual BoundingBoxd findBounds() const = 0;
};

template <typename SDFPrimitive>
class SDFEntity : public SDFBase
{
    SDFPrimitive sdf_primitive_;

  public:
    explicit SDFEntity(const std::string &name, const SDFPrimitive &sdf_primitive)
        : SDFBase(name), sdf_primitive_(sdf_primitive) {};
    virtual ~SDFEntity() {};
    SDFPrimitive &getSDFPrimitive() { return sdf_primitive_; };
    virtual Real operator()(const Vec3d &point) const override { return sdf_primitive_(point); };
    virtual BoundingBoxd findBounds() const override { return sdf_primitive_.findBounds(); };
};

using SDFPrimitiveAndOp = std::pair<SDFBase *, GeometricOps>;

class SDFShape : public Shape
{
    UniquePtrsKeeper<SDFBase> sdf_ptrs_;
    Real finest_grid_spacing_;

  public:
    explicit SDFShape(Real finest_grid_spacing, const std::string &shape_name);
    virtual ~SDFShape() {};
    EntityManager &getSDFManager() { return sdf_manager_; };

    template <typename SDFPrimitive>
    SDFShape &insertSDFPrimitive(
        const std::string &primitive_name, const SDFPrimitive &sdf_primitive, const GeometricOps &op)
    {
        SDFEntity<SDFPrimitive> *sdf_entity = sdf_ptrs_.createPtr<
            SDFEntity<SDFPrimitive>>(primitive_name, sdf_primitive);
        sdf_manager_.addEntity<SDFEntity<SDFPrimitive>>(sdf_entity);
        primitives_and_ops_.push_back(SDFPrimitiveAndOp(sdf_entity, op));
        return *this;
    };
    /** Only reliable when the probe point is close to the shape surface.
     * Need to be combined with level set shape and sign correction to avoid artifacts
     * when probe distance is far from the surface. */
    virtual bool checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vec3d findClosestPoint(const Vec3d &probe_point) override;
    virtual BoundingBoxd findBounds() override;

  protected:
    EntityManager sdf_manager_;
    StdVec<SDFPrimitiveAndOp> primitives_and_ops_;

    Real probeSignedDistance(const Vec3d &probe_point);
    Vecd probeNormalDirection(const Vec3d &probe_point);
};
} // namespace SPH
#endif // SDF_SHAPE_H
