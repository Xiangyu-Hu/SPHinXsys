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
 * @file 	  triangle_mesh_shape.h
 * @brief   Here, we define the 3D geometries based on the poly mesh.
 * @details The idea is to define complex geometry by passing stl, obj or other poly mesh files.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef TRIANGULAR_MESH_SHAPE_H
#define TRIANGULAR_MESH_SHAPE_H

#include "TriangleMeshDistance.h"
#include "base_geometry.h"
#include "simtk_wrapper.h"
#include "stl_reader.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

namespace SPH
{
class SPHSystem;
/**
 * @class TriangleMeshShape
 * @brief Derived class for triangle shape processing.
 */
class TriangleMeshShape : public Shape
{
  public:
    explicit TriangleMeshShape(const std::string &shape_name);
    /** Only reliable when the probe point is close to the shape surface.
     * Need to be combined with level set shape and sign correction to avoid artifacts
     * when probe distance is far from the surface. */
    virtual bool checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vec3d findClosestPoint(const Vec3d &probe_point) override;
    virtual BoundingBoxd findBounds() override;
    StdVec<std::array<Real, 3>> &getVertices() { return vertices_; }
    StdVec<std::array<int, 3>> &getFaces() { return faces_; }
    void writeMeshToFile(SPHSystem &sph_system, Transform transform = Transform());

  protected:
    StdVec<std::array<Real, 3>> vertices_;
    StdVec<std::array<int, 3>> faces_;
    tmd::TriangleMeshDistance triangle_mesh_distance_;
    void initializeFromPolygonalMesh(const SimTK::PolygonalMesh &poly_mesh);
    void initializeFromSTLMesh(const std::string &file_path_name, Vec3d translation, Real scale_factor);
};

/**
 * @class TriangleMeshShapeBrick
 * @brief Generate a brick triangle mesh using SIMBODy default shape.
 */
class TriangleMeshShapeBrick : public TriangleMeshShape
{
  public:
    class ShapeParameters
    {
      public:
        ShapeParameters() : halfsize_(Vec3d::Zero()), translation_(Vec3d::Zero()), resolution_(0) {};
        Vec3d halfsize_;
        Vec3d translation_;
        int resolution_;
    };
    explicit TriangleMeshShapeBrick(Vec3d halfsize, int resolution, Vec3d translation,
                                    const std::string &shape_name = "TriangleMeshShapeBrick");
    explicit TriangleMeshShapeBrick(const TriangleMeshShapeBrick::ShapeParameters &shape_parameters,
                                    const std::string &shape_name = "TriangleMeshShapeBrick");
    virtual ~TriangleMeshShapeBrick() {};
};

/**
 * @class TriangleMeshShapeSphere
 * @brief Generate a sphere triangle mesh using SIMBODy default shape.
 */
class TriangleMeshShapeSphere : public TriangleMeshShape
{
  public:
    explicit TriangleMeshShapeSphere(Real radius, int resolution, Vec3d translation,
                                     const std::string &shape_name = "TriangleMeshShapeSphere");
    virtual ~TriangleMeshShapeSphere() {};
};

/**
 * @class TriangleMeshShapeCylinder
 * @brief Generate a cylinder triangle mesh using SIMBODy default shape.
 */
class TriangleMeshShapeCylinder : public TriangleMeshShape
{
  public:
    explicit TriangleMeshShapeCylinder(Vec3d axis, Real radius,
                                       Real halflength, int resolution, Vec3d translation,
                                       const std::string &shape_name = "TriangleMeshShapeCylinder");
    virtual ~TriangleMeshShapeCylinder() {};
};

/**
 * @class TriangleMeshShapeSTL
 * @brief Input triangle mesh with stl file.
 */
class TriangleMeshShapeSTL : public TriangleMeshShape
{
  public:
    explicit TriangleMeshShapeSTL(const std::string &file_path_name, Vec3d translation, Real scale_factor,
                                  const std::string &shape_name = "TriangleMeshShapeSTL");
    virtual ~TriangleMeshShapeSTL() {};
};
} // namespace SPH

#endif // TRIANGULAR_MESH_SHAPE_H
