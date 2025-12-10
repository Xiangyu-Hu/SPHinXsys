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
 * @file base_geometry.h
 * @brief Define the base classes Shape, BinaryShape and Edge,
 * which are the base classes for all geometries.
 * @author	Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */

#ifndef BASE_GEOMETRY_H
#define BASE_GEOMETRY_H

#include "base_data_type_package.h"
#include "sphinxsys_containers.h"

#include <spdlog/spdlog.h>
#include <string>

namespace SPH
{
/**
 * @class 	ShapeBooleanOps
 * @brief 	Boolean operation for generate complex shapes
 * @details Note that, for 2d multi polygons, all four operations are implemented.
 * But for binary shapes and complex shapes,
 * only add and sub boolean operation have been defined for right now.
 * Also after operations all surfaces of all shapes should be still surfaces.
 */
enum class ShapeBooleanOps
{
    add,
    sub,
    sym_diff,
    intersect
};

/**
 * @class Shape
 * @brief Base class for all volumetric geometries
 * @details Several pure virtual functions are defined here.
 * (a) closest point on surface: to find the closest point on shape surface to a given point.
 * (b) check if a point contained by the shape,
 * (c) and find the lower and upper bounds.
 */
class Shape
{
  public:
    BoundingBoxd bounding_box_;

    explicit Shape(const std::string &shape_name);
    virtual ~Shape() {};

    std::string getName() { return name_; };
    void setName(const std::string &name) { name_ = name; };
    BoundingBoxd getBounds();
    virtual bool isValid() { return true; };
    virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) = 0;
    virtual Vecd findClosestPoint(const Vecd &probe_point) = 0;

    bool checkNotFar(const Vecd &probe_point, Real threshold);
    bool checkNearSurface(const Vecd &probe_point, Real threshold);
    /** Signed distance is negative for point within the shape. */
    Real findSignedDistance(const Vecd &probe_point);
    /** Normal direction point toward outside of the shape. */
    Vecd findNormalDirection(const Vecd &probe_point);
    virtual BoundingBoxd findBounds() = 0;

  protected:
    std::string name_;
    bool is_bounds_found_;
    std::shared_ptr<spdlog::logger> logger_;
};

using SubShapeAndOp = std::pair<Shape *, ShapeBooleanOps>;
/**
 * @class BinaryShapes
 * @brief A collections of shapes with binary operations.
 * This class has ownership of all shapes by using a unique pointer vector.
 * In this way, add or subtract a shape will call the shape's constructor other than
 * passing the shape pointer.
 */
class BinaryShapes : public Shape
{
  public:
    BinaryShapes() : Shape("BinaryShapes") {};
    explicit BinaryShapes(const std::string &shape_name) : Shape(shape_name) {};
    virtual ~BinaryShapes() {};

    void add(Shape *sub_shape)
    {
        SubShapeAndOp sub_shape_and_op(sub_shape, ShapeBooleanOps::add);
        sub_shapes_and_ops_.push_back(sub_shape_and_op);
    };

    template <class SubShapeType, typename... Args>
    void add(Args &&...args)
    {
        Shape *sub_shape = sub_shape_ptrs_keeper_.createPtr<SubShapeType>(std::forward<Args>(args)...);
        add(sub_shape);
    };

    void subtract(Shape *sub_shape)
    {
        SubShapeAndOp sub_shape_and_op(sub_shape, ShapeBooleanOps::sub);
        sub_shapes_and_ops_.push_back(sub_shape_and_op);
    };

    template <class SubShapeType, typename... Args>
    void subtract(Args &&...args)
    {
        Shape *sub_shape = sub_shape_ptrs_keeper_.createPtr<SubShapeType>(std::forward<Args>(args)...);
        subtract(sub_shape);
    };

    virtual bool isValid() override;
    virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;
    virtual BoundingBoxd findBounds() override;
    Shape *getSubShapeByName(const std::string &name);
    SubShapeAndOp *getSubShapeAndOpByName(const std::string &name);
    size_t getSubShapeIndexByName(const std::string &name);

  protected:
    UniquePtrsKeeper<Shape> sub_shape_ptrs_keeper_;
    StdVec<SubShapeAndOp> sub_shapes_and_ops_;
};

/**
 * @class Edge
 * @brief template base class of linear structure only with topology information.
 * Note that a edge is defined together with a structure which is composed of edges.
 * Such structure should have an interface function ContainerSize() returning
 * the curent total amount of edges.
 */
template <typename InEdgeType, typename OutEdgeType>
class Edge
{
  public:
    /** constructor without specifying a leading-in edge */
    template <class EdgeStructureType>
    explicit Edge(EdgeStructureType *structure)
        : id_(structure->ContainerSize()), in_edge_(MaxSize_t){};
    /** constructor with specifying a leading-in edge */
    template <class EdgeStructureType>
    Edge(InEdgeType in_edge, EdgeStructureType *structure)
        : id_(structure->ContainerSize()), in_edge_(in_edge){};
    virtual ~Edge() {};

    size_t id_;            /**< id of this edge */
    InEdgeType in_edge_;   /**< id(s) of parent edge(s) */
    OutEdgeType out_edge_; /**< id(s) of child edge(s) */
};
} // namespace SPH
#endif // BASE_GEOMETRY_H