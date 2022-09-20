/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file base_geometry.h
 * @brief Shape is the base class for all geometries.
 * @details Several pure virtual functions
 * are defined here. (a) closest point on surface: to find the closest point on shape
 * surface to a given point. (b) find the lower and upper bounds.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef BASE_GEOMETRY_H
#define BASE_GEOMETRY_H

#include "base_data_package.h"
#include "sph_data_containers.h"

#include <string>

namespace SPH
{
	/**
	 * @class ShapeBooleanOps
	 * @brief Boolean operation for generate complex shapes
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
	 * Note that checkContain and findClosest point are  basic function,
	 * They should not call other functions in shape.
	 */
	class Shape
	{
	public:
		explicit Shape(const std::string &shape_name)
			: name_(shape_name), is_bounds_found_(false){};
		virtual ~Shape(){};

		std::string getName() { return name_; };
		void setName(const std::string &name) { name_ = name; };
		BoundingBox getBounds();
		virtual bool isValid() { return true; };
		virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) = 0;
		virtual Vecd findClosestPoint(const Vecd &input_pnt) = 0;

		bool checkNotFar(const Vecd &input_pnt, Real threshold);
		bool checkNearSurface(const Vecd &input_pnt, Real threshold);
		/** Signed distance is negative for point within the complex shape. */
		Real findSignedDistance(const Vecd &input_pnt);
		/** Normal direction point toward outside of the complex shape. */
		Vecd findNormalDirection(const Vecd &input_pnt);

	protected:
		std::string name_;
		BoundingBox bounding_box_;
		bool is_bounds_found_;

		virtual BoundingBox findBounds() = 0;
	};

	using ShapeAndOp = std::pair<Shape *, ShapeBooleanOps>;
	/**
	 * @class BinaryShapes
	 * @brief a collections of shapes with binary operations
	 * This class so that it has ownership of all shapes by using a unique pointer vector.
	 * In this way, add or subtract a shape will call the shape's constructor other than
	 * passing the shape pointer.
	 */
	class BinaryShapes : public Shape
	{
	public:
		BinaryShapes() : Shape("BinaryShapes"){};
		explicit BinaryShapes(const std::string &shapes_name) : Shape(shapes_name){};
		virtual ~BinaryShapes(){};

		template <class ShapeType, typename... Args>
		void add(Args &&...args)
		{
			Shape *shape = shapes_ptr_keeper_.createPtr<ShapeType>(std::forward<Args>(args)...);
			ShapeAndOp shape_and_op(shape, ShapeBooleanOps::add);
			shapes_and_ops_.push_back(shape_and_op);
		};

		void add(Shape& shape)
		{
			shapes_and_ops_.push_back({&shape, ShapeBooleanOps::add});
		}

		template <class ShapeType, typename... Args>
		void subtract(Args &&...args)
		{
			Shape *shape = shapes_ptr_keeper_.createPtr<ShapeType>(std::forward<Args>(args)...);
			ShapeAndOp shape_and_op(shape, ShapeBooleanOps::sub);
			shapes_and_ops_.push_back(shape_and_op);
		};

		void subtract(Shape& shape)
		{
			shapes_and_ops_.push_back({&shape, ShapeBooleanOps::sub});
		}

		virtual bool isValid() override;
		virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vecd findClosestPoint(const Vecd &input_pnt) override;
		Shape *getShapeByName(const std::string &shape_name);
		ShapeAndOp *getShapeAndOpByName(const std::string &shape_name);
		size_t getShapeIndexByName(const std::string &shape_name);

	protected:
		UniquePtrKeepers<Shape> shapes_ptr_keeper_;
		StdVec<ShapeAndOp> shapes_and_ops_;

		virtual BoundingBox findBounds() override;
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
		virtual ~Edge(){};

		size_t id_;			   /**< id of this edge */
		InEdgeType in_edge_;   /**< id(s) of parent edge(s) */
		OutEdgeType out_edge_; /**< id(s) of child edge(s) */
	};
}
#endif // BASE_GEOMETRY_H