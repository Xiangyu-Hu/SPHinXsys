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
* are defined here. (a) closet point on surface: to find the closet point on shape
* surface to a given point. (b) find the lower and upper bounds.
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef BASE_GEOMETRY_H
#define BASE_GEOMETRY_H

#include "base_particles.h"
#include "base_data_package.h"
#include "sph_data_containers.h"

#include <string>

namespace SPH
{
	class Tree;
	class Neighborhood;
	/**
	 * @class ShapeBooleanOps
	 * @brief Boolian operation for generate complex shapes
	 * @details Note that, for 3D applications, only add and sub boolean operation have been defined right now
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
	 */
	class Shape
	{
	public:
		explicit Shape(const std::string &shape_name) : name_(shape_name){};
		virtual ~Shape(){};

		std::string getName() { return name_; };
		void setName(const std::string &name) { name_ = name; };
		virtual BoundingBox findBounds() = 0;
		virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) = 0;
		virtual Vecd findClosestPoint(const Vecd &input_pnt) = 0;

		virtual bool checkNotFar(const Vecd &input_pnt, Real threshold);
		virtual bool checkNearSurface(const Vecd &input_pnt, Real threshold);
		/** Signed distance is negative for point within the complex shape. */
		virtual Real findSignedDistance(const Vecd &input_pnt);
		/** Normal direction point toward outside of the complex shape. */
		virtual Vecd findNormalDirection(const Vecd &input_pnt);

	protected:
		std::string name_;
	};
	using ShapeAndOp = std::pair<Shape *, ShapeBooleanOps>;

	/**
	 * @class BinaryShapes
	 * @brief a collections of shapes with binary operations
	 * //TODO: I will reformulate this class so that it has ownship of all shapes 
	 * by using a unique pointer vector.
	 * In this way, add or substract a shape will call the shape's constructor other than 
	 * passing the shape pointer.
	 */
	class BinaryShapes : public Shape
	{
	private:
		UniquePtrVectorKeeper<Shape> shapes_ptr_keeper_;

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

		template <class ShapeType, typename... Args>
		void substract(Args &&...args)
		{
			Shape *shape = shapes_ptr_keeper_.createPtr<ShapeType>(std::forward<Args>(args)...);
			ShapeAndOp shape_and_op(shape, ShapeBooleanOps::sub);
			shapes_and_ops_.push_back(shape_and_op);
		};

		virtual BoundingBox findBounds() override;
		virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vecd findClosestPoint(const Vecd &input_pnt) override;
		Shape *getShapeByName(const std::string &shape_name);
		ShapeAndOp *getShapeAndOpByName(const std::string &shape_name);

	protected:
		StdVec<ShapeAndOp> shapes_and_ops_;
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
#endif //BASE_GEOMETRY_H