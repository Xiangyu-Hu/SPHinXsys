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
	enum class ShapeBooleanOps { add, sub, sym_diff, intersect };

	/**
	 * @class Shape
	 * @brief Base class for all geometries
	 */
	class Shape
	{
	public:
		Shape(std::string shape_name) : name_(shape_name) {};
		virtual ~Shape() {};

		std::string getName() { return name_; };
		virtual BoundingBox findBounds() = 0;
	protected:
		std::string name_;
	};

	/**
	 * @class Edge
	 * @brief template base class of linear structure only with topology information.
	 * Note that a edge is defined together with a structure which is composed of edges.
	 * Such structure should have an interface function ContainerSize() returning 
	 * the curent total amount of edges.
	 */
	template<typename InEdgeType, typename OutEdgeType>
	class Edge
	{
	public:
		/** constructor without specifying a leading-in edge */
		template<class EdgeStructureType>
		Edge(EdgeStructureType *structure) 
			: id_(structure->ContainerSize()), in_edge_(MaxSize_t) {}; 
		 /** constructor with specifying a leading-in edge */	
		template<class EdgeStructureType>
		Edge(InEdgeType in_edge, EdgeStructureType *structure)
			: id_(structure->ContainerSize()), in_edge_(in_edge) {};
		virtual ~Edge() {};

		size_t id_;					/**< id of this edge */
		InEdgeType in_edge_;		/**< id(s) of parent edge(s) */
		OutEdgeType out_edge_;		/**< id(s) of child edge(s) */
	};
}
#endif //BASE_GEOMETRY_H