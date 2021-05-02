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

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"

#include <string>

namespace SPH
{
	class Tree;

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
	 * @brief template base class of linear structure
	 * only with topology information
	 */
	template<typename InEdgeType, typename OutEdgeType>
	class Edge
	{
	public:
		Edge() : id_(0) {}; /**< constructor without specifying a leading-in edge */
		Edge(InEdgeType in_edge) : Edge() /**< constructor with specifying a leading-in edge */
		{
			in_edge_ = in_edge;
		};
		virtual ~Edge() {};

		size_t id_;					/**< id of this edge */
		InEdgeType in_edge_;		/**< id(s) of parent edge(s) */
		OutEdgeType out_edge_;		/**< id(s) of child edge(s) */
	};

	/**
	 * @class Branch
	 * @brief Each branch has a parent and several children, 
	 * and geometric information. It is a decorated and realized edge. 
	 * Many connected branches compose a tree.
	 */
	class Branch : public Edge<size_t, IndexVector>
	{
	public:
		/** construct the first or leading-in branch of a tree */
		Branch(Vecd init_point, Vecd auxillary_point, Tree* tree);
		/** construct a branch and connect with its parent */
		Branch(size_t parent_id, Tree* tree);
		virtual ~Branch() {};

		Vecd end_direction_;	/**< direction of the last segment of the branch.*/
		/** Inner point indexes of this branch. 
		 * The first is the last inner point from the parent or init_point,
		 * and the last is the first inner point of all its child branches.
		 */
		IndexVector inner_points_;	
		bool is_end_;	/**< whether is an end branch or not */
	};

	/**
	 * @class Structure
	 * @brief Base class for all structures
	 */
	class Structure
	{
	public:
		Structure() {};
		virtual ~Structure() {};

		StdVec<Vecd> points_;				/**< list of the global points containing the coordinates */
		IndexVector edge_locations_;		/**< in which edges are the points located */
		IndexVector face_locations_;		/**< in which faces are the points located */

		size_t EdgeLocation(size_t point_index);
	};

	/**
	 * @class Tree
	 * @brief tree structure
	 */
	class Tree : public Structure
	{
	public:
		Tree(Vecd init_point, Vecd auxillary_point);
		virtual ~Tree();

		size_t last_branch_id_;
		StdVec<Branch*> branches_;	/**< list of all branches */
		void addANewBranch(Branch* branch);
		void addANewBranchInnerVecd(Branch* branch, Vecd new_point, Vecd end_direction);
	};
}
