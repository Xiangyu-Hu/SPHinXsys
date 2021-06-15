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



#include "base_data_package.h"
#include "sph_data_conainers.h"

#include <string>

namespace SPH
{
	class Tree;
	class Neighborhood;
	class BaseParticles;
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
		/** generalized particle search algorithm */
		template<typename GetParticleIndex, typename GetSearchRange, typename GetNeighborRelation>
		void searchNeighborsByParticles(size_t number_of_particles, BaseParticles& source_particles, 
			StdLargeVec<Neighborhood>& particle_configuration, GetParticleIndex& get_particle_index,
			GetSearchRange& get_search_range, GetNeighborRelation& get_neighbor_relation) 
		{
			size_t particle_id;
			size_t parent_branch_id;
			size_t child_branch_id;
			size_t num_ele;
			std::vector<size_t> neighboring_ids;
			std::vector<size_t> child_ids;
			/** First branch
			* Note that the first branc has only one point, accordingly, one particle generated.
			* Find the neibors in child branch, the first branch only have one child, id = 1.
			*/
			particle_id = branches_[0]->inner_points_.front();
			neighboring_ids.clear();
			neighboring_ids.push_back(branches_[1]->inner_points_[0]);
			neighboring_ids.push_back(branches_[1]->inner_points_[1]);
			/** Build configuration. */
			Neighborhood& neighborhood = particle_configuration[particle_id];
			for (size_t n = 0; n != neighboring_ids.size(); ++n)
			{
				Vecd displacement = source_particles.pos_n_[particle_id] - source_particles.pos_n_[neighboring_ids[n]];
				get_neighbor_relation(neighborhood, displacement, particle_id, neighboring_ids[n]);
			}
			/** Second branch. 
		 	* The second branch has special parent branch, branch 0, consisting only one point.
		 	* The child branch are two normal branch. 
		 	*/
			num_ele = branches_[1]->inner_points_.size();
			child_ids.clear();
			for(int k = 0; k < branches_[1]->out_edge_.size(); ++k)
			{
				child_ids.push_back(branches_[1]->out_edge_[k]);
			}

			for(int i = 0; i != num_ele; i++)
			{
				neighboring_ids.clear();
				particle_id = branches_[1]->inner_points_.front() + i;
				if(i == 0)
				{
					neighboring_ids.push_back(branches_[0]->inner_points_.front());
					neighboring_ids.push_back(particle_id + 1);
					neighboring_ids.push_back(particle_id + 2);
				} else if(i == 1){
					neighboring_ids.push_back(branches_[0]->inner_points_.front());
					neighboring_ids.push_back(particle_id - 1);
					neighboring_ids.push_back(particle_id + 1);
					neighboring_ids.push_back(particle_id + 2);
				}else if(2 <= i &&  i <= (num_ele - 3))
				{
					neighboring_ids.push_back(particle_id - 1);
					neighboring_ids.push_back(particle_id - 2);
					neighboring_ids.push_back(particle_id + 1);
					neighboring_ids.push_back(particle_id + 2);
				} else if(i == (num_ele - 2)) 
				{
					neighboring_ids.push_back(particle_id - 2);
					neighboring_ids.push_back(particle_id - 1);
					neighboring_ids.push_back(particle_id + 1);

					for(size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
					{
						child_branch_id = branches_[1]->out_edge_[k];
						neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front());
					}
				}else if(i == (num_ele - 1))
				{
					neighboring_ids.push_back(particle_id - 1);
					neighboring_ids.push_back(particle_id - 2);

					for(int k = 0; k < branches_[1]->out_edge_.size(); ++k)
					{
						child_branch_id = branches_[1]->out_edge_[k];
						neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front());
						neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front() + 1);
					}
				}

				Neighborhood& neighborhood = particle_configuration[particle_id];
				for (int n = 0; n != neighboring_ids.size(); ++n)
				{
					Vecd displacement = source_particles.pos_n_[particle_id] - source_particles.pos_n_[neighboring_ids[n]];
					get_neighbor_relation(neighborhood, displacement, particle_id, neighboring_ids[n]);
				}
			}
			/** Other branches. 
		 	* They are may normal branch (fully growed, has child and parent) or non-fully growed branch (no child or less segmetns)
		 	*/
			for (size_t branch_idx = 2; branch_idx != branches_.size(); ++branch_idx)
			{
				num_ele = branches_[branch_idx]->inner_points_.size();
				parent_branch_id = branches_[branch_idx]->in_edge_;
				if(!branches_[branch_idx]->is_end_)
				{
					/** This branch is fully growed. */
					for(size_t i = 0; i != num_ele; i++)
					{
						neighboring_ids.clear();
						particle_id = branches_[branch_idx]->inner_points_.front() + i;
						if(i == 0)
						{
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back());
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back() - 1);
								
							neighboring_ids.push_back(particle_id + 1);
							neighboring_ids.push_back(particle_id + 2);
						} else if(i == 1){
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back());
							neighboring_ids.push_back(particle_id - 1);
							neighboring_ids.push_back(particle_id + 1);
							neighboring_ids.push_back(particle_id + 2);
						}else if(2 <= i &&  i <= (num_ele - 3))
						{
							neighboring_ids.push_back(particle_id - 1);
							neighboring_ids.push_back(particle_id - 2);
							neighboring_ids.push_back(particle_id + 1);
							neighboring_ids.push_back(particle_id + 2);
						} else if(i == (num_ele - 2)) 
						{
							neighboring_ids.push_back(particle_id - 2);
							neighboring_ids.push_back(particle_id - 1);
							neighboring_ids.push_back(particle_id + 1);

							for(int k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
							{
								child_branch_id = branches_[branch_idx]->out_edge_[k];
								neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front());
							}
						}else if(i == (num_ele - 1))
						{
							neighboring_ids.push_back(particle_id - 1);
							neighboring_ids.push_back(particle_id - 2);

							for(int k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
							{
								child_branch_id = branches_[branch_idx]->out_edge_[k];
								neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front());
								if(branches_[child_branch_id]->inner_points_.size() >= 2)
								{
									neighboring_ids.push_back(branches_[child_branch_id]->inner_points_.front() + 1);

								}
							}
						}

						Neighborhood& neighborhood = particle_configuration[particle_id];
						for (size_t n = 0; n != neighboring_ids.size(); ++n)
						{
							Vecd displacement = source_particles.pos_n_[particle_id] - source_particles.pos_n_[neighboring_ids[n]];
							get_neighbor_relation(neighborhood, displacement, particle_id, neighboring_ids[n]);
						}
					}
				}else
				{
					/** This branch is not fully growed. */
					for(size_t i = 0; i != num_ele; i++)
					{
						neighboring_ids.clear();
						particle_id = branches_[branch_idx]->inner_points_.front() + i;
						if(i == 0){
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back());
							if(branches_[parent_branch_id]->inner_points_.size() >= 2)
								neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back() - 1);
						} else if(i == 1){
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_points_.back());
							neighboring_ids.push_back(particle_id - 1);
						}else{
							neighboring_ids.push_back(particle_id - 1);
							neighboring_ids.push_back(particle_id - 2);
						}

						if(i + 1 < num_ele)
							neighboring_ids.push_back(particle_id + 1);
						if(i + 2 < num_ele)
							neighboring_ids.push_back(particle_id + 2);
							
						Neighborhood& neighborhood = particle_configuration[particle_id];
						for (size_t n = 0; n != neighboring_ids.size(); ++n)
						{
							Vecd displacement = source_particles.pos_n_[particle_id] - source_particles.pos_n_[neighboring_ids[n]];
							get_neighbor_relation(neighborhood, displacement, particle_id, neighboring_ids[n]);
						}
					}
				}
			}
			/** Delete the data in tree. */
			//~Tree();
		};
	};
}
#endif //BASE_GEOMETRY_H