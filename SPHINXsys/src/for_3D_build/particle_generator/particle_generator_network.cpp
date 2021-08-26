/**
 * @file 	particle_generator_network.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */
#include "sph_system.h"
#include "particle_generator_network.h"
#include "cell_linked_list.h"
#include "level_set.h"
#include "base_body.h"
#include "base_particles.h"
#include "in_output.h"
 //=================================================================================================//
namespace SPH 
{
	//=================================================================================================//
	ParticleGeneratorNetwork::
		ParticleGeneratorNetwork(Vecd starting_pnt, Vecd second_pnt, int iterator, Real grad_factor)
		: ParticleGenerator(), starting_pnt_(starting_pnt), second_pnt_(second_pnt),
		n_it_(iterator), fascicles_(true), segments_in_branch_(10), segment_length_(0), 
		grad_factor_(grad_factor), body_shape_(nullptr), tree_(nullptr){}
	//=================================================================================================//
	void ParticleGeneratorNetwork::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		cell_linked_list_ = dynamic_cast<RealBody*>(sph_body)->cell_linked_list_;
		segment_length_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		body_shape_ = sph_body_->body_shape_;
		tree_ = new GenerativeTree(sph_body_);
		sph_body_->generative_structure_ = tree_;
		Vecd displacement = second_pnt_ - starting_pnt_;
		Vecd end_direction = displacement / (displacement.norm() + TinyReal);
		//add particle to the first branch of the tree
		tree_->growAParticleOnBranch(tree_->root_, starting_pnt_, end_direction);
		cell_linked_list_->InsertACellLinkedListDataEntry(0, tree_->pos_n_[0]);
	}
	//=================================================================================================//
	Vecd ParticleGeneratorNetwork::getGradientFromNearestPoints(Vecd pt, Real delta)
	{
		Vecd upgrad(0), downgrad(0);
		Vecd shift(delta);
		for (int i = 0; i != Dimensions; i++) {
			Vecd upwind = pt;
			Vecd downwind = pt;
			upwind[i] -= shift[i];
			downwind[i] += shift[i];
			ListData up_nearest_list = cell_linked_list_->findNearestListDataEntry(upwind);
			ListData down_nearest_list = cell_linked_list_->findNearestListDataEntry(downwind);
			upgrad[i] = up_nearest_list.first != MaxSize_t ? (upwind - up_nearest_list.second).norm() / 2.0 * delta : 1.0;
			downgrad[i] = down_nearest_list.first != MaxSize_t ? (downwind - down_nearest_list.second).norm() / 2.0 * delta : 1.0;
		}
		return downgrad - upgrad;
	}
	//=================================================================================================//
	Vecd ParticleGeneratorNetwork::createATentativeNewBranchPoint(Vecd init_point, Vecd dir)
	{
		Vecd pnt_to_project = init_point + dir * segment_length_;

		Real phi = body_shape_->findSignedDistance(pnt_to_project);
		Vecd unit_normal = body_shape_->findNormalDirection(pnt_to_project);
		unit_normal /= unit_normal.norm() + TinyReal;
		Vecd new_point = pnt_to_project - phi * unit_normal;
		return new_point;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::
		isCollision(Vecd& new_point, ListData& nearest_neighbor, size_t parent_id)
	{
		bool collision = false;
		bool is_family = false;

		collision = extraCheck(new_point);

		size_t edge_location = tree_->BranchLocation(nearest_neighbor.first);
		if (edge_location == parent_id) is_family = true;
		for (const size_t& brother_branch : tree_->branches_[parent_id]->out_edge_)
		{
			if (edge_location == brother_branch) is_family = true;
		}

		if (!is_family)
		{
			Real min_distance = (new_point - nearest_neighbor.second).norm();
			if (min_distance < 5.0 * segment_length_) collision = true;
		}

		return collision;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::
		createABranchIfValid(SPHBody* sph_body, size_t parent_id, Real angle,
							 Real repulsivity, size_t number_segments)
	{
		bool is_valid = false;
		GenerativeTree::Branch* parent_branch = tree_->branches_[parent_id];
		IndexVector& parent_elements = parent_branch->inner_particles_;
		StdLargeVec<Vecd> &tree_points = tree_->pos_n_;
	
		Vecd init_point = tree_points[parent_elements.back()];
		Vecd init_direction = parent_branch->end_direction_;


		Vecd surface_norm = body_shape_->findNormalDirection(init_point);
		surface_norm /= surface_norm.norm() + TinyReal;
		Vecd in_plane = - SimTK::cross(init_direction, surface_norm);

		Real delta = grad_factor_ * segment_length_;
		Vecd grad = getGradientFromNearestPoints(init_point, delta);
		Vecd dir = cos(angle) * init_direction + sin(angle) * in_plane;
		dir /= dir.norm() + TinyReal;
		Vecd end_direction = (repulsivity * grad + dir) / ((repulsivity * grad + dir).norm() + TinyReal);
		Vecd end_point = init_point;

		Vecd new_point = createATentativeNewBranchPoint(end_point, end_direction);
		ListData nearest_neighbor = cell_linked_list_->findNearestListDataEntry(new_point);
		if (!isCollision(new_point, nearest_neighbor, parent_id)) 
		{
			is_valid = true;
			GenerativeTree::Branch* new_branch = new GenerativeTree::Branch(parent_id, tree_);
			tree_->growAParticleOnBranch(new_branch, new_point, end_direction);

			for (size_t i = 1; i < number_segments; i++)
			{
				surface_norm = body_shape_->findNormalDirection(new_point);
				surface_norm /= surface_norm.norm() + TinyReal;
				/** Project grad to surface. */
				grad = getGradientFromNearestPoints(new_point, delta);
				grad -= dot(grad, surface_norm) * surface_norm;
				dir = (repulsivity * grad + end_direction) / ((repulsivity * grad + end_direction).norm() + TinyReal);
				end_direction = dir;
				end_point = new_point;

				new_point = createATentativeNewBranchPoint(end_point, end_direction);
				ListData nearest_neighbor = cell_linked_list_->findNearestListDataEntry(new_point);
				if (isCollision(new_point, nearest_neighbor, parent_id))
				{
					new_branch->is_terminated_ = true;
					std::cout << "Branch Collision Detected, Break! " << std::endl;
					break;
				}
				/** This constraint imposed to avoid too small time step size. */
				if((new_point - end_point).norm() < 0.5 * segment_length_)
				{
					new_branch->is_terminated_ = true;
					std::cout << "New branch point is too close, Break! " << std::endl;
					break;					
				}
				tree_->growAParticleOnBranch(new_branch, new_point, end_direction);

			}

			for (const size_t& particle_idx : new_branch->inner_particles_)
			{
				cell_linked_list_->InsertACellLinkedListDataEntry(particle_idx, tree_points[particle_idx]);
			}
		}

		return is_valid;
	}
	//=================================================================================================//
	void ParticleGeneratorNetwork::createBaseParticles(BaseParticles* base_particles)
	{
		In_Output* in_output = sph_body_->getSPHSystem().in_output_;
		BodyStatesRecordingToVtu 	write_states(*in_output, { sph_body_ });

		std::cout << "Now creating Particles on network... " << "\n" << std::endl;

		//the second branch
		bool is_valid = createABranchIfValid(sph_body_, 0, 0.0, 0.0, segments_in_branch_);

		size_t ite = 0;
		sph_body_->setNewlyUpdated();
		write_states.writeToFile(0);

		IndexVector branches_to_grow;
		IndexVector new_branches_to_grow;
        if (is_valid) branches_to_grow.push_back(tree_->last_branch_id_);
	
		if (fascicles_)
		{
			/** Set vertices in family branch. */
			branches_to_grow.clear();
			for (size_t i = 0; i != 2; i++)
			{
				/** Creating a new branch. */
				Real  angle_to_use = fascicle_angles_[i];
				size_t fascicles_segments = int(fascicle_ratio_ * segments_in_branch_);
				bool is_valid = createABranchIfValid(sph_body_, 1, angle_to_use, 0.0, fascicles_segments);
				if (is_valid) branches_to_grow.push_back(tree_->last_branch_id_);
			}

			ite++;
			sph_body_->setNewlyUpdated();
			write_states.writeToFile(ite);

		}

		for(size_t i = 0; i != n_it_; i++)
        {
			new_branches_to_grow.clear();
			random_shuffle(branches_to_grow.begin(), branches_to_grow.end());
            for(size_t j = 0; j != branches_to_grow.size(); j++)
            {
                size_t grow_id = branches_to_grow[j];
				Real rand_num = ((Real)rand() / (RAND_MAX)) - 0.5;
				Real angle_to_use = angle_ + rand_num * 0.05;
                for(size_t k = 0; k != 2; k++)
                {   
                    /** Creating a new branch with fixed number of segments. */
                    size_t random_number_segments = segments_in_branch_;
					bool is_valid = createABranchIfValid(sph_body_, grow_id, angle_to_use, repulsivity_,
						random_number_segments);

					if(is_valid && !tree_->LastBranch()->is_terminated_)
                    {
                        new_branches_to_grow.push_back(tree_->last_branch_id_);
                    }

                    angle_to_use *= -1.0;
                }
            }
            branches_to_grow = new_branches_to_grow;

			ite++;
			sph_body_->setNewlyUpdated();
			write_states.writeToFile(ite);
		}
        
		std::cout << base_particles->total_real_particles_ << " Particles has been successfully created!" << "\n" << std::endl;
	}
	//=================================================================================================//
}
//=================================================================================================//