/**
 * @file 	particle_generator_network.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */
#include "sph_system.h"
#include "particle_generator_network.h"
#include "mesh_cell_linked_list.h"
#include "level_set.h"
#include "base_body.h"
#include "base_particles.h"
#include "in_output.h"
 //=================================================================================================//
namespace SPH 
{
	//=================================================================================================//
	ParticleGeneratorNetwork::ParticleGeneratorNetwork(Vecd starting_pnt, Vecd second_pnt, int iterator, Real grad_factor)
		: ParticleGenerator(), starting_pnt_(starting_pnt), second_pnt_(second_pnt),
		n_it_(iterator), fascicles_(true), segments_in_branch_(10), segment_length_(0), grad_factor_(grad_factor),
		body_shape_(NULL){}
	//=================================================================================================//
	void ParticleGeneratorNetwork::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		real_body_ = dynamic_cast<RealBody*>(sph_body);
		segment_length_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		body_shape_ = sph_body_->body_shape_;
		sph_body_->tree_ = new Tree(starting_pnt_, second_pnt_);
	}
	//=================================================================================================//
	Vecd ParticleGeneratorNetwork::
		getGradientFromNearestPoints(Vecd pt, Real delta, BaseMeshCellLinkedList* mesh_cell_linked_list)
	{
		Vecd upgrad(0), downgrad(0);
		Vecd shift(delta);
		for (int i = 0; i < pt.size(); i++) {
			Vecd upwind = pt;
			Vecd downwind = pt;
			upwind[i] -= shift[i];
			downwind[i] += shift[i];
			ListData up_nearest_list = mesh_cell_linked_list->findNearestListDataEntry(upwind);
			ListData down_nearest_list = mesh_cell_linked_list->findNearestListDataEntry(downwind);
			upgrad[i] = up_nearest_list.first != MaxSize_t ? (upwind - up_nearest_list.second).norm() / 2.0 * delta : 1.0;
			downgrad[i] = down_nearest_list.first != MaxSize_t ? (downwind - down_nearest_list.second).norm() / 2.0 * delta : 1.0;
		}
		return downgrad - upgrad;
	}
	//=================================================================================================//
	Vecd ParticleGeneratorNetwork::creatATentativeNewBranchVecd(Vecd init_point, Vecd dir)
	{
		Vecd pnt_to_project = init_point + dir * segment_length_;

		Real phi = body_shape_->findSignedDistance(pnt_to_project);
		Vecd unit_normal = body_shape_->findNormalDirection(pnt_to_project);
		unit_normal /= unit_normal.norm() + TinyReal;
		Vecd new_point = pnt_to_project - phi * unit_normal;
		return new_point;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::isCollision(Vecd& new_point, 
		ListData& nearest_neighbor, size_t parent_id, Tree* tree)
	{
		bool collision = false;
		bool is_family = false;

		collision = extraCheck(new_point);

		size_t edge_location = tree->EdgeLocation(nearest_neighbor.first);
		if (edge_location == parent_id) is_family = true;
		IndexVector& brother_branches = tree->branches_[parent_id]->out_edge_;
		for (size_t i = 0; i < brother_branches.size(); i++)
		{
			if (edge_location == brother_branches[i]) is_family = true;
		}

		if (!is_family)
		{
			Real min_distance = (new_point - nearest_neighbor.second).norm();
			if (min_distance < 5.0 * segment_length_) collision = true;
		}

		return collision;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::createABranchIfValid(SPHBody* sph_body, size_t parent_id, Real angle,
		Real repulsivity, size_t number_segments, Tree* tree)
	{
		bool is_valid = false;
		Branch* parent_branch = tree->branches_[parent_id];
		IndexVector& parent_elements = parent_branch->inner_points_;
		StdVec<Vecd>& points = tree->points_;
	
		Vecd init_point = points[parent_elements.back()];
		Vecd init_direction = parent_branch->end_direction_;

		BaseMeshCellLinkedList* mesh_cell_linked_list = real_body_->mesh_cell_linked_list_;

		Real delta = grad_factor_ * segment_length_;

		Vecd surface_norm = body_shape_->findNormalDirection(init_point);
		surface_norm /= surface_norm.norm() + TinyReal;
		Vecd in_plane = - SimTK::cross(init_direction, surface_norm);

		Vecd dir = cos(angle) * init_direction + sin(angle) * in_plane;
		dir /= dir.norm() + TinyReal;
		Vecd grad = getGradientFromNearestPoints(init_point, delta, mesh_cell_linked_list);
		Vecd end_direction = (repulsivity * grad + dir) / ((repulsivity * grad + dir).norm() + TinyReal);
		Vecd end_point = init_point;

		Vecd new_point = creatATentativeNewBranchVecd(end_point, end_direction);
		ListData nearest_neighbor = mesh_cell_linked_list->findNearestListDataEntry(new_point);
		if (!isCollision(new_point, nearest_neighbor, parent_id, tree)) {
			is_valid = true;
			Branch* new_branch = new Branch(parent_id, tree);
			tree->addANewBranch(new_branch);
			tree->addANewBranchInnerVecd(new_branch, new_point, end_direction);

			for (size_t i = 1; i < number_segments; i++)
			{
				surface_norm = body_shape_->findNormalDirection(new_point);
				surface_norm /= surface_norm.norm() + TinyReal;
				/** Project grad to surface. */
				grad = getGradientFromNearestPoints(new_point, delta, mesh_cell_linked_list);
				grad -= dot(grad, surface_norm) * surface_norm;
				dir = (repulsivity * grad + end_direction) / ((repulsivity * grad + end_direction).norm() + TinyReal);
				end_direction = dir;
				end_point = new_point;

				new_point = creatATentativeNewBranchVecd(end_point, end_direction);
				ListData nearest_neighbor = mesh_cell_linked_list->findNearestListDataEntry(new_point);
				if (isCollision(new_point, nearest_neighbor, parent_id, tree))
				{
					new_branch->is_end_ = true;
					std::cout << "Branch Collision Detected, Break! " << std::endl;
					break;
				}
				/** This constraint imposed to avoid too small time step size. */
				if((new_point - end_point).norm() < 0.5 * segment_length_)
				{
					new_branch->is_end_ = true;
					break;					
				}
				tree->addANewBranchInnerVecd(new_branch, new_point, end_direction);
			}

			for (size_t i = 0; i != new_branch->inner_points_.size(); i++)
			{
				mesh_cell_linked_list->InsertACellLinkedListDataEntry(new_branch->inner_points_[i], points[new_branch->inner_points_[i]]);
				sph_body->base_particles_->initializeABaseParticle(points[new_branch->inner_points_[i]], segment_length_);
			}
		}

		return is_valid;
	}
	//=================================================================================================//
	void ParticleGeneratorNetwork::createBaseParticles(BaseParticles* base_particles)
	{
		In_Output* in_output = sph_body_->getSPHSystem().in_output_;
		WriteBodyStatesToVtu 	write_states(*in_output, { sph_body_ });

		StdVec<Branch*>& branches = sph_body_->tree_->branches_;
		StdVec<Vecd>& points = sph_body_->tree_->points_;

		base_particles->initializeABaseParticle(starting_pnt_, segment_length_);
		BaseMeshCellLinkedList* mesh_cell_linked_list = real_body_->mesh_cell_linked_list_;
		mesh_cell_linked_list->InsertACellLinkedListDataEntry(0, points[0]);
		std::cout << "Now creating Particles on network... " << "\n" << std::endl;

		//the first branch
		bool is_valid = createABranchIfValid(sph_body_, 0, 0.0, 0.0, segments_in_branch_, sph_body_->tree_);

		size_t ite = 0;
		sph_body_->setNewlyUpdated();
		base_particles->total_real_particles_ = base_particles->pos_n_.size();
		write_states.WriteToFile(0.0);

		IndexVector edges_to_grow;
		IndexVector new_edges_to_grow;
        if (is_valid) edges_to_grow.push_back(sph_body_->tree_->last_branch_id_);
	
		if (fascicles_)
		{
			/** Set vertices in family edge. */
			edges_to_grow.clear();
			for (int i = 0; i < 2; i++)
			{
				/** Creating a new edge. */
				Real  angle_to_use = fascicle_angles_[i];
				size_t fascicles_segments = int(fascicle_ratio_ * segments_in_branch_);
				bool is_valid = createABranchIfValid(sph_body_, 1, angle_to_use, 0.0, fascicles_segments, sph_body_->tree_);
				if (is_valid) edges_to_grow.push_back(sph_body_->tree_->last_branch_id_);
			}

			ite++;
			sph_body_->setNewlyUpdated();
			base_particles->total_real_particles_ = base_particles->pos_n_.size();
			write_states.WriteToFile(Real(ite));

		}

		for(size_t i = 0; i != n_it_; i++)
        {
			new_edges_to_grow.clear();
			random_shuffle(edges_to_grow.begin(), edges_to_grow.end());
            for(size_t j = 0; j != edges_to_grow.size(); j++)
            {
                size_t grow_id = edges_to_grow[j];
				Real rand_num = ((double)rand() / (RAND_MAX)) - 0.5;
				Real angle_to_use = angle_ + rand_num * 0.05;
                for(size_t k = 0; k != 2; k++)
                {   
                    /** Creating a new edge. */
                    size_t random_number_segments = segments_in_branch_;// + rand() % 10 + 1;
					bool is_valid = createABranchIfValid(sph_body_, grow_id, angle_to_use, repulsivity_,
						random_number_segments, sph_body_->tree_);

					if(is_valid && !branches[sph_body_->tree_->last_branch_id_]->is_end_)
                    {
                        new_edges_to_grow.push_back(sph_body_->tree_->last_branch_id_);
                    }

                    angle_to_use *= -1.0;
                }
            }
            edges_to_grow = new_edges_to_grow;

			ite++;
			sph_body_->setNewlyUpdated();
			base_particles->total_real_particles_ = base_particles->pos_n_.size();
			write_states.WriteToFile(Real(ite));
		}
        
		base_particles->total_real_particles_  = base_particles->pos_n_.size();
		std::cout << base_particles->total_real_particles_ << " Particles has been successfully created!" << "\n" << std::endl;
	}
	//=================================================================================================//
}
//=================================================================================================//