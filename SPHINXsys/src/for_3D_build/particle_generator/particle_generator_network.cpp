/**
 * @file 	particle_generator_network.cpp
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.2
 *			In this version, a network generator is added. -- Chi ZHANG
 */
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
	ParticleGeneratorNetwork::ParticleGeneratorNetwork(Vecd starting_pnt, Vecd second_pnt)
		: ParticleGenerator(), starting_pnt_(starting_pnt), second_pnt_(second_pnt),
		n_it_(25), fascicles_(true), segments_in_branch_(10), segment_length_(0),
		body_shape_(NULL){}
	//=================================================================================================//
	void ParticleGeneratorNetwork::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		segment_length_ = sph_body_->particle_spacing_;
		body_shape_ = sph_body_->body_shape_;
	}
	//=================================================================================================//
	Vecd ParticleGeneratorNetwork::
		getGradientFromNearestPoints(Point pt, Real delta, BaseMeshCellLinkedList* mesh_cell_linked_list)
	{
		Point upgrad(0), downgrad(0);
		Point shift(delta);
		for (int i = 0; i < pt.size(); i++) {
			Point upwind = pt;
			Point downwind = pt;
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
	Point ParticleGeneratorNetwork::creatNewBranchPoint(Point init_point, Vecd dir)
	{
		Vecd pnt_to_project = init_point + dir * segment_length_;

		Real phi = body_shape_->findSignedDistance(pnt_to_project);
		Vecd unit_normal = body_shape_->findNormalDirection(pnt_to_project);
		unit_normal /= unit_normal.norm() + TinyReal;

		return pnt_to_project - phi * unit_normal;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::isCollision(Point& new_point, 
		ListData& nearest_neighbor, size_t parent_id, Tree* tree)
	{
		bool collision = false;
		bool is_family = false;

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
			if (min_distance < segment_length_) collision = true;
		}

		return collision;
	}
	//=================================================================================================//
	bool ParticleGeneratorNetwork::createABranchIfValid(SPHBody* sph_body, size_t parent_id, Real angle,
		Real repulsivity, size_t number_segments, Tree* tree)
	{
		bool is_valid = false;
		Branch* parent_branch = tree->branches_[parent_id];
		IndexVector& parent_elements = parent_branch->elements_;
		StdVec<Point>& points = tree->points_;
	
		Point init_point = points[parent_elements.back()];
		Vecd init_direction = parent_branch->end_direction_;

		BaseMeshCellLinkedList* mesh_cell_linked_list = sph_body->mesh_cell_linked_list_;

		Real delta = 5.0 * segment_length_;

		Vecd surface_norm = body_shape_->findNormalDirection(init_point);
		surface_norm /= surface_norm.norm() + TinyReal;
		Vecd in_plane = - SimTK::cross(init_direction, surface_norm);

		Vecd dir = cos(angle) * init_direction + sin(angle) * in_plane;
		dir /= dir.norm() + TinyReal;
		Vecd grad = getGradientFromNearestPoints(init_point, delta, mesh_cell_linked_list);
		Vecd end_direction = (repulsivity * grad + dir) / ((repulsivity * grad + dir).norm() + TinyReal);
		Point end_point = init_point;

		Point new_point = creatNewBranchPoint(end_point, end_direction);
		ListData nearest_neighbor = mesh_cell_linked_list->findNearestListDataEntry(new_point);
		if (!isCollision(new_point, nearest_neighbor, parent_id, tree)) {
			is_valid = true;
			Branch* new_branch = new Branch(parent_id, tree);
			tree->addANewBranch(new_branch);
			tree->addANewBranchElement(new_branch, new_point, end_direction);

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

				new_point = creatNewBranchPoint(end_point, end_direction);
				ListData nearest_neighbor = mesh_cell_linked_list->findNearestListDataEntry(new_point);
				if (isCollision(new_point, nearest_neighbor, parent_id, tree))
				{
					new_branch->is_end_ = true;
					std::cout << "Branch Collision Detected, Break! " << std::endl;
					break;
				}
				tree->addANewBranchElement(new_branch, new_point, end_direction);
			}

			for (size_t i = 0; i != new_branch->elements_.size(); i++)
			{
				mesh_cell_linked_list->InsertACellLinkedListDataEntry(new_branch->elements_[i], points[new_branch->elements_[i]]);
				sph_body->base_particles_->initializeABaseParticle(points[new_branch->elements_[i]], segment_length_);
			}
		}

		return is_valid;
	}
	//=================================================================================================//
	void ParticleGeneratorNetwork::CreateBaseParticles(BaseParticles* base_particles)
	{
		In_Output in_output(sph_body_->getSPHSystem());
		WriteBodyStatesToVtu 	write_states(in_output, { sph_body_ });

		Tree my_tree(starting_pnt_, second_pnt_);
		StdVec<Branch*>& branches = my_tree.branches_;
		StdVec<Point>& points = my_tree.points_;

		base_particles->initializeABaseParticle(starting_pnt_, segment_length_);
		BaseMeshCellLinkedList* mesh_cell_linked_list = sph_body_->mesh_cell_linked_list_;
		mesh_cell_linked_list->InsertACellLinkedListDataEntry(0, points[0]);
		std::cout << "Now creating Particles on network... " << "\n" << std::endl;

		//the first branch
		bool is_valid = createABranchIfValid(sph_body_, 0, 0.0, 0.0, segments_in_branch_, &my_tree);

		size_t ite = 0;
		sph_body_->setNewlyUpdated();
		sph_body_->number_of_particles_ = base_particles->pos_n_.size();
		write_states.WriteToFile(0.0);

		IndexVector edges_to_grow;
		IndexVector new_edges_to_grow;
        if (is_valid) edges_to_grow.push_back(my_tree.last_branch_id_);
	
		if (fascicles_)
		{
			/** Set vertices in family edge. */
			edges_to_grow.clear();
			for (int i = 0; i < 2; i++)
			{
				/** Creating a new edge. */
				Real  angle_to_use = fascicle_angles_[i];
				size_t fascicles_segments = int(fascicle_ratio_ * segments_in_branch_);
				bool is_valid = createABranchIfValid(sph_body_, 1, angle_to_use, 0.0, fascicles_segments, &my_tree);
				if (is_valid) edges_to_grow.push_back(my_tree.last_branch_id_);
			}

			ite++;
			sph_body_->setNewlyUpdated();
			sph_body_->number_of_particles_ = base_particles->pos_n_.size();
			write_states.WriteToFile(Real(ite));

		}

		for(size_t i = 0; i != n_it_; i++)
        {
			
			new_edges_to_grow.clear();
			random_shuffle(edges_to_grow.begin(), edges_to_grow.end());
            for(size_t j = 0; j != edges_to_grow.size(); j++)
            {
                size_t grow_id = edges_to_grow[j];
				Real angle_to_use = -2.0 * angle_ * (((double)rand() / (RAND_MAX)) - 0.5);
                for(size_t k = 0; k != 2; k++)
                {   
                    /** Creating a new edge. */
                    size_t random_number_segments = segments_in_branch_ + rand() % 10 + 1;
					bool is_valid = createABranchIfValid(sph_body_, grow_id, angle_to_use, repulsivity_,
						random_number_segments, &my_tree);

					if(is_valid && !branches[my_tree.last_branch_id_]->is_end_)
                    {
                        new_edges_to_grow.push_back(my_tree.last_branch_id_);
                    }

                    angle_to_use *= -1.0;
                }
            }
            edges_to_grow = new_edges_to_grow;

			ite++;
			sph_body_->setNewlyUpdated();
			sph_body_->number_of_particles_ = base_particles->pos_n_.size();
			write_states.WriteToFile(Real(ite));
		}
        
		sph_body_->number_of_particles_  = base_particles->pos_n_.size();
		std::cout << sph_body_->number_of_particles_ << " Particles has been successfully created!" << "\n" << std::endl;
	}
	//=================================================================================================//
}
//=================================================================================================//