/**
 * @file 	generative_structures.cpp
 * @author	Xiangyu Hu
 */

#include "generative_structures.h"

#include "base_body.h"
#include "base_particles.h"
#include "adaptation.h"

namespace SPH
{
	//=================================================================================================//
	GenerativeStructure::GenerativeStructure(SPHBody *sph_body)
		: sph_body_(sph_body),
		  spacing_ref_(sph_body_->sph_adaptation_->ReferenceSpacing()),
		  base_particles_(sph_body->base_particles_),
		  neighbor_relation_inner_(sph_body),
		  pos_n_(base_particles_->pos_n_),
		  Vol_(base_particles_->Vol_){};
	//=================================================================================================//
	GenerativeTree::GenerativeTree(SPHBody *sph_body)
		: GenerativeStructure(sph_body), last_branch_id_(0)
	{
		root_ =  branches_ptr_keeper_.createPtr<Branch>(this);
	}
	//=================================================================================================//
	void GenerativeTree::buildParticleConfiguration(BaseParticles &base_particles,
													ParticleConfiguration &particle_configuration)
	{
		size_t particle_id;
		size_t parent_branch_id;
		size_t child_branch_id;
		size_t num_ele;
		std::vector<size_t> neighboring_ids;
		std::vector<size_t> child_ids;
		/** First branch
		 * Note that the first branch has only one particle.
		 * Find the neighbors in child branch, the first branch only have one child, id = 1.
		 */
		particle_id = branches_[0]->inner_particles_.front();
		neighboring_ids.clear();
		neighboring_ids.push_back(branches_[1]->inner_particles_[0]);
		neighboring_ids.push_back(branches_[1]->inner_particles_[1]);
		/** Build configuration. */
		Neighborhood &neighborhood = particle_configuration[particle_id];
		for (size_t n = 0; n != neighboring_ids.size(); ++n)
		{
			Vecd displacement = base_particles.pos_n_[particle_id] - base_particles.pos_n_[neighboring_ids[n]];
			neighbor_relation_inner_(neighborhood, displacement, particle_id, neighboring_ids[n]);
		}
		/** Second branch. 
		 * The second branch has special parent branch, branch 0, consisting only one point.
		 * The child branch are two normal branch. 
		 */
		num_ele = branches_[1]->inner_particles_.size();
		child_ids.clear();
		for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
		{
			child_ids.push_back(branches_[1]->out_edge_[k]);
		}

		for (size_t i = 0; i != num_ele; i++)
		{
			neighboring_ids.clear();
			particle_id = branches_[1]->inner_particles_.front() + i;
			if (i == 0)
			{
				neighboring_ids.push_back(branches_[0]->inner_particles_.front());
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (i == 1)
			{
				neighboring_ids.push_back(branches_[0]->inner_particles_.front());
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (2 <= i && i <= (num_ele - 3))
			{
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id - 2);
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (i == (num_ele - 2))
			{
				neighboring_ids.push_back(particle_id - 2);
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id + 1);

				for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
				{
					child_branch_id = branches_[1]->out_edge_[k];
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
				}
			}
			else if (i == (num_ele - 1))
			{
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id - 2);

				for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
				{
					child_branch_id = branches_[1]->out_edge_[k];
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
				}
			}

			Neighborhood &neighborhood = particle_configuration[particle_id];
			for (size_t n = 0; n != neighboring_ids.size(); ++n)
			{
				Vecd displacement = base_particles.pos_n_[particle_id] - base_particles.pos_n_[neighboring_ids[n]];
				neighbor_relation_inner_(neighborhood, displacement, particle_id, neighboring_ids[n]);
			}
		}
		/** Other branches. 
		 * They are may normal branch (fully growed, has child and parent) or non-fully growed branch
		 */
		for (size_t branch_idx = 2; branch_idx != branches_.size(); ++branch_idx)
		{
			num_ele = branches_[branch_idx]->inner_particles_.size();
			parent_branch_id = branches_[branch_idx]->in_edge_;
			if (!branches_[branch_idx]->is_terminated_)
			{
				/** This branch is fully growed. */
				for (size_t i = 0; i != num_ele; i++)
				{
					neighboring_ids.clear();
					particle_id = branches_[branch_idx]->inner_particles_.front() + i;
					if (i == 0)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back() - 1);

						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (i == 1)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (2 <= i && i <= (num_ele - 3))
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);
						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (i == (num_ele - 2))
					{
						neighboring_ids.push_back(particle_id - 2);
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id + 1);

						for (size_t k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
						{
							child_branch_id = branches_[branch_idx]->out_edge_[k];
							neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
						}
					}
					else if (i == (num_ele - 1))
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);

						for (size_t k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
						{
							child_branch_id = branches_[branch_idx]->out_edge_[k];
							neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
							if (branches_[child_branch_id]->inner_particles_.size() >= 2)
							{
								neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
							}
						}
					}

					Neighborhood &neighborhood = particle_configuration[particle_id];
					for (size_t n = 0; n != neighboring_ids.size(); ++n)
					{
						Vecd displacement = base_particles.pos_n_[particle_id] - base_particles.pos_n_[neighboring_ids[n]];
						neighbor_relation_inner_(neighborhood, displacement, particle_id, neighboring_ids[n]);
					}
				}
			}
			else
			{
				/** This branch is not fully growed. */
				for (size_t i = 0; i != num_ele; i++)
				{
					neighboring_ids.clear();
					particle_id = branches_[branch_idx]->inner_particles_.front() + i;
					if (i == 0)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						if (branches_[parent_branch_id]->inner_particles_.size() >= 2)
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back() - 1);
					}
					else if (i == 1)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(particle_id - 1);
					}
					else
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);
					}

					if (i + 1 < num_ele)
						neighboring_ids.push_back(particle_id + 1);
					if (i + 2 < num_ele)
						neighboring_ids.push_back(particle_id + 2);

					Neighborhood &neighborhood = particle_configuration[particle_id];
					for (size_t n = 0; n != neighboring_ids.size(); ++n)
					{
						Vecd displacement = base_particles.pos_n_[particle_id] - base_particles.pos_n_[neighboring_ids[n]];
						neighbor_relation_inner_(neighborhood, displacement, particle_id, neighboring_ids[n]);
					}
				}
			}
		}
	}
	//=================================================================================================//
	void GenerativeTree::
		growAParticleOnBranch(Branch *branch, const Vecd &new_point, const Vecd &end_direction)
	{
		base_particles_->initializeABaseParticle(new_point, spacing_ref_);
		branch_locations_.push_back(branch->id_);
		branch->inner_particles_.push_back(pos_n_.size() - 1);
		branch->end_direction_ = end_direction;
	}
	//=================================================================================================//
	size_t GenerativeTree::BranchLocation(size_t particle_idx)
	{
		return particle_idx < pos_n_.size() ? branch_locations_[particle_idx] : MaxSize_t;
	}
	//=================================================================================================//
	GenerativeTree::Branch::Branch(GenerativeTree *tree)
		: Edge<size_t, IndexVector>(tree), is_terminated_(false)
	{
		tree->branches_.push_back(this);
		tree->last_branch_id_ = id_;
	}
	//=================================================================================================//
	GenerativeTree::Branch::Branch(size_t parent_id, GenerativeTree *tree)
		: Edge<size_t, IndexVector>(parent_id, tree), is_terminated_(false)
	{
		tree->branches_[parent_id]->out_edge_.push_back(id_);
		tree->branches_.push_back(this);
		tree->last_branch_id_ = id_;
	}
	//=================================================================================================//
}
