#include "tree_body.h"

#include "adaptation.h"
#include "base_material.h"
#include "base_particle_dynamics.h"
#include "base_particles.hpp"
#include "neighborhood.h"

namespace SPH
{
//=================================================================================================//
void TreeBody::buildParticleConfiguration(ParticleConfiguration &particle_configuration)
{
    /** First branch
     * Note that the first branch has only one particle.
     * Find the neighbors in child branch, the first branch only have one child, id = 1.
     */
    size_t particle_id = branches_[0]->inner_particles_.front();
    std::vector<size_t> neighboring_ids;
    neighboring_ids.push_back(branches_[1]->inner_particles_[0]);
    neighboring_ids.push_back(branches_[1]->inner_particles_[1]);
    /** Build configuration. */
    Vecd *pos = base_particles_->ParticlePositions();
    NeighborBuilderInner neighbor_relation_inner(*this);
    for (size_t n = 0; n != neighboring_ids.size(); ++n)
    {
        size_t index_j = neighboring_ids[n];
        ListData list_data_j = std::make_pair(index_j, pos[index_j]);
        Neighborhood &neighborhood = particle_configuration[particle_id];
        neighbor_relation_inner(neighborhood, pos[particle_id], particle_id, list_data_j);
    }
    /** Second branch.
     * The second branch has special parent branch, branch 0, consisting only one point.
     * The child branch are two normal branch.
     */
    size_t num_ele = branches_[1]->inner_particles_.size();
    std::vector<size_t> child_ids;
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
                size_t child_branch_id = branches_[1]->out_edge_[k];
                neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
            }
        }
        else if (i == (num_ele - 1))
        {
            neighboring_ids.push_back(particle_id - 1);
            neighboring_ids.push_back(particle_id - 2);

            for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
            {
                size_t child_branch_id = branches_[1]->out_edge_[k];
                neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
                neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
            }
        }

        for (size_t n = 0; n != neighboring_ids.size(); ++n)
        {
            size_t index_j = neighboring_ids[n];
            ListData list_data_j = std::make_pair(index_j, pos[index_j]);
            Neighborhood &neighborhood = particle_configuration[particle_id];
            neighbor_relation_inner(neighborhood, pos[particle_id], particle_id, list_data_j);
        }
    }
    /** Other branches.
     * They are may normal branch (fully grown, has child and parent) or non-fully grown branch
     */
    for (size_t branch_idx = 2; branch_idx != branches_.size(); ++branch_idx)
    {
        num_ele = branches_[branch_idx]->inner_particles_.size();
        size_t parent_branch_id = branches_[branch_idx]->in_edge_;
        if (!branches_[branch_idx]->is_terminated_)
        {
            /** This branch is fully grown. */
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
                        size_t child_branch_id = branches_[branch_idx]->out_edge_[k];
                        neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
                    }
                }
                else if (i == (num_ele - 1))
                {
                    neighboring_ids.push_back(particle_id - 1);
                    neighboring_ids.push_back(particle_id - 2);

                    for (size_t k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
                    {
                        size_t child_branch_id = branches_[branch_idx]->out_edge_[k];
                        neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
                        if (branches_[child_branch_id]->inner_particles_.size() >= 2)
                        {
                            neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
                        }
                    }
                }

                for (size_t n = 0; n != neighboring_ids.size(); ++n)
                {
                    size_t index_j = neighboring_ids[n];
                    ListData list_data_j = std::make_pair(index_j, pos[index_j]);
                    Neighborhood &neighborhood = particle_configuration[particle_id];
                    neighbor_relation_inner(neighborhood, pos[particle_id], particle_id, list_data_j);
                }
            }
        }
        else
        {
            /** This branch is not fully grown. */
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

                for (size_t n = 0; n != neighboring_ids.size(); ++n)
                {
                    size_t index_j = neighboring_ids[n];
                    ListData list_data_j = std::make_pair(index_j, pos[index_j]);
                    Neighborhood &neighborhood = particle_configuration[particle_id];
                    neighbor_relation_inner(neighborhood, pos[particle_id], particle_id, list_data_j);
                }
            }
        }
    }
}
//=================================================================================================//
size_t TreeBody::BranchLocation(size_t total_particles, size_t particle_idx)
{
    return particle_idx < total_particles ? branch_locations_[particle_idx] : MaxUnsignedInt;
}
//=================================================================================================//
TreeBody::Branch::Branch(TreeBody *tree)
    : Edge<size_t, IndexVector>(tree), is_terminated_(false)
{
    tree->branches_.push_back(this);
    tree->last_branch_id_ = id_;
}
//=================================================================================================//
TreeBody::Branch::Branch(size_t parent_id, TreeBody *tree)
    : Edge<size_t, IndexVector>(parent_id, tree), is_terminated_(false)
{
    tree->branches_[parent_id]->out_edge_.push_back(id_);
    tree->branches_.push_back(this);
    tree->last_branch_id_ = id_;
}
//=================================================================================================//
TreeTerminates::TreeTerminates(SPHBody &sph_body)
    : BodyPartByParticle(sph_body),
      tree_(DynamicCast<TreeBody>(this, sph_body))
{
    for (const auto *branch : tree_.branches_)
    {
        if (branch->is_terminated_)
        {
            size_t particle_index = branch->inner_particles_.back();
            body_part_particles_.push_back(particle_index);
        }
    }
}
//=================================================================================================//
} // namespace SPH
