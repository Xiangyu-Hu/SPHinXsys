#include "particle_generator_network.h"
#include "all_io.h"
#include "base_body.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "level_set.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<BaseParticles, Network>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, const Vecd &starting_pnt,
                      const Vecd &second_pnt, int iterator, Real grad_factor)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles),
      starting_pnt_(starting_pnt), second_pnt_(second_pnt),
      n_it_(iterator), fascicles_(true), segments_in_branch_(10),
      segment_length_(sph_body.getSPHAdaptation().ReferenceSpacing()),
      grad_factor_(grad_factor), sph_body_(sph_body), initial_shape_(sph_body.getInitialShape()),
      cell_linked_list_(DynamicCast<RealBody>(this, sph_body).getCellLinkedList()),
      tree_(DynamicCast<TreeBody>(this, &sph_body))
{
    Vecd displacement = second_pnt_ - starting_pnt_;
    Vecd end_direction = displacement / (displacement.norm() + TinyReal);
    /** Add initial particle to the first branch of the tree. */
    growAParticleOnBranch(tree_->root_, starting_pnt_, end_direction);
    cell_linked_list_.InsertListDataEntry(0, position_[0]);
}
//=================================================================================================//
void ParticleGenerator<BaseParticles, Network>::
    growAParticleOnBranch(TreeBody::Branch *branch, const Vecd &new_point, const Vecd &end_direction)
{
    addPositionAndVolumetricMeasure(new_point, segment_length_);
    tree_->branch_locations_.push_back(branch->id_);
    branch->inner_particles_.push_back(position_.size() - 1);
    branch->end_direction_ = end_direction;
}
//=================================================================================================//
Vecd ParticleGenerator<BaseParticles, Network>::getGradientFromNearestPoints(Vecd pt, Real delta)
{
    Vecd up_grad = Vecd::Zero();
    Vecd down_grad = Vecd::Zero();
    Vecd shift = delta * Vecd::Ones();

    for (int i = 0; i != Dimensions; i++)
    {
        Vecd upwind = pt;
        Vecd downwind = pt;
        upwind[i] -= shift[i];
        downwind[i] += shift[i];
        ListData up_nearest_list = cell_linked_list_.findNearestListDataEntry(upwind);
        ListData down_nearest_list = cell_linked_list_.findNearestListDataEntry(downwind);
        up_grad[i] = std::get<0>(up_nearest_list) != MaxUnsignedInt
                         ? (upwind - std::get<1>(up_nearest_list)).norm() / 2.0 * delta
                         : 1.0;
        down_grad[i] = std::get<0>(down_nearest_list) != MaxUnsignedInt
                           ? (downwind - std::get<1>(down_nearest_list)).norm() / 2.0 * delta
                           : 1.0;
    }
    return down_grad - up_grad;
}
//=================================================================================================//
Vecd ParticleGenerator<BaseParticles, Network>::createATentativeNewBranchPoint(Vecd init_point, Vecd dir)
{
    Vecd pnt_to_project = init_point + dir * segment_length_;

    Real phi = initial_shape_.findSignedDistance(pnt_to_project);
    Vecd unit_normal = initial_shape_.findNormalDirection(pnt_to_project);
    unit_normal /= unit_normal.norm() + TinyReal;
    Vecd new_point = pnt_to_project - phi * unit_normal;
    return new_point;
}
//=================================================================================================//
bool ParticleGenerator<BaseParticles, Network>::
    isCollision(const Vecd &new_point, const ListData &nearest_neighbor, size_t parent_id)
{
    bool collision = false;
    bool is_family = false;

    collision = extraCheck(new_point);

    size_t edge_location = tree_->BranchLocation(position_.size(), std::get<0>(nearest_neighbor));
    if (edge_location == parent_id)
        is_family = true;
    for (const size_t &brother_branch : tree_->branches_[parent_id]->out_edge_)
    {
        if (edge_location == brother_branch)
        {
            is_family = true;
        }
    }

    if (!is_family)
    {
        Real min_distance = (new_point - std::get<1>(nearest_neighbor)).norm();
        if (min_distance < 5.0 * segment_length_)
            collision = true;
    }

    return collision;
}
//=================================================================================================//
bool ParticleGenerator<BaseParticles, Network>::
    createABranchIfValid(size_t parent_id, Real angle, Real repulsivity, size_t number_segments)
{
    bool is_valid = false;
    TreeBody::Branch *parent_branch = tree_->branches_[parent_id];
    IndexVector &parent_elements = parent_branch->inner_particles_;

    Vecd init_point = position_[parent_elements.back()];
    Vecd init_direction = parent_branch->end_direction_;

    Vecd surface_norm = initial_shape_.findNormalDirection(init_point);
    surface_norm /= surface_norm.norm() + TinyReal;
    Vecd in_plane = -init_direction.cross(surface_norm);

    Real delta = grad_factor_ * segment_length_;
    Vecd grad = getGradientFromNearestPoints(init_point, delta);
    Vecd dir = cos(angle) * init_direction + sin(angle) * in_plane;
    dir /= dir.norm() + TinyReal;
    Vecd end_direction = (repulsivity * grad + dir) / ((repulsivity * grad + dir).norm() + TinyReal);
    Vecd end_point = init_point;

    Vecd new_point = createATentativeNewBranchPoint(end_point, end_direction);
    if (!isCollision(new_point, cell_linked_list_.findNearestListDataEntry(new_point), parent_id))
    {
        is_valid = true;
        TreeBody::Branch *new_branch = tree_->createANewBranch(parent_id);
        growAParticleOnBranch(new_branch, new_point, end_direction);

        for (size_t i = 1; i < number_segments; i++)
        {
            surface_norm = initial_shape_.findNormalDirection(new_point);
            surface_norm /= surface_norm.norm() + TinyReal;
            /** Project grad to surface. */
            grad = getGradientFromNearestPoints(new_point, delta);
            grad -= grad.dot(surface_norm) * surface_norm;
            dir = (repulsivity * grad + end_direction) / ((repulsivity * grad + end_direction).norm() + TinyReal);
            end_direction = dir;
            end_point = new_point;

            new_point = createATentativeNewBranchPoint(end_point, end_direction);
            if (isCollision(new_point, cell_linked_list_.findNearestListDataEntry(new_point), parent_id))
            {
                new_branch->is_terminated_ = true;
                std::cout << "Branch Collision Detected, Break! " << std::endl;
                break;
            }
            /** This constraint imposed to avoid too small time step size. */
            if ((new_point - end_point).norm() < 0.5 * segment_length_)
            {
                new_branch->is_terminated_ = true;
                std::cout << "New branch point is too close, Break! " << std::endl;
                break;
            }
            growAParticleOnBranch(new_branch, new_point, end_direction);
        }

        for (const size_t &particle_idx : new_branch->inner_particles_)
        {
            cell_linked_list_.InsertListDataEntry(particle_idx, position_[particle_idx]);
        }
    }

    return is_valid;
}
//=================================================================================================//
void ParticleGenerator<BaseParticles, Network>::prepareGeometricData()
{
    ParticleGenerationRecordingToVtp write_particle_generation(sph_body_, position_);

    std::cout << "Now creating Particles on network... " << std::endl;

    size_t ite = 0;
    sph_body_.setNewlyUpdated();
    write_particle_generation.writeToFile(0);

    IndexVector branches_to_grow;
    IndexVector new_branches_to_grow;

    // the second branch
    if (createABranchIfValid(0, 0.0, 0.0, segments_in_branch_))
    {
        branches_to_grow.push_back(tree_->last_branch_id_);
    }

    if (fascicles_)
    {
        /** Set vertices in family branch. */
        branches_to_grow.clear();
        for (size_t i = 0; i != 2; i++)
        {
            /** Creating a new branch. */
            Real angle_to_use = fascicle_angles_[i];
            size_t fascicles_segments = int(fascicle_ratio_ * segments_in_branch_);
            if (createABranchIfValid(1, angle_to_use, 0.0, fascicles_segments))
            {
                branches_to_grow.push_back(tree_->last_branch_id_);
            }
        }

        ite++;
        sph_body_.setNewlyUpdated();
        write_particle_generation.writeToFile(ite);
    }
    std::mt19937_64 random_engine;
    for (size_t i = 0; i != n_it_; i++)
    {
        new_branches_to_grow.clear();
        std::shuffle(branches_to_grow.begin(), branches_to_grow.end(), random_engine);
        for (size_t j = 0; j != branches_to_grow.size(); j++)
        {
            size_t grow_id = branches_to_grow[j];
            Real rand_num = rand_uniform(-0.5, 0.5);
            Real angle_to_use = angle_ + rand_num * 0.05;
            for (size_t k = 0; k != 2; k++)
            {
                /** Creating a new branch with fixed number of segments. */
                size_t random_number_segments = segments_in_branch_;
                if (createABranchIfValid(grow_id, angle_to_use, repulsivity_, random_number_segments) &&
                    !tree_->LastBranch()->is_terminated_)
                {
                    new_branches_to_grow.push_back(tree_->last_branch_id_);
                }

                angle_to_use *= -1.0;
            }
        }
        branches_to_grow = new_branches_to_grow;

        ite++;
        sph_body_.setNewlyUpdated();
        write_particle_generation.writeToFile(ite);
    }

    std::cout << position_.size() << " particles has been successfully created!" << std::endl;
}
//=================================================================================================//
} // namespace SPH
