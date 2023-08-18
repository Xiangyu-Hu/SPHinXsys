/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	complex_body.h
 * @brief 	A complex body is characterized with a secondary structure,
 * 			which can be imported externally or created according to specific rules.
 * 			The secondary structure will be used or even created by the corresponding
 * 			particle generator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef COMPLEX_BODY_H
#define COMPLEX_BODY_H

#include "base_body.h"
#include "base_body_part.h"
#include "base_geometry.h"

namespace SPH
{
/**
 * @class SecondaryStructure
 * @brief Abstract class as interface for all secondary structures.
 * Currently, it provides interface on building inner configuration.
 * The interface can be extended.
 */
class SecondaryStructure
{
  public:
    explicit SecondaryStructure(){};
    virtual ~SecondaryStructure(){};

    virtual void buildParticleConfiguration(ParticleConfiguration &particle_configuration) = 0;
};

/**
 * @class TreeBody
 * @brief The tree is composed of a root (the first branch)
 * and other branch generated sequentially.
 */
class TreeBody : public SecondaryStructure, public RealBody
{
  public:
    class Branch;

  private:
    UniquePtrsKeeper<Branch> branches_ptr_keeper_;

  public:
    StdVec<Branch *> branches_;    /**< Container of all branches */
    IndexVector branch_locations_; /**< in which branch are the particles located */
    size_t last_branch_id_;
    Branch *root_;

    template <typename... ConstructorArgs>
    TreeBody(ConstructorArgs &&...args)
        : SecondaryStructure(), RealBody(std::forward<ConstructorArgs>(args)...), last_branch_id_(0)
    {
        root_ = branches_ptr_keeper_.createPtr<Branch>(this);
    };
    virtual ~TreeBody(){};

    Branch *createANewBranch(size_t parent_id)
    {
        return branches_ptr_keeper_.createPtr<Branch>(parent_id, this);
    };
    size_t BranchLocation(size_t particle_idx);
    Branch *LastBranch() { return branches_[last_branch_id_]; };

    virtual void buildParticleConfiguration(ParticleConfiguration &particle_configuration) override;
    size_t ContainerSize() { return branches_.size(); };
};

/**
 * @class TreeBody::Branch
 * @brief Each branch (except the root) has a parent and several children, and geometric information.
 * It is a realized edge and has multi inner particles.
 * The first is the last particle from the parent or root,
 * and the last is the first particle of all its child branches.
 * Many connected branches compose a tree. */
class TreeBody::Branch : public Edge<size_t, IndexVector>
{
  public:
    /** construct the root branch  */
    explicit Branch(TreeBody *tree);
    /** construct an branch connecting with its parent */
    Branch(size_t parent_id, TreeBody *tree);
    virtual ~Branch(){};

    Vecd end_direction_; /**< the direction pointing to the last particle */
    /** The indexes of particle within this branch.
     * The first is the last particle from the parent or root,
     * and the last is the first of all its child branches. */
    IndexVector inner_particles_;
    bool is_terminated_; /**< whether is an terminate branch or not */
};

/**
 * @class TreeTerminates
 * @brief A  body part with the collection of particles as the terminates of the tree.
 */
class TreeTerminates : public BodyPartByParticle
{
  public:
    TreeBody &tree_;

    explicit TreeTerminates(SPHBody &sph_body);
    virtual ~TreeTerminates(){};
};
} // namespace SPH
#endif // COMPLEX_BODY_H