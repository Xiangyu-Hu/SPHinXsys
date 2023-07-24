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
 * @file 	base_body_relation.h
 * @brief 	Base classes on body and particle topology relations.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_BODY_RELATION_H
#define BASE_BODY_RELATION_H

#include "base_geometry.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "complex_body.h"
#include "neighborhood.h"

namespace SPH
{
/** a small functor for obtaining search range for the simplest case */
struct SearchDepthSingleResolution
{
    int operator()(size_t particle_index) const { return 1; };
};

/** @brief a small functor for obtaining search depth across resolution
 * @details Note that the search depth is defined on the target cell linked list.
 */
struct SearchDepthContact
{
    int search_depth_;
    SearchDepthContact(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
        : search_depth_(1)
    {
        Real inv_grid_spacing_ = 1.0 / target_cell_linked_list->GridSpacing();
        Kernel *kernel_ = sph_body.sph_adaptation_->getKernel();
        search_depth_ = 1 + (int)floor(kernel_->CutOffRadius() * inv_grid_spacing_);
    };
    int operator()(size_t particle_index) const { return search_depth_; };
};

/** @brief a small functor for obtaining search depth for variable smoothing length
 * @details Note that the search depth is defined on the target cell linked list.
 */
struct SearchDepthAdaptive
{
    Real inv_grid_spacing_;
    Kernel *kernel_;
    StdLargeVec<Real> &h_ratio_;
    SearchDepthAdaptive(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
        : inv_grid_spacing_(1.0 / target_cell_linked_list->GridSpacing()),
          kernel_(sph_body.sph_adaptation_->getKernel()),
          h_ratio_(*sph_body.getBaseParticles().getVariableByName<Real>("SmoothingLengthRatio")){};
    int operator()(size_t particle_index) const
    {
        return 1 + (int)floor(kernel_->CutOffRadius(h_ratio_[particle_index]) * inv_grid_spacing_);
    };
};

/** @brief a small functor for obtaining search depth for variable smoothing length
 * @details Note that this is only for building contact neighbor relation.
 */
struct SearchDepthAdaptiveContact
{
    Real inv_grid_spacing_;
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    SearchDepthAdaptiveContact(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
        : inv_grid_spacing_(1.0 / target_cell_linked_list->GridSpacing()),
          sph_adaptation_(*sph_body.sph_adaptation_),
          kernel_(*sph_body.sph_adaptation_->getKernel()){};
    int operator()(size_t particle_index) const
    {
        return 1 + (int)floor(kernel_.CutOffRadius(sph_adaptation_.SmoothingLengthRatio(particle_index)) * inv_grid_spacing_);
    };
};

/** Transfer body parts to real bodies. **/
RealBodyVector BodyPartsToRealBodies(BodyPartVector body_parts);

/**
 * @class SPHRelation
 * @brief The abstract class for all relations within a SPH body or with its contact SPH bodies
 */
class SPHRelation
{
  protected:
    SPHBody &sph_body_;

  public:
    BaseParticles &base_particles_;
    SPHBody &getSPHBody() { return sph_body_; };
    explicit SPHRelation(SPHBody &sph_body);
    virtual ~SPHRelation(){};

    void subscribeToBody() { sph_body_.body_relations_.push_back(this); };
    virtual void resizeConfiguration() = 0;
    virtual void updateConfiguration() = 0;
};

/**
 * @class BaseInnerRelation
 * @brief The abstract relation within a SPH body
 */
class BaseInnerRelation : public SPHRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize();

  public:
    RealBody *real_body_;
    ParticleConfiguration inner_configuration_; /**< inner configuration for the neighbor relations. */
    explicit BaseInnerRelation(RealBody &real_body);
    virtual ~BaseInnerRelation(){};

    virtual void resizeConfiguration() override;
};

/**
 * @class BaseContactRelation
 * @brief The base relation between a SPH body and its contact SPH bodies
 */
class BaseContactRelation : public SPHRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize();

  public:
    RealBodyVector contact_bodies_;
    StdVec<ParticleConfiguration> contact_configuration_; /**< Configurations for particle interaction between bodies. */

    BaseContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies);
    BaseContactRelation(SPHBody &sph_body, BodyPartVector contact_body_parts)
        : BaseContactRelation(sph_body, BodyPartsToRealBodies(contact_body_parts)){};
    virtual ~BaseContactRelation(){};

    virtual void resizeConfiguration() override;
};
} // namespace SPH
#endif // BASE_BODY_RELATION_H