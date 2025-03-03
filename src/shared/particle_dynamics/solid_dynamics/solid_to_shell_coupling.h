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
 * @file solid_to_shell_coupling.h
 * @brief Constraints for solid to shell coupling.
 * @author	Xiangyu Hu, Weiyi Kong (Virtonomy GmbH)
 */

#ifndef SOLID_TO_SHELL_CONSTRAINT_H
#define SOLID_TO_SHELL_CONSTRAINT_H

#include "force_prior.h"
#include "general_constraint.h"

namespace SPH
{
namespace solid_dynamics
{
/**@class CouplingPart
 * @brief Find the overlapping particles and set is_coupled to 1.
 * The overlapping criteria is set to distance < factor * average_spacing.
 */
class CouplingPart : public BodyPartByParticle,
                     public DataDelegateContact
{
  private:
    Real distance_factor_;
    int *is_coupled_;

  public:
    CouplingPart(BaseContactRelation &contact_relation, const std::string &body_part_name, Real distance_factor = 0.5);

  private:
    void tagManually(size_t index_i);
};

/**@class NearestNeighborSolidVelocityConstraint
 * @brief Constrain solid velocity by the nearest neighbor shell particle velocity.
 * When there is more than one contact body, use the nearest neighbor shell particle velocity among all contact bodies.
 */
class NearestNeighborSolidVelocityConstraint : public MotionConstraint<BodyPartByParticle>,
                                               public DataDelegateContact
{
  private:
    /**@class NearestNeighborSolidVelocityConstraint::SearchNearestParticle
     * @brief Find the id of the nearest neighbor shell particle among all contact bodies.
     */
    class SearchNearestParticle : public BaseLocalDynamics<BodyPartByParticle>,
                                  public DataDelegateContact
    {
      private:
        StdVec<int *> contact_is_coupled_;
        StdLargeVec<std::pair<size_t, size_t>> coupling_ids_;

      public:
        SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
        void update(size_t index_i, Real dt = 0.0);
        std::pair<size_t, size_t> get_nearest_id(size_t index_i) { return coupling_ids_[index_i]; };
    };

  public:
    NearestNeighborSolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
    void init() { search_nearest_particle_.exec(); }
    void update(size_t index_i, Real dt = 0.0);
    std::pair<size_t, size_t> get_nearest_id(size_t index_i) { return search_nearest_particle_.get_nearest_id(index_i); };

  private:
    SimpleDynamics<SearchNearestParticle> search_nearest_particle_;
    StdVec<Vecd *> contact_vel_;
};

/**@class NearestNeighborShellForceConstraint
 * @brief Apply the coupling force on shell from the nearest neighbor solid particles.
 * When there is more than one contact body, sum up the coupling forces from all contact bodies.
 */
class NearestNeighborShellForceConstraint : public BaseForcePrior<BodyPartByParticle>,
                                            public DataDelegateContact
{
    /**@class NearestNeighborShellForceConstraint::SearchNearestParticle
     * @brief Find the id of the nearest neighbor solid particle for each contact body.
     */
    class SearchNearestParticle : public BaseLocalDynamics<BodyPartByParticle>,
                                  public DataDelegateContact
    {
      private:
        StdVec<int *> contact_is_coupled_;
        StdVec<StdLargeVec<size_t>> coupling_ids_;

      public:
        SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
        void update(size_t index_i, Real dt = 0.0);
        size_t get_nearest_id(size_t index_i, size_t k) { return coupling_ids_[k][index_i]; };
    };

  public:
    NearestNeighborShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
    void init() { search_nearest_particle_.exec(); }
    void interaction(size_t index_i, Real dt = 0.0);
    size_t get_nearest_id(size_t index_i, size_t k) { return search_nearest_particle_.get_nearest_id(index_i, k); };

  private:
    SimpleDynamics<SearchNearestParticle> search_nearest_particle_;
    StdVec<Vecd *> contact_force_;
};

/**@class InterpolationSolidVelocityConstraint
 * @brief Constrain solid velocity by the interpolation from coupled shell particle velocities.
 * When there is more than one contact body, interpolate from all contact bodies.
 */
class InterpolationSolidVelocityConstraint : public MotionConstraint<BodyPartByParticle>,
                                             public DataDelegateContact
{
  private:
    StdVec<int *> contact_is_coupled_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_;

    StdLargeVec<Real> weight_; // weight only needs to be calculated once

  public:
    InterpolationSolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
    void update(size_t index_i, Real dt = 0.0);
    void init() { compute_weight(); }

  private:
    void compute_weight();
};

/**@class InterpolationShellForceConstraint
 * @brief Apply the coupling force on shell by the interpolation from coupled shell particle internal accelerations.
 * When there is more than one contact body, interpolate force from each contact body and sum up.
 */
class InterpolationShellForceConstraint : public BaseForcePrior<BodyPartByParticle>,
                                          public DataDelegateContact
{
  private:
    StdVec<int *> contact_is_coupled_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_force_;

    StdVec<StdLargeVec<Real>> weight_; // weight of each body only needs to be calculated once

  public:
    InterpolationShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation);
    void interaction(size_t index_i, Real dt = 0.0);
    void init() { compute_weight(); }

  private:
    void compute_weight();
};
} // namespace solid_dynamics
} // namespace SPH

#endif // SOLID_TO_SHELL_CONSTRAINT_H