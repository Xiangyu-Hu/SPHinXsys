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
 * @file 	ghost_bounding.h
 * @brief This is the particle dynamics for domain bounding
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef GHOST_BOUNDING_H
#define GHOST_BOUNDING_H

#include "domain_bounding.h"

namespace SPH
{
template <>
class Ghost<PeriodicAlongAxis> : public PeriodicAlongAxis
{
  public:
    Ghost(BoundingBox bounding_bounds, int axis)
        : PeriodicAlongAxis(bounding_bounds, axis){};
    virtual ~Ghost(){};

    void reserveGhostParticle(BaseParticles &base_particles, Real particle_spacing)
    {
        size_t ghost_size = getGhostSize(particle_spacing);
        lower_ghost_bound_.first = base_particles.addGhostParticles(ghost_size);
        upper_ghost_bound_.first = base_particles.addGhostParticles(ghost_size);
    };

    std::pair<size_t, size_t> &LowerGhostBound() { return lower_ghost_bound_; };
    std::pair<size_t, size_t> &UpperGhostBound() { return upper_ghost_bound_; };

  private:
    std::pair<size_t, size_t> lower_ghost_bound_;
    std::pair<size_t, size_t> upper_ghost_bound_;

    size_t getGhostSize(Real particle_spacing)
    {
        int next_axis = NextAxis(axis_);
        Real bound_size = bounding_bounds_.second_[next_axis] - bounding_bounds_.first_[next_axis];
        Real ghost_width = 4.0;
        return std::ceil(2.0 * ghost_width * ABS(bound_size) / particle_spacing);
    };
};

/**
 * @class PeriodicConditionUsingGhostParticles
 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
 *	The first step is carried out before update cell linked list and
 *	the second and third after the updating.
 *  Note that, currently, one should use this class for periodic condition in single direction.
 *  More work is required so that the class works for periodic condition in combined directions,
 *  such as periodic condition in both x and y directions.
 */
class PeriodicConditionUsingGhostParticles : public BasePeriodicCondition<execution::ParallelPolicy>
{
  protected:
    /**
     * @class CreatPeriodicGhostParticles
     * @brief create ghost particles in an axis direction
     */
    class CreatPeriodicGhostParticles : public PeriodicBounding
    {
      protected:
        std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
        StdVec<std::pair<size_t, size_t>> &ghost_ranges_;
        StdLargeVec<Real> &Vol_;
        virtual void setupDynamics(Real dt = 0.0) override;
        virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
        virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

      public:
        CreatPeriodicGhostParticles(Vecd &periodic_translation,
                                    StdVec<CellLists> &bound_cells_data,
                                    StdVec<std::pair<size_t, size_t>> &ghost_ranges,
                                    RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicBounding(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
              ghost_ranges_(ghost_ranges), Vol_(particles_->Vol_){};
        virtual ~CreatPeriodicGhostParticles(){};
    };

    /**
     * @class UpdatePeriodicGhostParticles
     * @brief update ghost particles in an axis direction
     */
    class UpdatePeriodicGhostParticles : public PeriodicBounding
    {
      protected:
        StdVec<std::pair<size_t, size_t>> &ghost_ranges_;
        void checkLowerBound(size_t index_i, Real dt = 0.0) override;
        void checkUpperBound(size_t index_i, Real dt = 0.0) override;

      public:
        UpdatePeriodicGhostParticles(Vecd &periodic_translation,
                                     StdVec<CellLists> &bound_cells_data,
                                     StdVec<std::pair<size_t, size_t>> &ghost_ranges,
                                     RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicBounding(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
              ghost_ranges_(ghost_ranges){};
        virtual ~UpdatePeriodicGhostParticles(){};

        virtual void exec(Real dt = 0.0) override;
    };

  public:
    PeriodicConditionUsingGhostParticles(RealBody &real_body, Ghost<PeriodicAlongAxis> &ghost_boundary)
        : BasePeriodicCondition<execution::ParallelPolicy>(real_body, ghost_boundary.BoundingBounds(), ghost_boundary.Axis()),
          bounding_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis),
          ghost_creation_(periodic_translation_, bound_cells_data_, ghost_ranges_, real_body, bounding_bounds, axis),
          ghost_update_(periodic_translation_, bound_cells_data_, ghost_ranges_, real_body, bounding_bounds, axis)
    {
        ghost_ranges_.resize(2);
    };

    virtual ~PeriodicConditionUsingGhostParticles(){};

    PeriodicBounding bounding_;
    CreatPeriodicGhostParticles ghost_creation_;
    UpdatePeriodicGhostParticles ghost_update_;
};
} // namespace SPH
#endif // GHOST_BOUNDING_H
