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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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
#include "particle_reserve.h"

namespace SPH
{
template <>
class Ghost<PeriodicAlongAxis> : public Ghost<Base>, public PeriodicAlongAxis
{
  public:
    Ghost(BoundingBoxd bounding_bounds, int axis);
    virtual ~Ghost() {};
    void reserveGhostParticles(BaseParticles &base_particles, Real particle_spacing);
    ParticlesBound &LowerGhostBound() { return lower_ghost_bound_; };
    ParticlesBound &UpperGhostBound() { return upper_ghost_bound_; };

  private:
    ParticlesBound lower_ghost_bound_ = {0, 0};
    ParticlesBound upper_ghost_bound_ = {0, 0};
    size_t calculateGhostSize(Real particle_spacing);
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
        Ghost<PeriodicAlongAxis> &ghost_boundary_;
        std::pair<size_t, size_t> &lower_ghost_bound_;
        std::pair<size_t, size_t> &upper_ghost_bound_;
        BaseCellLinkedList &cell_linked_list_;
        Real *Vol_;

        virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
        virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

      public:
        CreatPeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                    Ghost<PeriodicAlongAxis> &ghost_boundary);
        virtual ~CreatPeriodicGhostParticles() {};

        virtual void setupDynamics(Real dt = 0.0) override;
        virtual void exec(Real dt = 0.0) override;
    };

    /**
     * @class UpdatePeriodicGhostParticles
     * @brief update ghost particles in an axis direction
     */
    class UpdatePeriodicGhostParticles : public PeriodicBounding
    {
      protected:
        Ghost<PeriodicAlongAxis> &ghost_boundary_;
        std::pair<size_t, size_t> &lower_ghost_bound_;
        std::pair<size_t, size_t> &upper_ghost_bound_;
        UnsignedInt *sorted_id_;

        void checkLowerBound(size_t index_i, Real dt = 0.0) override;
        void checkUpperBound(size_t index_i, Real dt = 0.0) override;

      public:
        UpdatePeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                     Ghost<PeriodicAlongAxis> &ghost_boundary);
        virtual ~UpdatePeriodicGhostParticles() {};

        virtual void exec(Real dt = 0.0) override;
    };

  public:
    PeriodicConditionUsingGhostParticles(RealBody &real_body, Ghost<PeriodicAlongAxis> &ghost_boundary);
    virtual ~PeriodicConditionUsingGhostParticles() {};

    PeriodicBounding bounding_;
    CreatPeriodicGhostParticles ghost_creation_;
    UpdatePeriodicGhostParticles ghost_update_;
};
} // namespace SPH
#endif // GHOST_BOUNDING_H
