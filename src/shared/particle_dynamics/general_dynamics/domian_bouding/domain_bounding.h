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
 * @file 	domain_bounding.h
 * @brief This is the particle dynamics for domain bounding
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef DOMAIN_BOUNDING_H
#define DOMAIN_BOUNDING_H

#include "base_general_dynamics.h"

namespace SPH
{
/**
 * @class PeriodicAlongAxis
 * @brief Bounding particle position in along axis.
 * The axis must be 0, 1 for 2d and 0, 1, 2 for 3d
 */
struct PeriodicAlongAxis
{
  protected:
    BoundingBox bounding_bounds_; /**< lower and upper bound for checking. */
    const int axis_;              /**< the axis directions for bounding*/
    Vecd periodic_translation_;

    PeriodicAlongAxis(BoundingBox bounding_bounds, int axis)
        : bounding_bounds_(bounding_bounds), axis_(axis),
          periodic_translation_(Vecd::Zero())
    {
        periodic_translation_[axis] =
            bounding_bounds.second_[axis] - bounding_bounds.first_[axis];
    };
    virtual ~PeriodicAlongAxis(){};
};

/**
 * @class BasePeriodicCondition
 * @brief Base class for two different type periodic boundary conditions.
 */
template <class ExecutionPolicy>
class BasePeriodicCondition
{
  public:
    BasePeriodicCondition(RealBody &real_body, PeriodicAlongAxis &periodic_box)
    {
        bound_cells_data_.resize(2);
        BaseCellLinkedList &cell_linked_list = real_body.getCellLinkedList();
        cell_linked_list.tagBoundingCells(
            bound_cells_data_, periodic_box.bounding_bounds_, periodic_box.axis_);
    };
    virtual ~BasePeriodicCondition(){};

  protected:
    StdVec<CellLists> bound_cells_data_;

    /**
     * @class PeriodicBounding
     * @brief Periodic bounding particle position in an axis direction
     */
    class PeriodicBounding : public BaseDynamics<void>,
                             public LocalDynamics,
    {
      protected:
        BoundingBox bounding_bounds_;
        const int axis_;
        Vecd periodic_translation_;
        StdVec<CellLists> &bound_cells_data_;
        StdLargeVec<Vecd> &pos_;

        virtual void checkLowerBound(size_t index_i, Real dt = 0.0)
        {
            if (pos_[index_i][axis_] < bounding_bounds_.first_[axis_])
                pos_[index_i][axis_] += periodic_translation_[axis_];
        };

        virtual void checkUpperBound(size_t index_i, Real dt = 0.0)
        {
            if (pos_[index_i][axis_] > bounding_bounds_.second_[axis_])
                pos_[index_i][axis_] -= periodic_translation_[axis_];
        };

      public:
        PeriodicBounding(StdVec<CellLists> &bound_cells_data,
                         RealBody &real_body, PeriodicAlongAxis &periodic_box)
            : BaseDynamics<void>(real_body), LocalDynamics(real_body),
              bounding_bounds_(periodic_box.bounding_bounds_), axis_(periodic_box.axis_),
              periodic_translation_(periodic_box.periodic_translation_),
              bound_cells_data_(bound_cells_data),
              pos_(base_particles_.pos_){};
        virtual ~PeriodicBounding(){};

        virtual void exec(Real dt = 0.0) override
        {
            setupDynamics(dt);

            particle_for(ExecutionPolicy(), bound_cells_data_[0].first,
                         [&](size_t i)
                         { checkLowerBound(i, dt); });

            particle_for(ExecutionPolicy(), bound_cells_data_[1].first,
                         [&](size_t i)
                         { checkUpperBound(i, dt); });
        };
    };
};

/**
 * @class PeriodicConditionUsingCellLinkedList
 * @brief The method imposing periodic boundary condition in an axis direction.
 *	It includes two different steps, i.e. imposing periodic bounding and condition.
 *	The first step is carried out before update cell linked list and
 *	the second after the updating.
 *	If the exec or parallel_exec is called directly, error message will be given.
 */
class PeriodicConditionUsingCellLinkedList : public BasePeriodicCondition<execution::ParallelPolicy>
{
  protected:
    /**
     * @class PeriodicCellLinkedList
     * @brief Periodic boundary condition in an axis direction
     */
    class PeriodicCellLinkedList : public PeriodicBounding
    {
      protected:
        std::mutex mutex_cell_list_entry_; /**< mutex exclusion for memory conflict */
        StdVec<CellLists> &bound_cells_data_;
        BaseCellLinkedList &cell_linked_list_;
        Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */

        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0);
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0);

      public:
        PeriodicCellLinkedList(StdVec<CellLists> &bound_cells_data,
                               RealBody &real_body, PeriodicAlongAxis &periodic_box);
        virtual ~PeriodicCellLinkedList(){};

        virtual void exec(Real dt = 0.0) override;
    };

  public:
    PeriodicBounding bounding_;
    PeriodicCellLinkedList update_cell_linked_list_;

    PeriodicConditionUsingCellLinkedList(RealBody &real_body, PeriodicAlongAxis &periodic_box)
        : BasePeriodicCondition<execution::ParallelPolicy>(real_body, periodic_box),
          bounding_(bound_cells_data_, real_body, periodic_box),
          update_cell_linked_list_(bound_cells_data_, real_body, periodic_box){};
    virtual ~PeriodicConditionUsingCellLinkedList(){};
};
} // namespace SPH
#endif // DOMAIN_BOUNDING_H
