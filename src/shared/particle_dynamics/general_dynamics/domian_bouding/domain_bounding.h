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
  public:
    PeriodicAlongAxis(BoundingBoxd bounding_bounds, int axis)
        : bounding_bounds_(bounding_bounds), axis_(axis),
          periodic_translation_(Vecd::Zero())
    {
        periodic_translation_[axis] =
            bounding_bounds.upper_[axis] - bounding_bounds.lower_[axis];
    };
    virtual ~PeriodicAlongAxis() {};
    BoundingBoxd getBoundingBox() { return bounding_bounds_; };
    int getAxis() { return axis_; };
    Vecd getPeriodicTranslation() { return periodic_translation_; };

  protected:
    BoundingBoxd bounding_bounds_; /**< lower and upper bound for checking. */
    const int axis_;              /**< the axis directions for bounding*/
    Vecd periodic_translation_;
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
        cell_linked_list.tagBoundingCells(bound_cells_data_, periodic_box.getBoundingBox(),
                                          periodic_box.getAxis());
    };
    virtual ~BasePeriodicCondition() {};

  protected:
    StdVec<CellLists> bound_cells_data_;

    /**
     * @class PeriodicBounding
     * @brief Periodic bounding particle position in an axis direction
     */
    class PeriodicBounding : public LocalDynamics, public BaseDynamics<void>
    {
      protected:
        BoundingBoxd bounding_bounds_;
        const int axis_;
        Vecd periodic_translation_;
        Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */
        StdVec<CellLists> &bound_cells_data_;
        Vecd *pos_;

        virtual void checkLowerBound(size_t index_i, Real dt = 0.0)
        {
            if (pos_[index_i][axis_] < bounding_bounds_.lower_[axis_])
                pos_[index_i][axis_] += periodic_translation_[axis_];
        };

        virtual void checkUpperBound(size_t index_i, Real dt = 0.0)
        {
            if (pos_[index_i][axis_] > bounding_bounds_.upper_[axis_])
                pos_[index_i][axis_] -= periodic_translation_[axis_];
        };

      public:
        PeriodicBounding(StdVec<CellLists> &bound_cells_data,
                         RealBody &real_body, PeriodicAlongAxis &periodic_box)
            : LocalDynamics(real_body), BaseDynamics<void>(),
              bounding_bounds_(periodic_box.getBoundingBox()), axis_(periodic_box.getAxis()),
              periodic_translation_(periodic_box.getPeriodicTranslation()),
              cut_off_radius_max_(real_body.getSPHAdaptation().getKernel()->CutOffRadius()),
              bound_cells_data_(bound_cells_data),
              pos_(particles_->getVariableDataByName<Vecd>("Position")) {};
        virtual ~PeriodicBounding() {};

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
    class PeriodicCellLinkedList : public PeriodicBounding
    {
      protected:
        std::mutex mutex_cell_list_entry_; /**< mutex exclusion for memory conflict */
        BaseCellLinkedList &cell_linked_list_;

        void InsertListDataNearLowerBound(ListDataVector &cell_list_data, Real dt = 0.0);
        void InsertListDataNearUpperBound(ListDataVector &cell_list_data, Real dt = 0.0);

      public:
        PeriodicCellLinkedList(StdVec<CellLists> &bound_cells_data,
                               RealBody &real_body, PeriodicAlongAxis &periodic_box);
        virtual ~PeriodicCellLinkedList() {};

        virtual void exec(Real dt = 0.0) override;
    };

  public:
    PeriodicBounding bounding_;
    PeriodicCellLinkedList update_cell_linked_list_;

    PeriodicConditionUsingCellLinkedList(RealBody &real_body, PeriodicAlongAxis &periodic_box);
    virtual ~PeriodicConditionUsingCellLinkedList() {};
};
} // namespace SPH
#endif // DOMAIN_BOUNDING_H
