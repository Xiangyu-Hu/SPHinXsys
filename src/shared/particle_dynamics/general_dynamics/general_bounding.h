/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	general_bounding.h
 * @brief 	This is the particle dynamics for domain bounding
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_BOUNDING_H
#define GENERAL_BOUNDING_H

#include "general_dynamics.h"

namespace SPH
{

/**
 * @class BoundingAlongAxis
 * @brief Bounding particle position in along axis.
 * The axis must be 0, 1 for 2d and 0, 1, 2 for 3d
 */
class BoundingAlongAxis : public BaseDynamics<void>,
                          public LocalDynamics,
                          public GeneralDataDelegateSimple
{
  protected:
    const int axis_;              /**< the axis directions for bounding*/
    BoundingBox bounding_bounds_; /**< lower and upper bound for checking. */
    StdLargeVec<Vecd> &pos_;
    BaseCellLinkedList &cell_linked_list_;
    Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */

  public:
    BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis);
    virtual ~BoundingAlongAxis(){};
};

class PeriodicBoundary
{
  public:
    PeriodicBoundary(BoundingBox bounding_bounds, int axis)
        : bounding_bounds_(bounding_bounds), axis_(axis),
          translation_(bounding_bounds.second_[axis] - bounding_bounds.first_[axis]){};

    Vecd LowerBoundTranslation(Vecd pos)
    {
        pos[axis_] += translation_;
        return pos;
    };

    Vecd UpperBoundTranslation(Vecd pos)
    {
        pos[axis_] -= translation_;
        return pos;
    };

  protected:
    BoundingBox bounding_bounds_;
    int axis_;
    Real translation_;
};

class LeeEdwardsBoundary : public PeriodicBoundary
{
  public:
    LeeEdwardsBoundary(BoundingBox bounding_bounds, int axis, int shear_direction, Real shear_rate)
        : PeriodicBoundary(bounding_bounds, axis), shear_direction_(shear_direction),
          shear_center_(0.5 * (bounding_bounds.second_[shear_direction] + bounding_bounds.first_[shear_direction])),
          pos_increment_(bounding_bounds.second_[shear_direction] - bounding_bounds.first_[shear_direction]),
          vel_increment_(pos_increment_ * shear_rate){};

    Vecd LowerBoundTranslation(Vecd pos)
    {
        pos[axis_] += translation_;
        flipPosition(pos);
        return pos;
    };

    Vecd UpperBoundTranslation(Vecd pos)
    {
        pos[axis_] -= translation_;
        flipPosition(pos);
        return pos;
    };

    Vecd flipVelocity(const Vecd &pos, Vecd vel)
    {
        bool is_positive_side = pos[shear_direction_] > shear_center_;
        vel[shear_direction_] += is_positive_side ? -vel_increment_ : vel_increment_;
        return vel;
    };

  private:
    int shear_direction_; // upper bound moving toward positive side
    Real shear_center_;
    Real pos_increment_;
    Real vel_increment_;

    void flipPosition(Vecd &pos)
    {
        bool is_positive_side = pos[shear_direction_] > shear_center_;
        pos[shear_direction_] += is_positive_side ? -pos_increment_ : pos_increment_;
    };
};

/**
 * @class BasePeriodicCondition
 * @brief Base class for two different type periodic boundary conditions.
 */
template <class ExecutionPolicy>
class BasePeriodicCondition
{
  protected:
    Vecd periodic_translation_;
    StdVec<CellLists> bound_cells_data_;
    Vecd setPeriodicTranslation(BoundingBox &bounding_bounds, int axis)
    {
        Vecd periodic_translation = Vecd::Zero();
        periodic_translation[axis] =
            bounding_bounds.second_[axis] - bounding_bounds.first_[axis];
        return periodic_translation;
    };

    /**
     * @class PeriodicBounding
     * @brief Periodic bounding particle position in an axis direction
     */
    class PeriodicBounding : public BoundingAlongAxis
    {
      protected:
        Vecd &periodic_translation_;
        StdVec<CellLists> &bound_cells_data_;

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
        PeriodicBounding(Vecd &periodic_translation,
                         StdVec<CellLists> &bound_cells_data,
                         RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : BoundingAlongAxis(real_body, bounding_bounds, axis),
              periodic_translation_(periodic_translation),
              bound_cells_data_(bound_cells_data){};
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

    /**
     * @class PeriodicImage
     * @brief Creating periodic image for boundary condition in an axis direction
     */
    class PeriodicImage : public BoundingAlongAxis
    {
      protected:
        std::mutex mutex_cell_list_entry_; /**< mutex exclusion for memory conflict */
        Vecd &periodic_translation_;
        StdVec<CellLists> &bound_cells_data_;
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) = 0;
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) = 0;

      public:
        PeriodicImage(Vecd &periodic_translation,
                      StdVec<CellLists> &bound_cells_data,
                      RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : BoundingAlongAxis(real_body, bounding_bounds, axis),
              periodic_translation_(periodic_translation),
              bound_cells_data_(bound_cells_data){};
        ;
        virtual ~PeriodicImage(){};

        virtual void exec(Real dt = 0.0) override
        {
            setupDynamics(dt);

            particle_for(ExecutionPolicy(), bound_cells_data_[0].second,
                         [&](ListDataVector *cell_list)
                         { checkLowerBound(*cell_list, dt); });

            particle_for(ExecutionPolicy(), bound_cells_data_[1].second,
                         [&](ListDataVector *cell_list)
                         { checkUpperBound(*cell_list, dt); });
        };
    };

  public:
    BasePeriodicCondition(RealBody &real_body, BoundingBox bounding_bounds, int axis)
        : periodic_translation_(setPeriodicTranslation(bounding_bounds, axis))
    {
        bound_cells_data_.resize(2);
        BaseCellLinkedList &cell_linked_list = real_body.getCellLinkedList();
        cell_linked_list.tagBoundingCells(bound_cells_data_, bounding_bounds, axis);
        if (periodic_translation_.norm() < real_body.sph_adaptation_->ReferenceSpacing())
        {
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            std::cout << "\n Periodic bounding failure: bounds not defined!" << std::endl;
            exit(1);
        }
    };
    virtual ~BasePeriodicCondition(){};
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
     * @brief Periodic boundary condition in an axis direction using cell linked list image
     */
    class PeriodicCellLinkedList : public PeriodicImage
    {
      protected:
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) override;
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) override;

      public:
        PeriodicCellLinkedList(Vecd &periodic_translation,
                               StdVec<CellLists> &bound_cells_data,
                               RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicImage(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis){};
        virtual ~PeriodicCellLinkedList(){};
    };

  public:
    PeriodicConditionUsingCellLinkedList(RealBody &real_body, BoundingBox bounding_bounds, int axis)
        : BasePeriodicCondition<execution::ParallelPolicy>(real_body, bounding_bounds, axis),
          bounding_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis),
          update_cell_linked_list_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis){};
    virtual ~PeriodicConditionUsingCellLinkedList(){};

    PeriodicBounding bounding_;
    PeriodicCellLinkedList update_cell_linked_list_;
};

/**
 * @class PeriodicConditionUsingGhostParticles
 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
 *	The first step is carried out before update cell linked list and
 *	the second and third after the updating.
 *	If the exec or parallel_exec is called directly, error message will be given.
 *  Note that, currently, this class is not for periodic condition in combined directions,
 *  such as periodic condition in both x and y directions.
 */
class PeriodicConditionUsingGhostParticles : public BasePeriodicCondition<execution::ParallelPolicy>
{
  protected:
    StdVec<IndexVector> ghost_particles_;

    /**
     * @class PeriodicGhost
     * @brief create ghost particles in an axis direction
     */
    class PeriodicGhost : public PeriodicImage
    {
      protected:
        StdVec<IndexVector> &ghost_particles_;
        virtual void setupDynamics(Real dt = 0.0) override;
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) override;
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) override;

      public:
        PeriodicGhost(Vecd &periodic_translation,
                      StdVec<CellLists> &bound_cells_data,
                      StdVec<IndexVector> &ghost_particles,
                      RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicImage(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
              ghost_particles_(ghost_particles){};
        virtual ~PeriodicGhost(){};
    };

    /**
     * @class UpdatePeriodicGhost
     * @brief update ghost particles in an axis direction
     */
    class UpdatePeriodicGhost : public PeriodicBounding
    {
      protected:
        StdVec<IndexVector> &ghost_particles_;
        void checkLowerBound(size_t index_i, Real dt = 0.0) override;
        void checkUpperBound(size_t index_i, Real dt = 0.0) override;

      public:
        UpdatePeriodicGhost(Vecd &periodic_translation,
                            StdVec<CellLists> &bound_cells_data,
                            StdVec<IndexVector> &ghost_particles,
                            RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicBounding(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
              ghost_particles_(ghost_particles){};
        virtual ~UpdatePeriodicGhost(){};

        virtual void exec(Real dt = 0.0) override;
    };

  public:
    PeriodicConditionUsingGhostParticles(RealBody &real_body, BoundingBox bounding_bounds, int axis)
        : BasePeriodicCondition<execution::ParallelPolicy>(real_body, bounding_bounds, axis),
          bounding_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis),
          ghost_creation_(periodic_translation_, bound_cells_data_, ghost_particles_, real_body, bounding_bounds, axis),
          ghost_update_(periodic_translation_, bound_cells_data_, ghost_particles_, real_body, bounding_bounds, axis)
    {
        ghost_particles_.resize(2);
    };

    virtual ~PeriodicConditionUsingGhostParticles(){};

    PeriodicBounding bounding_;
    PeriodicGhost ghost_creation_;
    UpdatePeriodicGhost ghost_update_;
};
} // namespace SPH
#endif // GENERAL_BOUNDING_H
