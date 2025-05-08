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
 * @file base_solver.h
 * @brief Collection of particle dynamics.
 * @author	Xiangyu Hu
 */

#ifndef BASE_SOLVER_H
#define BASE_SOLVER_H

#include "base_particle_dynamics.h"
#include "simple_algorithms_ck.h"
#include "update_cell_linked_list.h"
#include "interaction_algorithms_ck.h"

namespace SPH
{
class SPHSolver
{
    UniquePtrsKeeper<BaseDynamics<void>> particle_dynamics_keeper_;

  public:
    SPHSolver() {};
    virtual ~SPHSolver() {};

    template <typename ExecutePolicy, class DynamicsIdentifier, typename... Args>
    BaseDynamics<void> &addCellLinkedListDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        BaseDynamics<void> *cell_linked_List_dynamics =
            particle_dynamics_keeper_.createPtr<
                UpdateCellLinkedList<ExecutePolicy, DynamicsIdentifier>>(
                identifier, std::forward<Args>(args)...);
        return *cell_linked_List_dynamics;
    };


        template <typename ExecutePolicy, class DynamicsIdentifier, typename... Args>
    BaseDynamics<void> &addInteractionDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        BaseDynamics<void> *interaction_dynamics =
            particle_dynamics_keeper_.createPtr<
                UpdateCellLinkedList<ExecutePolicy, DynamicsIdentifier>>(
                identifier, std::forward<Args>(args)...);
        return *cell_linked_List_dynamics;
    };
};
} // namespace SPH
#endif // BASE_SOLVER_H