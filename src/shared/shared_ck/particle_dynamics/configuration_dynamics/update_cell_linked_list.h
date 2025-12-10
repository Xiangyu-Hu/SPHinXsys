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
 * @file    update_cell_linked_list.h
 * @brief   Collection of dynamics for particle configuration.
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_CELL_LINKED_LIST_H
#define UPDATE_CELL_LINKED_LIST_H

#include "base_configuration_dynamics.h"

#include "all_bodies.h"
#include "all_particle_dynamics.h"

namespace SPH
{
template <typename... T>
class UpdateCellLinkedList;

template <class ExecutionPolicy, typename DynamicsIdentifier>
class UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>
    : public BaseLocalDynamics<DynamicsIdentifier>, public BaseDynamics<void>
{
    typedef UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier> EncloserType;
    using ParticleMask = typename DynamicsIdentifier::ListedParticleMask;

  protected:
    CellLinkedList &cell_linked_list_;
    Mesh mesh_;
    UnsignedInt total_number_of_cells_;
    UnsignedInt cell_offset_list_size_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_particle_index_;
    DiscreteVariable<UnsignedInt> *dv_cell_offset_;
    DiscreteVariable<UnsignedInt> *cell_dv_current_list_size_;

  public:
    UpdateCellLinkedList(DynamicsIdentifier &identifier);
    virtual ~UpdateCellLinkedList() {};

    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void clearAllLists(UnsignedInt index_i);
        void incrementCellSize(UnsignedInt index_i);
        void updateCellList(UnsignedInt index_i);

      protected:
        Mesh mesh_;
        ParticleMask particle_mask_;

        Vecd *pos_;
        UnsignedInt *particle_index_;
        UnsignedInt *cell_offset_;
        UnsignedInt *current_list_size_;
    };

    virtual void exec(Real dt = 0.0) override;

  protected:
    DynamicsIdentifier &identifier_;
    Implementation<ExecutionPolicy, EncloserType, ComputingKernel> kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_H
