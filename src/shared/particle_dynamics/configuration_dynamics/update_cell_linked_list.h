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
 * @file    update_cell_linked_list.h
 * @brief   Collection of dynamics for particle configuration.
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_CELL_LINKED_LIST_H
#define UPDATE_CELL_LINKED_LIST_H

#include "base_local_dynamics.h"

namespace SPH
{
template <typename... T>
class UpdateCellLinkedList;

template <typename... T>
struct AtomicUnsignedIntRef;

template <class ExecutionPolicy>
struct AtomicUnsignedIntRef<ExecutionPolicy>
{
};

template <typename CellLinkedListType>
class UpdateCellLinkedList<CellLinkedListType> : public LocalDynamics
{
  public:
    template <class ExecutionPolicy>
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};

    UnsignedInt *setParticleOffsetListUpperBound();
    void setParticleOffsetListUpperBound(const ExecutionPolicy &execution_policy);
    void setParticleOffsetListUpperBound(const ParallelDevicePolicy &par_device);

    void exclusiveScanParticleOffsetList(const ExecutionPolicy &execution_policy);
    void exclusiveScanParticleOffsetList(const ParallelDevicePolicy &par_device);
    class ComputingKernel
    {
      public:
        ComputingKernel(UpdateCellLinkedList<CellLinkedListType> &update_cell_linked_list);

        void clearOffsetLists(UnsignedInt linear_cell_index);
        void incrementCellSize(UnsignedInt particle_i);
        void updateCellLists(UnsignedInt particle_i);

      protected:
        friend class UpdateCellLinkedList<CellLinkedListType>;
        Mesh mesh_;
        UnsignedInt number_of_cells_, particles_bound_;

        Vecd *pos_;
        UnsignedInt *particle_id_list_;
        UnsignedInt *particle_offset_list_;
        UnsignedInt *current_size_list_;
        UnsignedInt *total_real_particles_;
    };

  protected:
    Mesh mesh_;
    UnsignedInt number_of_cells_, particles_bound_;
    SingularVariable<UnsignedInt> *v_total_real_particles_;
    DiscreteVariable<UnsignedInt> *v_particle_offset_list_;

    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;
    UnsignedInt *current_size_list_;
};

template <class CellLinkedListType, class ExecutionPolicy = ParallelPolicy>
class UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>
    : public UpdateCellLinkedList<CellLinkedListType>, public BaseDynamics<void>
{
    using LocalDynamicsType = typename UpdateCellLinkedList<CellLinkedListType>;
    using ComputingKernel = typename LocalDynamicsType::ComputingKernel;

  public:
    UpdateCellLinkedList(RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    Implementation<LocalDynamicsType, ExecutionPolicy> kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_H
