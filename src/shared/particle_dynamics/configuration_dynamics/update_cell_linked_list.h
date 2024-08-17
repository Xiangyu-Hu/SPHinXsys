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

#include "base_configuration_dynamics.h"

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
template <typename... T>
class UpdateCellLinkedList;

template <class MeshType>
class ParticleCellLinkedList : public MeshType
{
  protected:
    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;
    Real grid_spacing_squared_;

  public:
    ParticleCellLinkedList(const MeshType &mesh, Vecd *pos_,
                           UnsignedInt *particle_id_list, UnsignedInt *particle_offset_list);
    ~ParticleCellLinkedList(){};

    template <typename FunctionOnEach>
    void forEachNeighbor(UnsignedInt index_i, const Vecd *source_pos,
                         const FunctionOnEach &function) const;
};

template <typename MeshType>
class UpdateCellLinkedList<MeshType> : public LocalDynamics
{
  public:
    template <class ExecutionPolicy>
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};

    ParticleCellLinkedList<MeshType> getParticleCellLinkedList() const;
    UnsignedInt getNumberOfCellsPlusOne() { return number_of_cells_plus_one_; }
    template <class T>
    class ComputingKernel
    {
      public:
        ComputingKernel(UpdateCellLinkedList<MeshType> &update_cell_linked_list);
        void clearAllLists(UnsignedInt index_i);
        void incrementCellSize(UnsignedInt index_i);
        void updateCellLists(UnsignedInt index_i);

      protected:
        friend class UpdateCellLinkedList<MeshType>;
        Mesh mesh_;
        UnsignedInt number_of_cells_plus_one_;
        UnsignedInt particle_id_list_size_; // at least number_of_cells_pluse_one_

        Vecd *pos_;
        UnsignedInt *particle_id_list_;
        UnsignedInt *particle_offset_list_;
        UnsignedInt *current_size_list_;
    };

  protected:
    Mesh mesh_;
    UnsignedInt number_of_cells_plus_one_;
    UnsignedInt particle_id_list_size_; // at least number_of_cells_pluse_one_

    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;
    UnsignedInt *current_size_list_;
};

template <class MeshType, class ExecutionPolicy>
class UpdateCellLinkedList<MeshType, ExecutionPolicy>
    : public UpdateCellLinkedList<MeshType>, public BaseDynamics<void>
{
    using ComputingKernel = typename UpdateCellLinkedList<MeshType>::
        template ComputingKernel<ExecutionPolicy>;

  public:
    UpdateCellLinkedList(RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    Implementation<UpdateCellLinkedList<MeshType>, ExecutionPolicy> kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_H
