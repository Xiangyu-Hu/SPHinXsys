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

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.hpp"

#include "tbb/parallel_for.h"

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

template <class MeshType>
class ParticleCellLinkedList : public MeshType
{
  protected:
    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;

  public:
    ParticleCellLinkedList(
        const MeshType &mesh, Vecd *pos_,
        UnsignedInt *particle_id_list, UnsignedInt *particle_offset_list);
    ~ParticleCellLinkedList(){};

    template <typename NeighborhoodType, typename FunctionOnEach, typename... Args>
    void forEachNeighbor(UnsignedInt index_i, const Vecd *source_pos,
                         const NeighborhoodType &neighborhood,
                         const FunctionOnEach &function, Args &&...) const;
};

template <typename MeshType>
class UpdateCellLinkedList<MeshType> : public LocalDynamics
{
  public:
    template <class ExecutionPolicy>
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};

    template <class ExecutionPolicy>
    void clearParticleOffsetList(const ExecutionPolicy &execution_policy);
    void clearParticleOffsetList(const ParallelDevicePolicy &par_device);

    template <class ExecutionPolicy>
    void incrementCellSize(const ExecutionPolicy &execution_policy);
    void incrementCellSize(const ParallelDevicePolicy &par_device);

    template <class ExecutionPolicy>
    void exclusiveScanParticleOffsetList(const ExecutionPolicy &execution_policy);
    void exclusiveScanParticleOffsetList(const ParallelDevicePolicy &par_device);

    template <class ExecutionPolicy>
    void updateCellLists(const ExecutionPolicy &execution_policy);
    void updateCellLists(const ParallelDevicePolicy &par_device);

    ParticleCellLinkedList<MeshType> getParticleCellLinkedList() const;

  protected:
    const Mesh *mesh_;
    UnsignedInt number_of_cells_, particles_bound_;

    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;
    UnsignedInt *current_size_list_;
};

template <class MeshType, class ExecutionPolicy>
class UpdateCellLinkedList<MeshType, ExecutionPolicy>
    : public UpdateCellLinkedList<MeshType>, public BaseDynamics<void>
{
  public:
    UpdateCellLinkedList(RealBody &real_body);
    virtual ~UpdateCellLinkedList(){};
    virtual void exec(Real dt = 0.0) override;
};

template <typename T, typename Op>
exclusive_scan(T *first, T *last, T *d_first, Op op)
{
    // Exclusive scan is the same as inclusive, but shifted by one
    UnsignedInt scan_size = last - first;
    using range_type = tbb::blocked_range<UnsignedInt>;
    tbb::parallel_scan(
        range_type(0, scan_size), ZeroData<T>::value,
        [&](const range_type &r, T sum, bool is_final_scan) {
            T tmp = sum;
            for (UnsignedInt i = r.begin(); i < r.end(); ++i)
            {
                tmp = op(tmp, first[i]);
                if (is_final_scan)
                {
                    d_first[i + 1] = tmp;
                }
            }
            return tmp;
        },
        [&](const T &a, const T &b) {
            return op(a, b);
        });
}

} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_H
