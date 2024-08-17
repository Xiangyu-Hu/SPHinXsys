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
 * @file    update_body_relation.h
 * @brief   Collection of dynamics for particle configuration.
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_BODY_RELATION_H
#define UPDATE_BODY_RELATION_H

#include "update_cell_linked_list.h"

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
template <int DIMENSION>
struct NeighborSize
{
    static inline UnsignedInt value = 0;
};

template <>
struct NeighborSize<1>
{
    static inline UnsignedInt value = int(2.0 * 2.6);
};
template <>
struct NeighborSize<2>
{
    static inline UnsignedInt value = int(Pi * 2.6 * 2.6);
};
template <>
struct NeighborSize<3>
{
    static inline UnsignedInt value = int(4.0 * Pi * 2.6 * 2.6 * 2.6 / 3.0);
};

template <typename... T>
class Relation;

class ParticleNeighborList
{
  protected:
    UnsignedInt *neighbor_id_list_;
    UnsignedInt *neighbor_offset_list_;

  public:
    ParticleNeighborList(UnsignedInt *neighbor_id_list, UnsignedInt *neighbor_offset_list);
    ~ParticleNeighborList(){};

    template <typename FunctionOnEach>
    void forEachNeighbor(UnsignedInt index_i, const Vecd *source_pos,
                         const FunctionOnEach &function) const;
};

template <class ParticleCellLinkedListType>
class Relation<Inner<ParticleCellLinkedListType>> : public LocalDynamics
{

  public:
    template <class ExecutionPolicy>
    explicit Relation(const ExecutionPolicy &execution_policy, RealBody &real_body,
                      const ParticleCellLinkedListType &particle_cell_linked_list);
    virtual ~Relation(){};

    ParticleNeighborList getParticleNeighborList()
    {
        return ParticleNeighborList(neighbor_id_list_, neighbor_offset_list_);
    };

    template <class T>
    class ComputingKernel
    {
      public:
        ComputingKernel(Relation<Inner<ParticleCellLinkedListType>> &update_inner_relation);
        void clearAllLists(UnsignedInt index_i);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborIDList(UnsignedInt index_i);

      protected:
        friend class Relation<Inner<ParticleCellLinkedListType>>;
        ParticleCellLinkedListType particle_cell_linked_list_;
        UnsignedInt real_particle_bound_plus_one_;
        UnsignedInt neighbor_id_list_size_;

        Vecd *pos_;
        UnsignedInt *neighbor_id_list_;
        UnsignedInt *neighbor_offset_list_;
        UnsignedInt *neighbor_size_list_;
    };

  protected:
    ParticleCellLinkedListType particle_cell_linked_list_;
    UnsignedInt real_particle_bound_plus_one_;
    UnsignedInt neighbor_id_list_size_;

    Vecd *pos_;
    UnsignedInt *neighbor_id_list_;
    UnsignedInt *neighbor_offset_list_;
    UnsignedInt *neighbor_size_list_;
};

template <typename... T, class ExecutionPolicy>
class Relation<ExecutionPolicy, T...>
    : public Relation<T...>, public BaseDynamics<void>
{
    using ComputingKernel = typename Relation<T...>::
        template ComputingKernel<ExecutionPolicy>;

  public:
    Relation(RealBody &real_body, T &&...parameters);
    virtual ~Relation(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    Implementation<Relation<T...>, ExecutionPolicy> kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
