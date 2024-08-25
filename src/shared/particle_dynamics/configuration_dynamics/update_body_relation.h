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
template <typename... T>
class BodyRelationUpdate;

template <>
class BodyRelationUpdate<Inner<>> : public LocalDynamics
{

  public:
    explicit BodyRelationUpdate(BaseInnerRelation &inner_relation);
    virtual ~BodyRelationUpdate(){};
    template <class T>
    class ComputingKernel
    {
      public:
        ComputingKernel(BodyRelationUpdate<Inner<>> &update_inner_relation);
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
    BodyRelationType &body_relation_;
    UnsignedInt real_particle_bound_plus_one_;
    UnsignedInt neighbor_id_list_size_;

    Vecd *pos_;
    UnsignedInt *neighbor_id_list_;
    UnsignedInt *neighbor_offset_list_;
    UnsignedInt *neighbor_size_list_;
};

template <class RelationType, class ExecutionPolicy>
class UpdateRelation : public RelationType, public BaseDynamics<void>
{
    using ComputingKernel = typename RelationType::
        template ComputingKernel<ExecutionPolicy>;

  public:
    template <typename... Args>
    UpdateRelation(RealBody &real_body, Args &&...args);
    virtual ~UpdateRelation(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    Implementation<RelationType, ExecutionPolicy> kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
