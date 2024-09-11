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
 * @file relation_ck.h
 * @brief Base classes on body and particle topology relations.
 * @author Xiangyu Hu
 */

#ifndef RELATION_CK_H
#define RELATION_CK_H

#include "base_body.h"
#include "base_particles.h"
#include "execution.h"

namespace SPH
{

template <typename...>
class Relation;

template <>
class Relation<Base>
{
    UniquePtrsKeeper<BaseVariable> relation_variable_ptrs_;

  public:
    explicit Relation(SPHBody &sph_body);
    virtual ~Relation(){};
    SPHBody &getSPHBody() { return sph_body_; };
    UnsignedInt getParticleOffsetListSize() { return offset_list_size_; };

  protected:
    SPHBody &sph_body_;
    BaseParticles &particles_;
    UnsignedInt offset_list_size_;

    template <class DataType>
    DiscreteVariable<DataType> *addRelationVariable(const std::string &name, size_t data_size);
};

template <>
class Relation<Inner<>> : public Relation<Base>
{
  public:
    explicit Relation(RealBody &real_body);
    virtual ~Relation(){};
    RealBody &getRealBody() { return *real_body_; };
    CellLinkedList &getCellLinkedList() { return cell_linked_list_; };
    DiscreteVariable<UnsignedInt> *getNeighborIndex() { return dv_neighbor_index_; };
    DiscreteVariable<UnsignedInt> *getParticleOffset() { return dv_particle_offset_; };
    void registerComputingKernel(execution::Implementation<Base> *implementation);
    void resetComputingKernelUpdated();

  protected:
    RealBody *real_body_;
    CellLinkedList &cell_linked_list_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
    StdVec<execution::Implementation<Base> *> all_inner_computing_kernels_;
};

template <>
class Relation<Contact<>> : public Relation<Base>
{
  protected:
    RealBodyVector contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    StdVec<CellLinkedList *> target_cell_linked_lists_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
    StdVec<StdVec<execution::Implementation<Base> *>> all_contact_computing_kernels_;

  public:
    Relation(SPHBody &sph_body, RealBodyVector contact_bodies);
    virtual ~Relation(){};
    RealBodyVector getContactBodies() { return contact_bodies_; };
    StdVec<BaseParticles *> getContactParticles() { return contact_particles_; };
    StdVec<SPHAdaptation *> getContactAdaptations() { return contact_adaptations_; };
    StdVec<CellLinkedList *> getContactCellLinkedList() { return target_cell_linked_lists_; }
    StdVec<DiscreteVariable<UnsignedInt> *> getContactNeighborIndex() { return dv_contact_neighbor_index_; };
    StdVec<DiscreteVariable<UnsignedInt> *> getContactParticleOffset() { return dv_contact_particle_offset_; };
    void registerComputingKernel(execution::Implementation<Base> *implementation, UnsignedInt contact_index);
    void resetComputingKernelUpdated(UnsignedInt contact_index);
};
} // namespace SPH
#endif // RELATION_CK_H