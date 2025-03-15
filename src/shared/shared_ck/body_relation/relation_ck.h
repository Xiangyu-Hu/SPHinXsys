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
#include "implementation.h"

namespace SPH
{

template <typename...>
class Relation;

template <>
class Relation<Base>
{
    UniquePtrsKeeper<Entity> relation_variable_ptrs_;

  public:
    template <class SourceIdentifier, class TargetIdentifier>
    Relation(SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> target_identifiers);
    virtual ~Relation() {};
    SPHBody &getSPHBody() { return sph_body_; };
    DiscreteVariable<UnsignedInt> *getNeighborIndex(UnsignedInt target_index = 0);
    DiscreteVariable<UnsignedInt> *getParticleOffset(UnsignedInt target_index = 0);
    void registerComputingKernel(execution::Implementation<Base> *implementation, UnsignedInt target_index = 0);
    void resetComputingKernelUpdated(UnsignedInt target_index = 0);

    class NeighborList
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborList(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                     UnsignedInt target_index = 0);

      protected:
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
        inline UnsignedInt FirstNeighbor(UnsignedInt i) { return particle_offset_[i]; };
        inline UnsignedInt LastNeighbor(UnsignedInt i) { return particle_offset_[i + 1]; };
    };

  protected:
    SPHBody &sph_body_;
    BaseParticles &particles_;
    UnsignedInt offset_list_size_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_target_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_target_particle_offset_;
    StdVec<StdVec<execution::Implementation<Base> *>> registered_computing_kernels_;

    template <class DataType>
    DiscreteVariable<DataType> *addRelationVariable(const std::string &name, size_t data_size);
};

template <>
class Relation<Inner<>> : public Relation<Base>
{
  public:
    explicit Relation(RealBody &real_body);
    virtual ~Relation() {};
    RealBody &getRealBody() { return *real_body_; };

  protected:
    RealBody *real_body_;
};

template <class SourceIdentifier, class TargetIdentifier>
class Relation<Contact<SourceIdentifier, TargetIdentifier>> : public Relation<Base>
{
  protected:
    SourceIdentifier &source_identifier_;
    StdVec<TargetIdentifier *> contact_identifiers_;
    StdVec<RealBody *> contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;

  public:
    typedef SourceIdentifier SourceType;
    typedef TargetIdentifier TargetType;

    Relation(SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers);
    virtual ~Relation() {};
    SourceIdentifier &getSourceIdentifier() { return source_identifier_; };
    StdVec<TargetIdentifier *> getContactIdentifiers() { return contact_identifiers_; };
    TargetIdentifier &getContactIdentifier(UnsignedInt target_index) { return *contact_identifiers_[target_index]; };
    StdVec<RealBody *> getContactBodies() { return contact_bodies_; };
    StdVec<BaseParticles *> getContactParticles() { return contact_particles_; };
    StdVec<SPHAdaptation *> getContactAdaptations() { return contact_adaptations_; };
};
template <>
class Relation<Contact<>> : public Relation<Contact<SPHBody, RealBody>>
{
  public:
    Relation(SPHBody &sph_body, StdVec<RealBody *> contact_bodies)
        : Relation<Contact<SPHBody, RealBody>>(sph_body, contact_bodies) {}
    virtual ~Relation() {};
};
} // namespace SPH
#endif // RELATION_CK_H