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
 * @file relation_ck.h
 * @brief Base classes on body and particle topology relations.
 * @author Xiangyu Hu
 */

#ifndef RELATION_CK_H
#define RELATION_CK_H

#include "base_body.h"
#include "base_particles.h"
#include "implementation.h"
#include "neighbor_method.hpp"

namespace SPH
{
enum class ConfigType
{
    Eulerian,
    Lagrangian,
};

template <typename...>
class Relation;

class RelationBase
{
  public:
    virtual ~RelationBase() {};
};

template <typename SourceIdentifier, typename TargetIdentifier>
class Relation<SourceIdentifier, TargetIdentifier> : public RelationBase
{

    using SourceAdaptation = typename SourceIdentifier::BaseAdaptation;
    using TargetAdaptation = typename TargetIdentifier::BaseAdaptation;
    SharedPtrsKeeper<Entity> relation_variable_ptrs_;
    SharedPtrsKeeper<Neighbor<Base>> neighborhood_ptrs_;
    DiscreteVariable<Vecd> *assignConfigPosition(BaseParticles &particles, ConfigType config_type);

    template <class DataType>
    DiscreteVariable<DataType> *addRelationVariable(const std::string &name, size_t data_size);

  public:
    typedef SourceIdentifier SourceType;
    typedef TargetIdentifier TargetType;
    using NeighborhoodType = Neighbor<SourceAdaptation, TargetAdaptation>;
    Relation(SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers,
             ConfigType config_type = ConfigType::Eulerian);
    virtual ~Relation() {};
    SPHBody &getSPHBody() { return *sph_body_; };
    DiscreteVariable<Vecd> *dvSourcePosition() { return dv_source_pos_; };
    DiscreteVariable<UnsignedInt> *dvNeighborSize() { return dv_neighbor_size_; };
    DiscreteVariable<Vecd> *dvTargetPosition(UnsignedInt target_index = 0) { return dv_target_pos_[target_index]; };
    DiscreteVariable<UnsignedInt> *dvNeighborIndex(UnsignedInt target_index = 0) { return dv_target_neighbor_index_[target_index]; };
    DiscreteVariable<UnsignedInt> *dvParticleOffset(UnsignedInt target_index = 0) { return dv_target_particle_offset_[target_index]; };
    NeighborhoodType &getNeighborhood(UnsignedInt target_index = 0) { return *neighborhoods_[target_index]; }
    void registerComputingKernel(execution::Implementation<Base> *implementation, UnsignedInt target_index = 0);
    void resetComputingKernelUpdated(UnsignedInt target_index = 0);

    class NeighborList
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborList(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                     UnsignedInt target_index = 0);
        inline bool hasNeighor(UnsignedInt i) { return particle_offset_[i] != particle_offset_[i + 1]; };
        
      protected:
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
        inline UnsignedInt FirstNeighbor(UnsignedInt i) { return particle_offset_[i]; };
        inline UnsignedInt LastNeighbor(UnsignedInt i) { return particle_offset_[i + 1]; };
    };

  protected:
    SPHBody *sph_body_;
    BaseParticles *particles_;
    DiscreteVariable<Vecd> *dv_source_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_size_;
    StdVec<DiscreteVariable<Vecd> *> dv_target_pos_;
    UnsignedInt offset_list_size_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_target_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_target_particle_offset_;
    StdVec<NeighborhoodType *> neighborhoods_;
    StdVec<StdVec<execution::Implementation<Base> *>> registered_computing_kernels_;
};

template <typename DynamicsIdentifier>
class Inner<Relation<DynamicsIdentifier>> : public Relation<DynamicsIdentifier, DynamicsIdentifier>
{
  public:
    typedef DynamicsIdentifier SourceType;
    template <typename... Args>
    explicit Inner(DynamicsIdentifier &identifier, Args &&...args);
    virtual ~Inner() {};
    DynamicsIdentifier &getDynamicsIdentifier() { return *identifier_; };

  protected:
    DynamicsIdentifier *identifier_;
};

template <>
class Inner<> : public Inner<Relation<RealBody>>
{
  public:
    template <typename... Args>
    Inner(RealBody &real_body, Args &&...args)
        : Inner<Relation<RealBody>>(real_body, std::forward<Args>(args)...) {}
    virtual ~Inner() {};
};

template <typename SourceIdentifier, class TargetIdentifier>
class Contact<Relation<SourceIdentifier, TargetIdentifier>> : public Relation<SourceIdentifier, TargetIdentifier>
{
  public:
    Contact(SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> target_identifiers,
            ConfigType config_type = ConfigType::Eulerian);
    virtual ~Contact() {};
    SourceIdentifier &getSourceIdentifier() { return *source_identifier_; };
    StdVec<SPHBody *> getContactBodies() { return contact_bodies_; };
    StdVec<BaseParticles *> getContactParticles() { return contact_particles_; };
    StdVec<SPHAdaptation *> getContactAdaptations() { return contact_adaptations_; };
    StdVec<TargetIdentifier *> getContactIdentifiers() { return contact_identifiers_; };
    SPHBody &getContactBody(UnsignedInt target_index) { return *contact_bodies_[target_index]; };
    BaseParticles &getContactParticles(UnsignedInt target_index) { return *contact_particles_[target_index]; };
    SPHAdaptation &getContactAdaptation(UnsignedInt target_index) { return *contact_adaptations_[target_index]; };
    TargetIdentifier &getContactIdentifier(UnsignedInt target_index) { return *contact_identifiers_[target_index]; };

  protected:
    SourceIdentifier *source_identifier_;
    StdVec<SPHBody *> contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    StdVec<TargetIdentifier *> contact_identifiers_;
};

template <>
class Contact<> : public Contact<Relation<SPHBody, RealBody>>
{
  public:
    template <typename... Args>
    Contact(SPHBody &sph_body, StdVec<RealBody *> contact_bodies, Args &&...args)
        : Contact<Relation<SPHBody, RealBody>>(sph_body, contact_bodies, std::forward<Args>(args)...) {}
    virtual ~Contact() {};
};
} // namespace SPH
#endif // RELATION_CK_H
