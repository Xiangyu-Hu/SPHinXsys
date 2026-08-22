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
#include "base_implementation.h"
#include "base_particles.h"
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
    RelationBase(const std::string &src_name, const std::string &tgt_name)
        : name_(src_name == tgt_name ? src_name : src_name + tgt_name),
          reverse_name_(tgt_name == src_name ? tgt_name : tgt_name + src_name) {};
    virtual ~RelationBase() {};
    const std::string &Name() const { return name_; }
    const std::string &ReverseName() const { return reverse_name_; }

  protected:
    std::string name_;
    std::string reverse_name_;
};

template <typename SourceIdentifier, typename TargetIdentifier>
class Relation<SourceIdentifier, TargetIdentifier> : public RelationBase
{

    using SourceAdaptation = typename SourceIdentifier::Adaptation;
    using TargetAdaptation = typename TargetIdentifier::Adaptation;
    UniquePtrsKeeper<Quantity> relation_variable_ptrs_;
    UniquePtrKeeper<Neighbor<Base>> neighborhood_ptr_;
    DiscreteVariable<Vecd> *assignConfigPosition(BaseParticles &particles, ConfigType config_type);

    template <class DataType>
    DiscreteVariable<DataType> *addRelationVariable(const std::string &name, size_t data_size);

  public:
    typedef SourceIdentifier SourceType;
    typedef TargetIdentifier TargetType;
    using NeighborhoodType = Neighbor<SourceAdaptation, TargetAdaptation>;
    Relation(SourceIdentifier &src_identifier, TargetIdentifier &tgt_identifier,
             ConfigType config_type = ConfigType::Eulerian);
    virtual ~Relation() {};
    SPHBody &getSPHBody() { return *sph_body_; };
    DiscreteVariable<Vecd> *dvSourcePosition() { return dv_source_pos_; };
    DiscreteVariable<UnsignedInt> *dvNeighborSize() { return dv_neighbor_size_; };
    DiscreteVariable<Vecd> *dvTargetPosition() { return dv_target_pos_; };
    DiscreteVariable<UnsignedInt> *dvNeighborIndex() { return dv_target_neighbor_index_; };
    DiscreteVariable<UnsignedInt> *dvParticleOffset() { return dv_target_particle_offset_; };
    NeighborhoodType &getNeighborhood() { return *neighborhood_; }
    void registerComputingKernel(execution::Implementation<Base> *implementation);
    void resetComputingKernelUpdated();

    class NeighborList
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborList(const ExecutionPolicy &ex_policy, EncloserType &encloser);

      protected:
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
        inline UnsignedInt FirstNeighbor(UnsignedInt i) { return particle_offset_[i]; };
        inline UnsignedInt LastNeighbor(UnsignedInt i) { return particle_offset_[i + 1]; };
        inline bool hasNeighbor(UnsignedInt i) { return particle_offset_[i] != particle_offset_[i + 1]; };
    };

  protected:
    SPHBody *sph_body_;
    BaseParticles *particles_;
    DiscreteVariable<Vecd> *dv_source_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_size_;
    DiscreteVariable<Vecd> *dv_target_pos_;
    UnsignedInt offset_list_size_;
    DiscreteVariable<UnsignedInt> *dv_target_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_target_particle_offset_;
    NeighborhoodType *neighborhood_;
    StdVec<execution::Implementation<Base> *> registered_computing_kernels_;
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
    Contact(SourceIdentifier &src_identifier, TargetIdentifier &tgt_identifier,
            ConfigType config_type = ConfigType::Eulerian);
    virtual ~Contact() {};
    SourceIdentifier &getSourceIdentifier() { return *src_identifier_; };
    SPHBody &getContactBody() { return *contact_body_; };
    BaseParticles &getContactParticles() { return *contact_particles_; };
    SPHAdaptation &getContactAdaptation() { return *contact_adaptation_; };
    TargetIdentifier &getTargetIdentifier() { return *tgt_identifier_; };

  protected:
    SourceIdentifier *src_identifier_;
    TargetIdentifier *tgt_identifier_;
    SPHBody *contact_body_;
    BaseParticles *contact_particles_;
    SPHAdaptation *contact_adaptation_;
};

template <>
class Contact<> : public Contact<Relation<SPHBody, RealBody>>
{
  public:
    template <typename... Args>
    Contact(SPHBody &sph_body, RealBody &contact_body, Args &&...args)
        : Contact<Relation<SPHBody, RealBody>>(sph_body, contact_body, std::forward<Args>(args)...) {}
    virtual ~Contact() {};
};
} // namespace SPH
#endif // RELATION_CK_H
