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
 * @file 	base_particle_dynamics.h
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLE_DYNAMICS_H
#define BASE_PARTICLE_DYNAMICS_H

#include "all_body_relations.h"
#include "base_body.h"
#include "base_data_package.h"
#include "neighborhood.h"
#include "sph_data_containers.h"

#include <functional>

using namespace std::placeholders;

namespace SPH
{
/**
 * @class GlobalStaticVariables
 * @brief A place to put all global variables
 */
class GlobalStaticVariables
{
  public:
    explicit GlobalStaticVariables(){};
    virtual ~GlobalStaticVariables(){};

    /** the physical time is global value for all dynamics */
    static inline Real physical_time_ = 0.0;
};

/**
 * @class BaseDynamics
 * @brief The base class for all dynamics
 * This class contains only the interface functions available
 * for all dynamics. An specific implementation should be realized.
 */
template <class ReturnType = void>
class BaseDynamics : public GlobalStaticVariables
{
  public:
    BaseDynamics(SPHBody &sph_body)
        : sph_body_(sph_body), is_newly_updated_(false){};
    virtual ~BaseDynamics(){};
    bool checkNewlyUpdated() { return is_newly_updated_; };
    void setNotNewlyUpdated() { is_newly_updated_ = false; };

    void setUpdated()
    {
        sph_body_.setNewlyUpdated();
        is_newly_updated_ = true;
    };

    /** There is the interface functions for computing. */
    virtual ReturnType exec(Real dt = 0.0) = 0;

  private:
    SPHBody &sph_body_;
    bool is_newly_updated_;
};

/**
 * @class DataDelegateBase
 * @brief empty base class for mixin template.
 */
class DataDelegateEmptyBase
{
  public:
    explicit DataDelegateEmptyBase(SPHBody &sph_body){};
    virtual ~DataDelegateEmptyBase(){};
};

/**
 * @class DataDelegateSimple
 * @brief prepare data for simple particle dynamics.
 */
template <class ParticlesType = BaseParticles>
class DataDelegateSimple
{
  public:
    explicit DataDelegateSimple(SPHBody &sph_body)
        : particles_(DynamicCast<ParticlesType>(this, &sph_body.getBaseParticles())),
          sorted_id_(sph_body.getBaseParticles().sorted_id_),
          unsorted_id_(sph_body.getBaseParticles().unsorted_id_){};
    virtual ~DataDelegateSimple(){};
    ParticlesType *getParticles() { return particles_; };

  protected:
    ParticlesType *particles_;
    StdLargeVec<size_t> &sorted_id_;
    StdLargeVec<size_t> &unsorted_id_;
};

/**
 * @class DataDelegateInner
 * @brief prepare data for inner particle dynamics
 */
template <class ParticlesType = BaseParticles,
          class BaseDataDelegateType = DataDelegateSimple<ParticlesType>>
class DataDelegateInner : public BaseDataDelegateType
{
  public:
    explicit DataDelegateInner(BaseInnerRelation &inner_relation)
        : BaseDataDelegateType(inner_relation.getSPHBody()),
          inner_configuration_(inner_relation.inner_configuration_){};
    virtual ~DataDelegateInner(){};

  protected:
    /** inner configuration of the designated body */
    ParticleConfiguration &inner_configuration_;
};

/**
 * @class DataDelegateContact
 * @brief prepare data for contact particle dynamics
 */
template <class ParticlesType = BaseParticles,
          class ContactParticlesType = BaseParticles,
          class BaseDataDelegateType = DataDelegateSimple<ParticlesType>>
class DataDelegateContact : public BaseDataDelegateType
{
  public:
    explicit DataDelegateContact(BaseContactRelation &body_contact_relation);
    virtual ~DataDelegateContact(){};
    void addExtraContactRelation(SPHBody &this_body, BaseContactRelation &extra_contact_relation);

  protected:
    SPHBodyVector contact_bodies_;
    StdVec<ContactParticlesType *> contact_particles_;
    /** Configurations for particle interaction between bodies. */
    StdVec<ParticleConfiguration *> contact_configuration_;
};

/**
 * @class DataDelegateComplex
 * @brief prepare data for complex particle dynamics
 */
template <class ParticlesType = BaseParticles,
          class ContactParticlesType = BaseParticles>
class DataDelegateComplex : public DataDelegateInner<ParticlesType>,
                            public DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>
{
  public:
    explicit DataDelegateComplex(ComplexRelation &complex_relation)
        : DataDelegateInner<ParticlesType>(complex_relation.getInnerRelation()),
          DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>(complex_relation.getContactRelation()){};
    virtual ~DataDelegateComplex(){};
};

/**
 * @class BaseInteractionComplex
 * @brief Abstract base class for general complex interaction dynamics
 */
template <class InteractionInnerType, class ContactDataType>
class BaseInteractionComplex : public InteractionInnerType, public ContactDataType
{
  public:
    // template for different combination of constructing body relations
    template <typename... Args>
    BaseInteractionComplex(BaseContactRelation &contact_relation,
                           BaseInnerRelation &inner_relation, Args &&...args)
        : InteractionInnerType(inner_relation, std::forward<Args>(args)...),
          ContactDataType(contact_relation){};
    template <typename... Args>
    BaseInteractionComplex(ComplexRelation &complex_relation, Args &&...args)
        : BaseInteractionComplex(complex_relation.getContactRelation(),
                                 complex_relation.getInnerRelation(), std::forward<Args>(args)...){};
    template <typename... Args>
    BaseInteractionComplex(BaseContactRelation &extra_contact_relation,
                           ComplexRelation &complex_relation, Args &&...args)
        : BaseInteractionComplex(complex_relation, std::forward<Args>(args)...)
    {
        this->addExtraContactRelation(this->sph_body_, extra_contact_relation);
    };
    virtual ~BaseInteractionComplex(){};
};
} // namespace SPH
#endif // BASE_PARTICLE_DYNAMICS_H