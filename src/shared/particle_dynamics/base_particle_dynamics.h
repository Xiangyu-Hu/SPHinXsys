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
class DataDelegateSimple
{
  public:
    explicit DataDelegateSimple(SPHBody &sph_body)
        : sph_body_(sph_body),
          particles_(&sph_body.getBaseParticles()){};
    virtual ~DataDelegateSimple(){};
    SPHBody &getSPHBody() { return sph_body_; };
    BaseParticles *getParticles() { return particles_; };

  protected:
    SPHBody &sph_body_;
    BaseParticles *particles_;
};

/**
 * @class DataDelegateInner
 * @brief prepare data for inner particle dynamics
 */
template <class BaseDataDelegateType>
class BaseDataDelegateInner : public BaseDataDelegateType
{
    BaseInnerRelation &inner_relation_;

  public:
    explicit BaseDataDelegateInner(BaseInnerRelation &inner_relation)
        : BaseDataDelegateType(inner_relation.getSPHBody()),
          inner_relation_(inner_relation),
          inner_configuration_(inner_relation.inner_configuration_){};
    virtual ~BaseDataDelegateInner(){};
    BaseInnerRelation &getBodyRelation() { return inner_relation_; };

  protected:
    /** inner configuration of the designated body */
    ParticleConfiguration &inner_configuration_;
};
using DataDelegateInner = BaseDataDelegateInner<DataDelegateSimple>;
using DataDelegateInnerOnly = BaseDataDelegateInner<DataDelegateEmptyBase>;
/**
 * @class DataDelegateContact
 * @brief prepare data for contact particle dynamics
 */
template <class BaseDataDelegateType>
class BaseDataDelegateContact : public BaseDataDelegateType
{
    BaseContactRelation &contact_relation_;

  public:
    explicit BaseDataDelegateContact(BaseContactRelation &contact_relation)
        : BaseDataDelegateType(contact_relation.getSPHBody()),
          contact_relation_(contact_relation)
    {
        RealBodyVector contact_sph_bodies = contact_relation.contact_bodies_;
        for (size_t i = 0; i != contact_sph_bodies.size(); ++i)
        {
            contact_bodies_.push_back(contact_sph_bodies[i]);
            contact_particles_.push_back(&contact_sph_bodies[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
    };
    virtual ~BaseDataDelegateContact(){};
    BaseContactRelation &getBodyRelation() { return contact_relation_; };

  protected:
    SPHBodyVector contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    /** Configurations for particle interaction between bodies. */
    StdVec<ParticleConfiguration *> contact_configuration_;
};
using DataDelegateContact = BaseDataDelegateContact<DataDelegateSimple>;
/** Required for projection-based operator splitting (implicit) method,
 * the error and parameter need to computed for inner and contact particles together
 * other than can be handled by separately as other explicit methods. */
using DataDelegateContactOnly = BaseDataDelegateContact<DataDelegateEmptyBase>;
} // namespace SPH
#endif // BASE_PARTICLE_DYNAMICS_H
