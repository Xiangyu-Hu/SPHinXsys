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
#include "base_data_type_package.h"
#include "base_dynamics.h"
#include "neighborhood.h"
#include "sphinxsys_containers.h"

namespace SPH
{
/**
 * @class DataDelegateInner
 * @brief prepare data for inner particle dynamics
 */
class DataDelegateInner
{
    BaseInnerRelation &inner_relation_;

  public:
    explicit DataDelegateInner(BaseInnerRelation &inner_relation)
        : inner_relation_(inner_relation),
          inner_configuration_(inner_relation.inner_configuration_) {};
    virtual ~DataDelegateInner() {};
    BaseInnerRelation &getBodyRelation() { return inner_relation_; };

  protected:
    /** inner configuration of the designated body */
    ParticleConfiguration &inner_configuration_;
};

/**
 * @class DataDelegateContact
 * @brief prepare data for contact particle dynamics
 */
class DataDelegateContact
{
    BaseContactRelation &contact_relation_;

  public:
    explicit DataDelegateContact(BaseContactRelation &contact_relation)
        : contact_relation_(contact_relation)
    {
        RealBodyVector contact_sph_bodies = contact_relation.contact_bodies_;
        for (size_t i = 0; i != contact_sph_bodies.size(); ++i)
        {
            contact_bodies_.push_back(contact_sph_bodies[i]);
            contact_particles_.push_back(&contact_sph_bodies[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
    };
    virtual ~DataDelegateContact() {};
    BaseContactRelation &getBodyRelation() { return contact_relation_; };

  protected:
    SPHBodyVector contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    /** Configurations for particle interaction between bodies. */
    StdVec<ParticleConfiguration *> contact_configuration_;
};
} // namespace SPH
#endif // BASE_PARTICLE_DYNAMICS_H
