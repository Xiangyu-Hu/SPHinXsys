/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	base_particle_dynamics.hpp
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define  
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef BASE_PARTICLE_DYNAMICS_HPP
#define BASE_PARTICLE_DYNAMICS_HPP

#include "base_particle_dynamics.h"
#include "base_body.h"

namespace SPH
{
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, class BaseDataDelegateType>
	DataDelegateContact<ParticlesType, ContactParticlesType, BaseDataDelegateType>::
		DataDelegateContact(BaseContactRelation &body_contact_relation)
		: BaseDataDelegateType(body_contact_relation.getSPHBody())
	{
		RealBodyVector contact_sph_bodies = body_contact_relation.contact_bodies_;
		for (size_t i = 0; i != contact_sph_bodies.size(); ++i)
		{
			contact_bodies_.push_back(contact_sph_bodies[i]);
			contact_particles_.push_back(DynamicCast<ContactParticlesType>(this, &contact_sph_bodies[i]->getBaseParticles()));
			contact_configuration_.push_back(&body_contact_relation.contact_configuration_[i]);
		}
	}
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, class BaseDataDelegateType>
	void DataDelegateContact<ParticlesType, ContactParticlesType, BaseDataDelegateType>::
		addExtraContactRelation(SPHBody &this_body, BaseContactRelation &extra_contact_relation)
	{
		if (&this_body != &extra_contact_relation.getSPHBody())
		{
			std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		for (auto &extra_body : extra_contact_relation.contact_bodies_)
		{
			// here we first obtain the pointer to the most derived class and then implicitly downcast it to
			// the types defined in the base complex dynamics
			contact_bodies_.push_back(extra_body);
			contact_particles_.push_back(DynamicCast<ContactParticlesType>(this, &extra_body->getBaseParticles()));
		}

		for (size_t i = 0; i != extra_contact_relation.contact_bodies_.size(); ++i)
		{
			contact_configuration_.push_back(&extra_contact_relation.contact_configuration_[i]);
		}
	}
	//=================================================================================================//
}
//=================================================================================================//
#endif // BASE_PARTICLE_DYNAMICS_HPP