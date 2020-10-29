/**
* @file 	base_particle_dynamics.hpp
* @brief 	This is the implementation of the template class for 3D build
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_particle_dynamics.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template <class BodyType, 
			  class ParticlesType, 
		      class MaterialType,
			  class ContactBodyType, 
			  class ContactParticlesType, 
			  class ContactMaterialType>
	DataDelegateContact<BodyType, ParticlesType, MaterialType, ContactBodyType, ContactParticlesType, ContactMaterialType>
		::DataDelegateContact(SPHBodyContactRelation* body_contact_relation) :	
		body_(dynamic_cast<BodyType*>(body_contact_relation->sph_body_)),
		particles_(dynamic_cast<ParticlesType*>(body_->base_particles_)),
		material_(dynamic_cast<MaterialType*>(body_->base_particles_->base_material_)),
		sorted_id_(body_->base_particles_->sorted_id_),
		unsorted_id_(body_->base_particles_->unsorted_id_),
		contact_configuration_(body_contact_relation->contact_configuration_)
	{
		SPHBodyVector contact_sph_bodies = body_contact_relation->contact_sph_bodies_;
		for (size_t i = 0; i != contact_sph_bodies.size(); ++i) {
			contact_bodies_.push_back(dynamic_cast<ContactBodyType*>(contact_sph_bodies[i]));
			contact_particles_.push_back(dynamic_cast<ContactParticlesType*>(contact_sph_bodies[i]->base_particles_));
			contact_material_.push_back(dynamic_cast<ContactMaterialType*>(contact_sph_bodies[i]->base_particles_->base_material_));
		}
	}
	//=================================================================================================//
	template <class BodyType,
		class ParticlesType,
		class MaterialType,
		class ContactBodyType,
		class ContactParticlesType,
		class ContactMaterialType>
	DataDelegateComplex<BodyType, ParticlesType, MaterialType, ContactBodyType, ContactParticlesType, ContactMaterialType>
		::DataDelegateComplex(SPHBodyComplexRelation* body_complex_relation) :
		body_(dynamic_cast<BodyType*>(body_complex_relation->sph_body_)),
		particles_(dynamic_cast<ParticlesType*>(body_->base_particles_)),
		material_(dynamic_cast<MaterialType*>(body_->base_particles_->base_material_)),
		inner_configuration_(body_complex_relation->inner_configuration_),
		sorted_id_(body_->base_particles_->sorted_id_),
		unsorted_id_(body_->base_particles_->unsorted_id_),
		contact_configuration_(body_complex_relation->contact_configuration_)
	{
		SPHBodyVector contact_sph_bodies = body_complex_relation->contact_sph_bodies_;
		for (size_t i = 0; i != contact_sph_bodies.size(); ++i) {
			contact_bodies_.push_back(dynamic_cast<ContactBodyType*>(contact_sph_bodies[i]));
			contact_particles_.push_back(dynamic_cast<ContactParticlesType*>(contact_sph_bodies[i]->base_particles_));
			contact_material_.push_back(dynamic_cast<ContactMaterialType*>(contact_sph_bodies[i]->base_particles_->base_material_));
		}
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType>& reduce_functor, ReduceOperation& reduce_operation, Real dt)
	{
		for (size_t i = 0; i < number_of_particles; ++i)
		{
			temp = reduce_operation(temp, reduce_functor(i, dt));
		}
		return temp;
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType>& reduce_functor, ReduceOperation& reduce_operation, Real dt)
	{
		return parallel_reduce(blocked_range<size_t>(0, number_of_particles),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					temp0 = reduce_operation(temp0, reduce_functor(i, dt));
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y)->ReturnType {
				return reduce_operation(x, y);
			}
			);
	}
	//=================================================================================================//
}
//=================================================================================================//
