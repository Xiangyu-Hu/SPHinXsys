/**
* @file 	base_particle_dynamics.hpp
* @brief 	This is the implementation of the template class for 3D build
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef BASE_PARTICLE_DYNAMICS_HPP
#define BASE_PARTICLE_DYNAMICS_HPP



#include "base_particle_dynamics.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template <class BodyType, 
			  class ParticlesType, 
		      class MaterialType,
			  class ContactBodyType, 
			  class ContactParticlesType, 
			  class ContactMaterialType,
			  class BaseDataDelegateType>
	DataDelegateContact<BodyType, ParticlesType, MaterialType, ContactBodyType, ContactParticlesType, ContactMaterialType, BaseDataDelegateType>
		::DataDelegateContact(BaseBodyRelationContact &body_contact_relation) :
		BaseDataDelegateType(*body_contact_relation.sph_body_)
	{
		RealBodyVector contact_sph_bodies = body_contact_relation.contact_bodies_;
		for (size_t i = 0; i != contact_sph_bodies.size(); ++i) {
			contact_bodies_.push_back(DynamicCast<ContactBodyType>(this, contact_sph_bodies[i]));
			contact_particles_.push_back(DynamicCast<ContactParticlesType>(this, contact_sph_bodies[i]->base_particles_));
			contact_material_.push_back(DynamicCast<ContactMaterialType>(this, contact_sph_bodies[i]->base_particles_->base_material_));
			contact_configuration_.push_back(&body_contact_relation.contact_configuration_[i]);
		}
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t total_real_particles, ReturnType temp,
		ReduceFunctor<ReturnType>& reduce_functor, ReduceOperation& reduce_operation, Real dt)
	{
		for (size_t i = 0; i < total_real_particles; ++i)
		{
			temp = reduce_operation(temp, reduce_functor(i, dt));
		}
		return temp;
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t total_real_particles, ReturnType temp,
		ReduceFunctor<ReturnType>& reduce_functor, ReduceOperation& reduce_operation, Real dt)
	{
		return parallel_reduce(blocked_range<size_t>(0, total_real_particles),
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
	template<class ParticleDynamicsInnerType, class ContactDataType>
	ParticleDynamicsComplex<ParticleDynamicsInnerType, ContactDataType>:: 
		ParticleDynamicsComplex(ComplexBodyRelation &complex_relation, 
			BaseBodyRelationContact &extra_contact_relation) :
		ParticleDynamicsInnerType(complex_relation.inner_relation_),
		ContactDataType(complex_relation.contact_relation_)
	{
		if (complex_relation.sph_body_ != extra_contact_relation.sph_body_)
		{
			std::cout << "\n Error: the two body_realtions do not have the same source body!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		for (auto& extra_body : extra_contact_relation.contact_bodies_)
		{
			// here we first obtain the pointer to the most derived class and then implicitly downcast it to
			// the types defined in the base complex dynamics  
			this->contact_bodies_.push_back(extra_body->ThisObjectPtr());
			this->contact_particles_.push_back(extra_body->base_particles_->ThisObjectPtr());
			this->contact_material_.push_back(extra_body->base_particles_->base_material_->ThisObjectPtr());
		}

		for (size_t i = 0; i != extra_contact_relation.contact_bodies_.size(); ++i)
		{
			this->contact_configuration_.push_back(&extra_contact_relation.contact_configuration_[i]);
		}
	}
	//=================================================================================================//
}
//=================================================================================================//
#endif //BASE_PARTICLE_DYNAMICS_HPP