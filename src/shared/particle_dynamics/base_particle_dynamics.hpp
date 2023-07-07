/**
 * @file 	base_particle_dynamics.hpp
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef BASE_PARTICLE_DYNAMICS_HPP
#define BASE_PARTICLE_DYNAMICS_HPP

#include "base_body.h"
#include "base_particle_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class BaseDataDelegateType>
DataDelegateContact<ParticlesType, ContactParticlesType, BaseDataDelegateType>::
    DataDelegateContact(BaseContactRelation &body_contact_relation)
    : BaseDataDelegateType(body_contact_relation.getSPHBody())
{
    RealBodyVector contact_sph_bodies = body_contact_relation.contact_bodies_;
    contact_configuration_device_ = makeSharedDevice<StdSharedVec<NeighborhoodDevice*>>(contact_sph_bodies.size(),
        execution::executionQueue.getQueue());    
    for (size_t i = 0; i != contact_sph_bodies.size(); ++i)
    {
        contact_bodies_.push_back(contact_sph_bodies[i]);
        contact_particles_.push_back(DynamicCast<ContactParticlesType>(this, &contact_sph_bodies[i]->getBaseParticles()));
        contact_configuration_.push_back(&body_contact_relation.contact_configuration_[i]);
        contact_configuration_device_->at(i) = body_contact_relation.contact_configuration_device_.at(i).data();
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

    auto new_contact_configuration_device = makeSharedDevice<StdSharedVec<NeighborhoodDevice*>>(
        contact_configuration_device_->size() + extra_contact_relation.contact_bodies_.size(),
        execution::executionQueue.getQueue());
    std::copy(contact_configuration_device_->begin(), contact_configuration_device_->end(),
        new_contact_configuration_device->begin());
    for (size_t i = 0; i != extra_contact_relation.contact_bodies_.size(); ++i)
    {
        contact_configuration_.push_back(&extra_contact_relation.contact_configuration_[i]);
        new_contact_configuration_device->at(contact_configuration_device_->size() + i) =
                extra_contact_relation.contact_configuration_device_.at(i).data();
    }
    contact_configuration_device_ = std::move(new_contact_configuration_device);
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
#endif // BASE_PARTICLE_DYNAMICS_HPP
