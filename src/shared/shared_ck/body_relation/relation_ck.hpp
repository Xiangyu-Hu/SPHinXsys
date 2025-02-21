#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DataType>
DiscreteVariable<DataType> *Relation<Base>::
    addRelationVariable(const std::string &name, size_t data_size)
{
    return relation_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
}
//=================================================================================================//
template <typename DynamicsIdentifier, typename ContactIdentifier>
Relation<Contact<DynamicsIdentifier, ContactIdentifier>>::
    Relation(DynamicsIdentifier &identifier, StdVec<ContactIdentifier *> contact_identifier)
    : Relation<Contact<>>(identifier.getSPHBody()), contact_identifier_(contact_identifier)
{
    StdVec<RealBody *> contact_bodies;
    for (size_t k = 0; k != contact_identifier_.size(); ++k)
    {
        contact_bodies.push_back(DynamicCast<RealBody>(this, &contact_identifier_[k]->getSPHBody()));
    }
    initializeContactRelation(contact_bodies);
}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
