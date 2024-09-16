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
} // namespace SPH
#endif // RELATION_CK_HPP
