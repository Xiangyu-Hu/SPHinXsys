#ifndef SPH_SYSTEM_HPP
#define SPH_SYSTEM_HPP

#include "sph_system.h"

#include "adaptive_body.h"
#include "ownership.h"
#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *SPHSystem::registerSystemVariable(const std::string &name, DataType initial_value)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);

    return variable != nullptr
               ? variable
               : addVariableToAssemble<DataType>(
                     all_system_variables_, all_system_variable_ptrs_, name, initial_value);
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *SPHSystem::getSystemVariableByName(const std::string &name)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);
    checkPointer(variable, name, "system variable");
    return variable;
}
//=================================================================================================//
template <typename DataType>
DataType *SPHSystem::getSystemVariableDataByName(const std::string &name)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);
    checkPointer(variable, name, "system variable");
    return variable->Data();
}
//=================================================================================================//
template <class BodyType, typename... Args>
BodyType &SPHSystem::addBody(Args &&...args)
{
    return *sph_bodies_keeper_.createPtr<BodyType>(*this, std::forward<Args>(args)...);
}
//=================================================================================================//
template <class BaseBodyType, class AdaptationType, typename... Args>
auto &SPHSystem::addAdaptiveBody(const AdaptationType &adaptation, Args &&...args)
{
    return *sph_bodies_keeper_.createPtr<AdaptiveBody<AdaptationType, BaseBodyType>>(
        *this, adaptation, std::forward<Args>(args)...);
}
//=================================================================================================//
template <class ShapeType, typename... Args>
auto &SPHSystem::addShape(Args &&...args)
{
    return *shapes_keeper_.createPtr<ShapeType>(std::forward<Args>(args)...);
}
//=================================================================================================//
template <class DynamicIdentifier, typename... Args>
auto &SPHSystem::addInnerRelation(DynamicIdentifier &identifier, Args &&...args)
{
    Inner<Relation<DynamicIdentifier>> *relation = relations_keeper_.createPtr<
        Inner<Relation<DynamicIdentifier>>>(identifier, std::forward<Args>(args)...);
    relations_.push_back(relation);
    return *relation;
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Args>
auto &SPHSystem::addContactRelation(
    SourceIdentifier &src_identifier, StdVec<TargetIdentifier *> tar_identifiers, Args &&...args)
{
    Contact<Relation<SourceIdentifier, TargetIdentifier>> *relation =
        relations_keeper_.createPtr<Contact<Relation<SourceIdentifier, TargetIdentifier>>>(
            src_identifier, tar_identifiers, std::forward<Args>(args)...);
    relations_.push_back(relation);
    return *relation;
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Args>
auto &SPHSystem::addContactRelation(
    SourceIdentifier &src_identifier, TargetIdentifier &tar_identifiers, Args &&...args)
{
    return addContactRelation(src_identifier, StdVec<TargetIdentifier *>{&tar_identifiers},
                              std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename DerivedBodyType>
DerivedBodyType &SPHSystem::getBodyByName(const std::string &name)
{
    StdVec<DerivedBodyType *> collected_bodies = collectBodies<DerivedBodyType>();
    for (auto &body : collected_bodies)
    {
        if (body->getName() == name)
        {
            return *DynamicCast<DerivedBodyType>(this, body);
        }
    }
    throw std::runtime_error(
        std::string(type_name<DerivedBodyType>()) + ": " + name + " not found in SPHSystem.");
}
//=================================================================================================//
template <typename DerivedBodyType>
StdVec<DerivedBodyType *> SPHSystem::collectBodies()
{
    StdVec<DerivedBodyType *> collected_bodies;
    for (auto &sph_body : sph_bodies_)
    {
        if (auto casted_body = dynamic_cast<DerivedBodyType *>(sph_body))
        {
            collected_bodies.push_back(casted_body);
        }
    }
    return collected_bodies;
}
//=================================================================================================//
template <typename RelationType>
RelationType *SPHSystem::getRelationByName(const std::string &name)
{
    for (auto &relation : relations_)
    {
        if (relation->getSPHBody()->getName() == name)
        {
            return DynamicCast<RelationType>(this, relation);
        }
    }
    return nullptr;
}
//=================================================================================================//
} // namespace SPH
#endif // SPH_SYSTEM_HPP
