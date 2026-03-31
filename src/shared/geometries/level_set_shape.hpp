#ifndef LEVEL_SET_SHAPE_HPP
#define LEVEL_SET_SHAPE_HPP

#include "level_set_shape.h"

#include "level_set.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
LevelSetShape::LevelSetShape(
    const ExecutionPolicy &ex_policy, SPHSystem &sph_system, const SPHAdaptation &sph_adaptation,
    Shape &shape, Real refinement, UsageType usage_type)
    : LevelSetShape(sph_system, sph_adaptation, shape, refinement)
{
    finishInitialization(ex_policy, usage_type);
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSetShape::finishInitialization(const ExecutionPolicy &ex_policy, UsageType usage_type)
{
    level_set_.finishInitialization(ex_policy, usage_type);
}
//=================================================================================================//
template <typename DataType>
LevelSetShape &LevelSetShape::addPackageVariableToWrite(const std::string &variable_name)
{
    level_set_.addPackageVariableToWrite<DataType>(variable_name);
    return *this;
}
//=================================================================================================//
template <typename DataType>
LevelSetShape &LevelSetShape::addCellVariableToWrite(const std::string &variable_name)
{
    level_set_.addCellVariableToWrite<DataType>(variable_name);
    return *this;
}
//=================================================================================================//
} // namespace SPH
#endif // LEVEL_SET_SHAPE_HPP
