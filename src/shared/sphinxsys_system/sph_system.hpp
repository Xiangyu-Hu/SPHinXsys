#ifndef SPH_SYSTEM_HPP
#define SPH_SYSTEM_HPP

#include "sph_system.h"

#include "base_variable.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DataType *SPHSystem::registerSystemVariable(const std::string &name, DataType initial_value)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(all_system_variables_,
                                                 all_system_variable_ptrs_, name, initial_value)
                     ->ValueAddress();
}
//=================================================================================================//
template <typename DataType, class ExecutionPolicy>
SingularVariable<DataType> *SPHSystem::registerSystemVariable(const ExecutionPolicy &execution_policy,
                                                              const std::string &name, DataType initial_value)
{
    registerSystemVariable<DataType>(name, initial_value);
    return getSystemVariableByName<DataType>(name);
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *SPHSystem::registerSystemVariable(const ParallelDevicePolicy &execution_policy,
                                                              const std::string &name, DataType initial_value)
{
    registerSystemVariable<DataType>(name, initial_value);
    return getSystemVariableByName<DataType>(execution_policy, name);
}
//=================================================================================================//
template <typename DataType>
DataType *SPHSystem::getSystemVariableDataByName(const std::string &name)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the system variable '" << name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }

    return variable->ValueAddress();
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *SPHSystem::getSystemVariableByName(const std::string &name)
{
    SingularVariable<DataType> *variable =
        findVariableByName<DataType>(all_system_variables_, name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the system variable '" << name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }

    return variable;
}
//=================================================================================================//
template <typename DataType, class ExecutionPolicy>
SingularVariable<DataType> *SPHSystem::getSystemVariableByName(
    const ExecutionPolicy &execution_policy, const std::string &name)
{
    return getSystemVariableByName<DataType>(name);
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *SPHSystem::getSystemVariableByName(
    const ParallelDevicePolicy &execution_policy, const std::string &name)
{
    SingularVariable<DataType> *variable = getSystemVariableByName<DataType>(name);
    if (variable->isValueDelegated())
    {
        unique_system_variable_ptrs_.createPtr<SingularDeviceSharedVariable<DataType>>(variable);
    }

    return variable;
}
//=================================================================================================//
} // namespace SPH
#endif // SPH_SYSTEM_HPP
