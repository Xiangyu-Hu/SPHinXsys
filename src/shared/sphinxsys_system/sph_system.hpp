#ifndef SPH_SYSTEM_HPP
#define SPH_SYSTEM_HPP

#include "sph_system.h"

#include "sphinxsys_variable.h"

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

    if (variable == nullptr)
    {
        std::cout << "\nError: the system variable '" << name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }

    return variable;
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

    return variable->Data();
}
//=================================================================================================//
} // namespace SPH
#endif // SPH_SYSTEM_HPP
