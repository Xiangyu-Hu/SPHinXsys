#ifndef IO_BASE_HPP
#define IO_BASE_HPP

#include "io_base.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
BodyStatesRecording &BodyStatesRecording::addToWrite(SPHBody &sph_body, const std::string &name)
{
    if (isBodyIncluded(bodies_, &sph_body))
    {
        sph_body.getBaseParticles().addVariableToWrite<DataType>(name);
    }
    else
    {
        std::cout << "\n Error: the body:" << sph_body.getName()
                  << " is not in the recording list" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *this;
}
//=================================================================================================//
template <typename DerivedVariableMethod, typename DynamicsIdentifier, typename... Args>
BodyStatesRecording &BodyStatesRecording::addDerivedVariableRecording(
    DynamicsIdentifier &identifier, Args &&...args)
{
    SPHBody &sph_body = identifier.getSPHBody();
    if (isBodyIncluded(bodies_, &sph_body))
    {
        derived_variables_.push_back(
            derived_variables_keeper_.createPtr<DerivedVariableMethod>(
                identifier, std::forward<Args>(args)...));
    }
    else
    {
        std::cout << "\n Error: the body:" << sph_body.getName()
                  << " is not in the recording body list" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *this;
}
//=================================================================================================//
template <typename DataType>
void ReloadParticleIO::addToReload(SPHBody &sph_body, const std::string &name)
{
    if (isBodyIncluded(bodies_, &sph_body))
    {
        sph_body.getBaseParticles().addEvolvingVariable<DataType>(name);
    }
    else
    {
        std::cout << "\n Error: the body:" << sph_body.getName()
                  << " is not in the recording list" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // IO_BASE_HPP
