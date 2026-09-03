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
    checkBodyIncluded(bodies_, &sph_body);
    sph_body.getBaseParticles().addVariableToWrite<DataType>(name);
    return *this;
}
//=================================================================================================//
template <typename DerivedVariableMethod, typename DynamicsIdentifier, typename... Args>
BodyStatesRecording &BodyStatesRecording::addDerivedVariableRecording(
    DynamicsIdentifier &identifier, Args &&...args)
{
    SPHBody &sph_body = identifier.getSPHBody();
    checkBodyIncluded(bodies_, &sph_body);
    derived_variables_.push_back(derived_variables_keeper_.createPtr<DerivedVariableMethod>(
        identifier, std::forward<Args>(args)...));
    return *this;
}
//=================================================================================================//
template <typename DataType>
void ReloadParticleIO::addToReload(SPHBody &sph_body, const std::string &name)
{
    checkBodyIncluded(bodies_, &sph_body);
    sph_body.getBaseParticles().addEvolvingVariable<DataType>(name);
}
//=================================================================================================//
} // namespace SPH
#endif // IO_BASE_HPP
