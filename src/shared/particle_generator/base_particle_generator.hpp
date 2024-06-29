#ifndef BASE_PARTICLE_GENERATOR_HPP
#define BASE_PARTICLE_GENERATOR_HPP

#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
template <typename... BaseGeneratorType>
ParticleGenerator<Reload, BaseGeneratorType...>::
    ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name)
    : ParticleGenerator<BaseGeneratorType...>(sph_body)
{
    std::string reload_folder = sph_body.getSPHSystem().getIOEnvironment().reload_folder_;
    if (!fs::exists(reload_folder))
    {
        std::cout << "\n Error: the particle reload folder:" << reload_folder << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    file_path_ = reload_folder + "/" + reload_body_name + "_rld.xml";
}
//=================================================================================================//
template <typename... BaseGeneratorType>
void ParticleGenerator<Reload, BaseGeneratorType...>::prepareGeometricData()
{
    this->base_particles_.readReloadXmlFile(file_path_);
}
//=================================================================================================//
template <typename... BaseGeneratorType>
void ParticleGenerator<Reload, BaseGeneratorType...>::setAllParticleBounds()
{
    this->base_particles_.initializeAllParticlesBoundsFromReloadXml();
};
//=================================================================================================//
template <typename... BaseGeneratorType>
void ParticleGenerator<Reload, BaseGeneratorType...>::initializeParticleVariables()
{
    ParticleGenerator<BaseGeneratorType...>::initializeParticleVariablesFromReload();
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_HPP
