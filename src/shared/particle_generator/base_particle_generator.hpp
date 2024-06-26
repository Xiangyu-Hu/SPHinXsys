#ifndef BASE_PARTICLE_GENERATOR_HPP
#define BASE_PARTICLE_GENERATOR_HPP

#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
template <typename... BaseGeneratorParameters>
ParticleGenerator<Reload, BaseGeneratorParameters...>::
    ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name)
    : ParticleGenerator<BaseGeneratorParameters...>(sph_body)
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
template <typename... BaseGeneratorParameters>
void ParticleGenerator<Reload, BaseGeneratorParameters...>::prepareGeometricData()
{
    this->base_particles_.readReloadXmlFile(file_path_);
}
//=================================================================================================//
template <typename... BaseGeneratorParameters>
void ParticleGenerator<Reload, BaseGeneratorParameters...>::setAllParticleBounds()
{
    this->base_particles_.initializeAllParticlesBounds(file_path_);
};
//=================================================================================================//
template <typename... BaseGeneratorParameters>
void ParticleGenerator<Reload, BaseGeneratorParameters...>::initializeParticleVariables()
{
    ParticleGenerator<BaseGeneratorParameters...>::initializeParticleVariablesFromReload();
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_HPP
