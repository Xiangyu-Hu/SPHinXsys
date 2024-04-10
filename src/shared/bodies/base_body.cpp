#include "base_body.h"

#include "base_body_relation.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> initial_shape_ptr, const std::string &body_name)
    : sph_system_(sph_system), body_name_(body_name), newly_updated_(true), base_particles_(nullptr),
      initial_shape_(initial_shape_ptr_keeper_.assignPtr(initial_shape_ptr)),
      sph_adaptation_(sph_adaptation_ptr_keeper_.createPtr<SPHAdaptation>(*this)),
      base_material_(nullptr)
{
    sph_system_.sph_bodies_.push_back(this);
}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> initial_shape_ptr)
    : SPHBody(sph_system, initial_shape_ptr, initial_shape_ptr->getName()) {}
//=================================================================================================//
BoundingBox SPHBody::getSPHSystemBounds()
{
    return sph_system_.system_domain_bounds_;
}
//=================================================================================================//
SPHSystem &SPHBody::getSPHSystem()
{
    return sph_system_;
}
//=================================================================================================//
BaseParticles &SPHBody::getBaseParticles()
{
    if (base_particles_ == nullptr)
    {
        std::cout << "\n Error: BaseParticle not generated yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_particles_;
};
//=================================================================================================//
BaseMaterial &SPHBody::getBaseMaterial()
{
    if (base_material_ == nullptr)
    {
        std::cout << "\n Error: BaseMaterial not generated yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_material_;
};
//=================================================================================================//
void SPHBody::allocateConfigurationMemoriesForBufferParticles()
{
    for (size_t i = 0; i < body_relations_.size(); i++)
    {
        body_relations_[i]->resizeConfiguration();
    }
}
//=================================================================================================//
BoundingBox SPHBody::getBodyShapeBounds()
{
    return initial_shape_->getBounds();
}
//=================================================================================================//
void SPHBody::defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
{
    sph_adaptation_->resetAdaptationRatios(h_spacing_ratio, new_system_refinement_ratio);
}
//=================================================================================================//
void SPHBody::writeParticlesToVtuFile(std::ostream &output_file)
{
    base_particles_->writeParticlesToVtk(output_file);
}
//=================================================================================================//
void SPHBody::writeParticlesToVtpFile(std::ofstream &output_file)
{
    base_particles_->writeParticlesToVtk(output_file);
}
//=================================================================================================//
void SPHBody::writeSurfaceParticlesToVtuFile(std::ofstream &output_file, BodySurface &surface_particles)
{
    base_particles_->writeSurfaceParticlesToVtuFile(output_file, surface_particles);
}
//=================================================================================================//
void SPHBody::writeParticlesToPltFile(std::ofstream &output_file)
{
    base_particles_->writeParticlesToPltFile(output_file);
}
//=================================================================================================//
void SPHBody::writeParticlesToXmlForRestart(std::string &filefullpath)
{
    base_particles_->writeParticlesToXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::readParticlesFromXmlForRestart(std::string &filefullpath)
{
    base_particles_->readParticleFromXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::writeToXmlForReloadParticle(std::string &filefullpath)
{
    base_particles_->writeToXmlForReloadParticle(filefullpath);
}
//=================================================================================================//
void SPHBody::readFromXmlForReloadParticle(std::string &filefullpath)
{
    base_particles_->readFromXmlForReloadParticle(filefullpath);
}
//=================================================================================================//
} // namespace SPH
