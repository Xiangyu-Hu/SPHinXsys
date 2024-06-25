#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.h"
#include "io_all.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Base>::ParticleGenerator(SPHBody &sph_body)
    : base_particles_(sph_body.getBaseParticles()),
      particle_spacing_ref_(sph_body.sph_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
void ParticleGenerator<Base>::addParticlePosition(const Vecd &position)
{
    position_.push_back(position);
    base_particles_.total_real_particles_++;
}
//=================================================================================================//
void ParticleGenerator<Base>::generateParticlesWithGeometricVariables()
{
    prepareGeometricData();
    setAllParticleBounds();
    initializeGeometricParticleVariables();
}
//=================================================================================================//
void ParticleGenerator<Base>::setAllParticleBounds()
{
    base_particles_.initializeAllParticlesBounds(position_.size());
}
//=================================================================================================//
void ParticleGenerator<Base>::addPositionAndVolumetricMeasure(
    const Vecd &position, Real volumetric_measure)
{
    addParticlePosition(position);
    volumetric_measure_.push_back(volumetric_measure);
}
//=================================================================================================//
void ParticleGenerator<Base>::initializeGeometricParticleVariables()
{
    base_particles_.registerSharedVariableFrom<Vecd>("Position", position_);
    base_particles_.registerSharedVariableFrom<Real>("VolumetricMeasure", volumetric_measure_);
}
//=================================================================================================//
ParticleGenerator<Surface>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body) {}
//=================================================================================================//
void ParticleGenerator<Surface>::addSurfaceProperties(const Vecd &surface_normal, Real thickness)
{
    surface_normal_.push_back(surface_normal);
    surface_thickness_.push_back(thickness);
}
//=================================================================================================//
void ParticleGenerator<Surface>::initializeGeometricParticleVariables()
{
    ParticleGenerator<Base>::initializeGeometricParticleVariables();
    base_particles_.registerSharedVariableFrom<Vecd>("NormalDirection", surface_normal_);
    base_particles_.registerSharedVariableFrom<Real>("Thickness", surface_thickness_);
}
//=================================================================================================//
void ParticleGenerator<Observer>::prepareGeometricData()
{
    for (size_t i = 0; i < positions_.size(); ++i)
    {
        addPositionAndVolumetricMeasure(positions_[i], 0.0);
    }
}
//=================================================================================================//
ParticleGenerator<Reload>::ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name)
    : ParticleGenerator<Base>(sph_body), base_material_(sph_body.getBaseMaterial())
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
void ParticleGenerator<Reload>::prepareGeometricData()
{
    base_particles_.readReloadXmlFile(file_path_);
}
//=================================================================================================//
void ParticleGenerator<Reload>::setAllParticleBounds()
{
    base_particles_.initializeAllParticlesBounds(file_path_);
}
//=================================================================================================//
void ParticleGenerator<Reload>::initializeGeometricParticleVariables()
{
    base_particles_.readFromXmlForReloadParticle(file_path_);
}
//=================================================================================================//
} // namespace SPH
