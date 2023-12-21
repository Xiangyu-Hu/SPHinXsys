#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.h"
#include "io_all.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Base>::ParticleGenerator(SPHBody &sph_body)
    : base_particles_(sph_body.getBaseParticles()),
      base_material_(base_particles_.getBaseMaterial()),
      pos_(base_particles_.pos_), Vol_(base_particles_.Vol_),
      unsorted_id_(base_particles_.unsorted_id_) {}
//=================================================================================================//
void ParticleGenerator<Base>::initializePosition(const Vecd &position)
{
    pos_.push_back(position);
    unsorted_id_.push_back(base_particles_.total_real_particles_);
    base_particles_.total_real_particles_++;
}
//=================================================================================================//
void ParticleGenerator<Base>::generateParticlesWithBasicVariables()
{
    initializeGeometricVariables();
    // should be determined first before register other variables
    base_particles_.real_particles_bound_ = base_particles_.total_real_particles_;
    base_material_.registerReloadLocalParameters(&base_particles_);
}
//=================================================================================================//
void ParticleGenerator<Base>::initializePositionAndVolumetricMeasure(
    const Vecd &position, Real volumetric_measure)
{
    initializePosition(position);
    Vol_.push_back(volumetric_measure);
}
//=================================================================================================//
ParticleGenerator<Surface>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body),
      n_(*base_particles_.getVariableByName<Vecd>("NormalDirection")),
      thickness_(*base_particles_.getVariableByName<Real>("Thickness")) {}
//=================================================================================================//
void ParticleGenerator<Surface>::initializeSurfaceProperties(const Vecd &surface_normal, Real thickness)
{
    n_.push_back(surface_normal);
    thickness_.push_back(thickness);
}
//=================================================================================================//
void ParticleGenerator<Observer>::initializeGeometricVariables()
{
    for (size_t i = 0; i < positions_.size(); ++i)
    {
        initializePositionAndVolumetricMeasure(positions_[i], 0.0);
        initializePositionAndVolumetricMeasure(positions_[i], 0.0);
    }
}
//=================================================================================================//
ParticleGeneratorReload::ParticleGeneratorReload(SPHBody &sph_body, const std::string &reload_body_name)
    : ParticleGenerator(sph_body)
{
    std::string reload_folder = sph_body.getSPHSystem().io_environment_->reload_folder_;
    if (!fs::exists(reload_folder))
    {
        std::cout << "\n Error: the particle reload folder:" << reload_folder << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    file_path_ = reload_folder + "/" + reload_body_name + "_rld.xml";
}
//=================================================================================================//
void ParticleGeneratorReload::initializeGeometricVariables()
{
    base_particles_.readFromXmlForReloadParticle(file_path_);
}
//=================================================================================================//
void ParticleGeneratorReload::generateParticlesWithBasicVariables()
{
    base_material_.registerReloadLocalParameters(&base_particles_);
    initializeGeometricVariables();
    // should be determined first before register other variables
    base_particles_.real_particles_bound_ = base_particles_.total_real_particles_;
}
//=================================================================================================//
} // namespace SPH
