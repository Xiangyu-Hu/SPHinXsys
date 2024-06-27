#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.h"
#include "io_all.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<>::ParticleGenerator(SPHBody &sph_body)
    : base_particles_(sph_body.getBaseParticles()),
      particle_spacing_ref_(sph_body.sph_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
void ParticleGenerator<>::addParticlePosition(const Vecd &position)
{
    position_.push_back(position);
    base_particles_.total_real_particles_++;
}
//=================================================================================================//
void ParticleGenerator<>::generateParticlesWithGeometricVariables()
{
    prepareGeometricData();
    setAllParticleBounds();
    initializeParticleVariables();
}
//=================================================================================================//
void ParticleGenerator<>::setAllParticleBounds()
{
    base_particles_.initializeAllParticlesBounds(position_.size());
}
//=================================================================================================//
void ParticleGenerator<>::addPositionAndVolumetricMeasure(
    const Vecd &position, Real volumetric_measure)
{
    addParticlePosition(position);
    volumetric_measure_.push_back(volumetric_measure);
}
//=================================================================================================//
void ParticleGenerator<>::initializeParticleVariables()
{
    base_particles_.registerSharedVariableFrom<Vecd>("Position", position_);
    base_particles_.registerSharedVariableFrom<Real>("VolumetricMeasure", volumetric_measure_);
    base_particles_.addVariableToReload<Vecd>("Position");
    base_particles_.addVariableToReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void ParticleGenerator<>::initializeParticleVariablesFromReload()
{
    base_particles_.registerSharedVariableFromReload<Vecd>("Position");
    base_particles_.registerSharedVariableFromReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
ParticleGenerator<Surface>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<>(sph_body) {}
//=================================================================================================//
void ParticleGenerator<Surface>::addSurfaceProperties(const Vecd &surface_normal, Real thickness)
{
    surface_normal_.push_back(surface_normal);
    surface_thickness_.push_back(thickness);
}
//=================================================================================================//
void ParticleGenerator<Surface>::initializeParticleVariables()
{
    ParticleGenerator<>::initializeParticleVariables();
    base_particles_.registerSharedVariableFrom<Vecd>("NormalDirection", surface_normal_);
    base_particles_.registerSharedVariableFrom<Real>("Thickness", surface_thickness_);
    base_particles_.addVariableToReload<Vecd>("NormalDirection");
    base_particles_.addVariableToReload<Real>("Thickness");
}
//=================================================================================================//
void ParticleGenerator<Surface>::initializeParticleVariablesFromReload()
{
    ParticleGenerator<>::initializeParticleVariablesFromReload();
    base_particles_.registerSharedVariableFromReload<Vecd>("NormalDirection");
    base_particles_.registerSharedVariableFromReload<Real>("Thickness");
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
} // namespace SPH
