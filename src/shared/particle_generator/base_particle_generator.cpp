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
    initializeParticleVariables();
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
void ParticleGenerator<Base>::initializeParticleVariables()
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
void ParticleGenerator<Surface>::initializeParticleVariables()
{
    ParticleGenerator<Base>::initializeParticleVariables();
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
} // namespace SPH
