#include "base_particle_generator.h"

#include "adaptation.h"
#include "all_io.h"
#include "base_body.h"
#include "base_particles.h"

#include "io_vtk.h"
namespace SPH
{
//=================================================================================================//
ParticleGenerator<BaseParticles>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
    : sph_body_(sph_body), base_particles_(base_particles),
      particle_spacing_ref_(sph_body.getSPHAdaptation().ReferenceSpacing()) {}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::addParticlePosition(const Vecd &position)
{
    position_.push_back(position);
}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::generateParticlesWithGeometricVariables()
{
    prepareGeometricData();
    setAllParticleBounds();
    initializeDiscreteVariables();
}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::setAllParticleBounds()
{
    base_particles_.initializeAllParticlesBounds(position_.size());
}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::addPositionAndVolumetricMeasure(
    const Vecd &position, Real volumetric_measure)
{
    addParticlePosition(position);
    volumetric_measure_.push_back(volumetric_measure);
}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::initializeDiscreteVariables()
{
    base_particles_.registerPositionAndVolumetricMeasure(position_, volumetric_measure_);
}
//=================================================================================================//
void ParticleGenerator<BaseParticles>::initializeDiscreteVariablesFromReload()
{
    base_particles_.registerPositionAndVolumetricMeasureFromReload();
}
//=================================================================================================//
ParticleGenerator<SurfaceParticles>::
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
    : ParticleGenerator<BaseParticles>(sph_body, surface_particles),
      surface_particles_(surface_particles) {}
//=================================================================================================//
void ParticleGenerator<SurfaceParticles>::addSurfaceProperties(const Vecd &surface_normal, Real thickness)
{
    surface_normal_.push_back(surface_normal);
    surface_thickness_.push_back(thickness);
}
//=================================================================================================//
void ParticleGenerator<SurfaceParticles>::initializeDiscreteVariables()
{
    ParticleGenerator<BaseParticles>::initializeDiscreteVariables();
    surface_particles_.registerSurfaceProperties(surface_normal_, surface_thickness_);
}
//=================================================================================================//
void ParticleGenerator<SurfaceParticles>::initializeDiscreteVariablesFromReload()
{
    ParticleGenerator<BaseParticles>::initializeDiscreteVariablesFromReload();
    surface_particles_.registerSurfacePropertiesFromReload();
}
//=================================================================================================//
ParticleGenerator<ObserverParticles>::ParticleGenerator(
    SPHBody &sph_body, BaseParticles &base_particles, const StdVec<Vecd> &positions)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles), positions_(positions) {}
//=================================================================================================//
void ParticleGenerator<ObserverParticles>::prepareGeometricData()
{
    for (size_t i = 0; i < positions_.size(); ++i)
    {
        addPositionAndVolumetricMeasure(positions_[i], 0.0);
    }
}
//=================================================================================================//
ParticleGenerator<ObserverParticles>::~ParticleGenerator()
{
    BodyStatesRecordingToVtp write_observer(sph_body_);
    write_observer.writeToFile();
}
//=================================================================================================//
} // namespace SPH
