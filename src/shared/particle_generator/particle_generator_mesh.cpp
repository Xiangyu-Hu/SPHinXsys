#include "particle_generator_mesh.h"

#include "base_particle_dynamics.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<UnstructuredMesh>::GeneratingMethod(ANSYSMesh &ansys_mesh)
    : elements_centroids_(ansys_mesh.elements_centroids_),
      elements_volumes_(ansys_mesh.elements_volumes_) {}
//=================================================================================================//
ParticleGenerator<BaseParticles, UnstructuredMesh>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, ANSYSMesh &ansys_mesh)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles),
      GeneratingMethod<UnstructuredMesh>(ansys_mesh) {}
//=================================================================================================//
void ParticleGenerator<BaseParticles, UnstructuredMesh>::prepareGeometricData()
{
    for (size_t i = 0; i != elements_centroids_.size(); ++i)
    {
        addPositionAndVolumetricMeasure(elements_centroids_[i], elements_volumes_[i]);
    }
}
//=================================================================================================//
} // namespace SPH
