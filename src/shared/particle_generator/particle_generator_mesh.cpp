#include "particle_generator_mesh.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<UnstructuredMesh>::GeneratingMethod(ANSYSMesh &ansys_mesh)
    : elements_centroids_(ansys_mesh.elements_centroids_),
      elements_volumes_(ansys_mesh.elements_volumes_) {}
//=================================================================================================//
GeneratingMethod<UnstructuredMesh>::GeneratingMethod(ANSYSMesh_3d& ansys_mesh_3d)
    : elements_centroids_(ansys_mesh_3d.elements_centroids_), // Assume ANSYSMesh_3d has a similar interface
    elements_volumes_(ansys_mesh_3d.elements_volumes_) {}
//=================================================================================================//
ParticleGenerator<UnstructuredMesh>::ParticleGenerator(SPHBody &sph_body, ANSYSMesh &ansys_mesh)
    : ParticleGenerator<Base>(sph_body), GeneratingMethod<UnstructuredMesh>(ansys_mesh) {}
//=================================================================================================//
ParticleGenerator<UnstructuredMesh>::ParticleGenerator(SPHBody& sph_body, ANSYSMesh_3d& ansys_mesh_3d)
    : ParticleGenerator<Base>(sph_body), GeneratingMethod<UnstructuredMesh>(ansys_mesh_3d) {}
//=================================================================================================//
void ParticleGenerator<UnstructuredMesh>::initializeGeometricVariables()
{
    for (size_t i = 0; i != elements_centroids_.size(); ++i)
    {
        initializePositionAndVolumetricMeasure(elements_centroids_[i], elements_volumes_[i]);
    }
}
//=================================================================================================//
} // namespace SPH
