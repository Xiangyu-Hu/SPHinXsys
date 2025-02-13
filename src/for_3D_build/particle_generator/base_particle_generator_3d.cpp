#include "base_particle_generator.h"

#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<ObserverParticles>::ParticleGenerator(
    SPHBody &sph_body, BaseParticles &base_particles, TriangleMeshShape &triangle_mesh_shape)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles)
{
    TriangleMesh &triangle_mesh = *triangle_mesh_shape.getTriangleMesh();
    for (int i = 0; i != triangle_mesh.getNumVertices(); ++i)
    {
        positions_.push_back(SimTKToEigen(triangle_mesh.getVertexPosition(i)));
    }
}
//=================================================================================================//
} // namespace SPH
