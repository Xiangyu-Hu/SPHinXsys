#include "base_particle_generator.h"

#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<ObserverParticles>::ParticleGenerator(
    SPHBody &sph_body, BaseParticles &base_particles, TriangleMeshShape &triangle_mesh_shape)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles)
{
    StdVec<std::array<Real, 3>> &vertices = triangle_mesh_shape.getVertices();
    for (size_t i = 0; i != vertices.size(); ++i)
    {
        positions_.push_back(Vec3d(vertices[i][0], vertices[i][1], vertices[i][2]));
    }
}
//=================================================================================================//
} // namespace SPH
