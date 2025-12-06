#include "particle_generator_lattice.h"

#include "base_body.h"
#include "base_mesh.hpp"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
void ParticleGenerator<BaseParticles, Lattice>::prepareGeometricData()
{
    Mesh mesh(domain_bounds_, lattice_spacing_, 0);
    Real particle_volume = lattice_spacing_ * lattice_spacing_;
    Arrayi number_of_lattices = mesh.AllCells();
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
        {
            Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j));
            if (initial_shape_.checkNotFar(particle_position, lattice_spacing_))
            {
                if (initial_shape_.checkContain(particle_position))
                {
                    addPositionAndVolumetricMeasure(particle_position, particle_volume);
                }
            }
        }
}
//=================================================================================================//
void ParticleGenerator<SurfaceParticles, Lattice>::prepareGeometricData()
{
    // Calculate the total volume and
    // count the number of cells inside the body volume, where we might put particles.
    Mesh mesh(domain_bounds_, lattice_spacing_, 0);
    Arrayi number_of_lattices = mesh.AllCells();
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
        {
            Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j));
            if (initial_shape_.checkNotFar(particle_position, lattice_spacing_))
            {
                if (initial_shape_.checkContain(particle_position))
                {
                    all_cells_++;
                    total_volume_ += lattice_spacing_ * lattice_spacing_;
                }
            }
        }
    Real number_of_particles = total_volume_ / avg_particle_volume_ + 0.5;
    planned_number_of_particles_ = int(number_of_particles);

    // Calculate the interval based on the number of particles.
    Real interval = (Real)planned_number_of_particles_ / (Real)all_cells_;
    if (interval <= 0)
        interval = 1; // It has to be lager than 0.

    // initialize a uniform distribution between 0 (inclusive) and 1 (exclusive)
    std::mt19937_64 rng;
    std::uniform_real_distribution<Real> unif(0, 1);

    // Add a particle in each interval, randomly. We will skip the last intervals if we already reach the number of particles
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
        {
            Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j));
            if (initial_shape_.checkNotFar(particle_position, lattice_spacing_))
            {
                if (initial_shape_.checkContain(particle_position))
                {
                    Real random_real = unif(rng);
                    // If the random_real is smaller than the interval, add a particle, only if we haven't reached the max. number of particles
                    if (random_real <= interval && base_particles_.TotalRealParticles() < planned_number_of_particles_)
                    {
                        addPositionAndVolumetricMeasure(particle_position, avg_particle_volume_ / thickness_);
                        addSurfaceProperties(initial_shape_.findNormalDirection(particle_position), thickness_);
                    }
                }
            }
        }
}
//=================================================================================================//
} // namespace SPH
