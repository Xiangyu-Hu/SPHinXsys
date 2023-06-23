#include "particle_generator_lattice.h"

#include "base_body.h"
#include "base_mesh.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "complex_shape.h"

namespace SPH
{
//=================================================================================================//
void ParticleGeneratorLattice::initializeGeometricVariables()
{
    BaseMesh mesh(domain_bounds_, lattice_spacing_, 0);
    Real particle_volume = lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
    Arrayi number_of_lattices = mesh.AllCellsFromAllGridPoints(mesh.AllGridPoints());
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
            for (int k = 0; k < number_of_lattices[2]; ++k)
            {
                Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j, k));
                if (body_shape_.checkNotFar(particle_position, lattice_spacing_))
                {
                    if (body_shape_.checkContain(particle_position))
                    {
                        initializePositionAndVolumetricMeasure(particle_position, particle_volume);
                    }
                }
            }
}
//=================================================================================================//
void ThickSurfaceParticleGeneratorLattice::initializeGeometricVariables()
{
    // Calculate the total volume and
    // count the number of cells inside the body volume, where we might put particles.
    std::unique_ptr<BaseMesh> mesh(new BaseMesh(domain_bounds_, lattice_spacing_, 0));
    Arrayi number_of_lattices = mesh->AllCellsFromAllGridPoints(mesh->AllGridPoints());
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
            for (int k = 0; k < number_of_lattices[2]; ++k)
            {
                Vecd particle_position = mesh->CellPositionFromIndex(Arrayi(i, j, k));
                if (body_shape_.checkNotFar(particle_position, lattice_spacing_))
                {
                    if (body_shape_.checkContain(particle_position))
                    {
                        all_cells_++;
                        total_volume_ += lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
                    }
                }
            }
    Real number_of_particles = total_volume_ / avg_particle_volume_ + 0.5;
    planned_number_of_particles_ = int(number_of_particles);

    // Calculate the interval based on the number of particles.
    Real interval = planned_number_of_particles_ / (all_cells_ + TinyReal);
    if (interval <= 0)
        interval = 1; // It has to be lager than 0.
    // Add a particle in each interval, randomly. We will skip the last intervals if we already reach the number of particles.
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
            for (int k = 0; k < number_of_lattices[2]; ++k)
            {
                Vecd particle_position = mesh->CellPositionFromIndex(Arrayi(i, j, k));
                if (body_shape_.checkNotFar(particle_position, lattice_spacing_))
                {
                    if (body_shape_.checkContain(particle_position))
                    {
                        Real random_real = (Real)rand() / (RAND_MAX);
                        // If the random_real is smaller than the interval, add a particle, only if we haven't reached the max. number of particles.
                        if (random_real <= interval && base_particles_.total_real_particles_ < planned_number_of_particles_)
                        {
                            initializePositionAndVolumetricMeasure(particle_position, avg_particle_volume_ / global_avg_thickness_);
                            initializeSurfaceProperties(body_shape_.findNormalDirection(particle_position), global_avg_thickness_);
                        }
                    }
                }
            }
}
//=================================================================================================//
} // namespace SPH
