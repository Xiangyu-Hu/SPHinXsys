#include "particle_generator_lattice.h"

#include "adaptation.h"
#include "base_body.h"
#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<Lattice>::GeneratingMethod(SPHBody &sph_body)
    : sph_adaptation_(sph_body.getSPHAdaptation()), lattice_spacing_(sph_adaptation_.MinimumSpacing()),
      domain_bounds_(sph_body.getSPHSystemBounds()), initial_shape_(sph_body.getInitialShape())
{
    if (!initial_shape_.isValid())
    {
        std::cout << "\n BaseParticleGeneratorLattice Error: initial_shape_ is invalid." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
}
//=================================================================================================//
ParticleGenerator<BaseParticles, Lattice>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, Shape &target_shape)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles),
      GeneratingMethod<Lattice>(sph_body), target_shape_(target_shape) {}
//=================================================================================================//
ParticleGenerator<BaseParticles, Lattice>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
    : ParticleGenerator(sph_body, base_particles, sph_body.getInitialShape()) {}
//=================================================================================================//
void ParticleGenerator<BaseParticles, Lattice>::
    addPositionAndVolumetricMeasure(const Vecd &position, Real volume)
{
    Real local_particle_spacing = sph_adaptation_.getLocalSpacing(target_shape_, position);
    Real local_particle_volume_ratio = pow(lattice_spacing_ / local_particle_spacing, Dimensions);
    if (rand_uniform(0.0, 1.0) < local_particle_volume_ratio)
    {
        ParticleGenerator<BaseParticles>::addPositionAndVolumetricMeasure(
            position, volume / local_particle_volume_ratio);
    }
}
//=================================================================================================//
ParticleGenerator<SurfaceParticles, Lattice>::
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, Real thickness)
    : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
      GeneratingMethod<Lattice>(sph_body), total_volume_(0), thickness_(thickness),
      particle_spacing_(sph_body.getSPHAdaptation().ReferenceSpacing()),
      avg_particle_volume_(pow(particle_spacing_, Dimensions - 1) * thickness_),
      all_cells_(0), planned_number_of_particles_(0)
{
    lattice_spacing_ = thickness_ > particle_spacing_
                           ? 0.5 * particle_spacing_
                           : 0.5 * thickness_;
}
//=================================================================================================//
} // namespace SPH
