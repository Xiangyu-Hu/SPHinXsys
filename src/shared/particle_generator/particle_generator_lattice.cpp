#include "particle_generator_lattice.h"

#include "adaptation.h"
#include "base_body.h"
#include "complex_shape.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<Lattice>::GeneratingMethod(SPHBody &sph_body, Shape &contain_shape)
    : lattice_spacing_(sph_body.sph_adaptation_->ReferenceSpacing()),
      domain_bounds_(sph_body.getSPHSystemBounds()), contain_shape_(contain_shape)
{
    if (!contain_shape_.isValid())
    {
        std::cout << "\n BaseParticleGeneratorLattice Error: shape_ is invalid." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
}
//=================================================================================================//
ParticleGenerator<BaseParticles, Lattice>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, Shape &contain_shape)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles),
      GeneratingMethod<Lattice>(sph_body, contain_shape) {}
//=================================================================================================//
ParticleGenerator<BaseParticles, Lattice, Adaptive>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      Shape &contain_shape, Shape &adaptation_shape)
    : ParticleGenerator<BaseParticles, Lattice>(sph_body, base_particles, contain_shape),
      adaptation_shape_(adaptation_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, sph_body.sph_adaptation_))
{
    lattice_spacing_ = particle_adaptation_->MinimumSpacing();
}
//=================================================================================================//
ParticleGenerator<BaseParticles, Lattice, Adaptive>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, Shape &common_shape)
    : ParticleGenerator<BaseParticles, Lattice, Adaptive>(
          sph_body, base_particles, common_shape, common_shape) {}
//=================================================================================================//
void ParticleGenerator<BaseParticles, Lattice, Adaptive>::
    addPositionAndVolumetricMeasure(const Vecd &position, Real volume)
{
    Real local_particle_spacing = particle_adaptation_->getLocalSpacing(adaptation_shape_, position);
    Real local_particle_volume_ratio = pow(lattice_spacing_ / local_particle_spacing, Dimensions);
    if (rand_uniform(0.0, 1.0) < local_particle_volume_ratio)
    {
        ParticleGenerator<BaseParticles>::addPositionAndVolumetricMeasure(
            position, volume / local_particle_volume_ratio);
    }
}
//=================================================================================================//
ParticleGenerator<SurfaceParticles, Lattice>::
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                      Shape &contain_shape, Real thickness)
    : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
      GeneratingMethod<Lattice>(sph_body, contain_shape), total_volume_(0), thickness_(thickness),
      particle_spacing_(sph_body.sph_adaptation_->ReferenceSpacing()),
      avg_particle_volume_(pow(particle_spacing_, Dimensions - 1) * thickness_),
      all_cells_(0), planned_number_of_particles_(0)
{
    lattice_spacing_ = thickness_ > particle_spacing_
                           ? 0.5 * particle_spacing_
                           : 0.5 * thickness_;
}
//=================================================================================================//
} // namespace SPH
