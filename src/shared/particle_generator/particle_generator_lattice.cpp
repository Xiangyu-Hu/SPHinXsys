#include "particle_generator_lattice.h"

#include "adaptation.h"
#include "base_body.h"
#include "complex_shape.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<Lattice>::GeneratingMethod(SPHBody &sph_body)
    : lattice_spacing_(sph_body.sph_adaptation_->ReferenceSpacing()),
      domain_bounds_(sph_body.getSPHSystemBounds()),
      initial_shape_(sph_body.getInitialShape())
{
    if (!initial_shape_.isValid())
    {
        std::cout << "\n BaseParticleGeneratorLattice Error: initial_shape_ is invalid." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
}
//=================================================================================================//
ParticleGenerator<Lattice>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body), GeneratingMethod<Lattice>(sph_body) {}
//=================================================================================================//
ParticleGenerator<Lattice, Adaptive>::ParticleGenerator(SPHBody &sph_body, Shape &target_shape)
    : ParticleGenerator<Lattice>(sph_body), target_shape_(target_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, sph_body.sph_adaptation_))
{
    lattice_spacing_ = particle_adaptation_->MinimumSpacing();
}
//=================================================================================================//
ParticleGenerator<Lattice, Adaptive>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Lattice, Adaptive>(sph_body, sph_body.getInitialShape()) {}
//=================================================================================================//
void ParticleGenerator<Lattice, Adaptive>::
    initializePositionAndVolumetricMeasure(const Vecd &position, Real volume)
{
    Real local_particle_spacing = particle_adaptation_->getLocalSpacing(target_shape_, position);
    Real local_particle_volume_ratio = pow(lattice_spacing_ / local_particle_spacing, Dimensions);
    if (rand_uniform(0.0, 1.0) < local_particle_volume_ratio)
    {
        ParticleGenerator<Base>::initializePositionAndVolumetricMeasure(
            position, volume / local_particle_volume_ratio);
    }
}
//=================================================================================================//
ParticleGenerator<ThickSurface, Lattice>::ParticleGenerator(SPHBody &sph_body, Real thickness)
    : ParticleGenerator<Surface>(sph_body), GeneratingMethod<Lattice>(sph_body),
      total_volume_(0), thickness_(thickness),
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
