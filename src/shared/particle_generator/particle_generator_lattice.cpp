#include "particle_generator_lattice.h"

#include "adaptation.h"
#include "base_body.h"
#include "complex_shape.h"
#include "solid_particles.h"

namespace SPH
{
//=================================================================================================//
BaseParticleGeneratorLattice::BaseParticleGeneratorLattice(SPHBody &sph_body)
    : lattice_spacing_(sph_body.sph_adaptation_->ReferenceSpacing()),
      domain_bounds_(sph_body.getSPHSystemBounds()),
      body_shape_(*sph_body.body_shape_)
{
    if (!body_shape_.isValid())
    {
        std::cout << "\n BaseParticleGeneratorLattice Error: body_shape_ is invalid." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
}
//=================================================================================================//
ParticleGeneratorLattice::ParticleGeneratorLattice(SPHBody &sph_body)
    : BaseParticleGeneratorLattice(sph_body), ParticleGenerator(sph_body) {}
//=================================================================================================//
ParticleGeneratorMultiResolution::ParticleGeneratorMultiResolution(SPHBody &sph_body, Shape &target_shape)
    : ParticleGeneratorLattice(sph_body), target_shape_(target_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, sph_body.sph_adaptation_)),
      h_ratio_(*base_particles_.getVariableByName<Real>("SmoothingLengthRatio"))
{
    lattice_spacing_ = particle_adaptation_->MinimumSpacing();
}
//=================================================================================================//
ParticleGeneratorMultiResolution::ParticleGeneratorMultiResolution(SPHBody &sph_body)
    : ParticleGeneratorMultiResolution(sph_body, *sph_body.body_shape_) {}
//=================================================================================================//
void ParticleGeneratorMultiResolution::
    initializePositionAndVolumetricMeasure(const Vecd &position, Real volume)
{
    Real local_particle_spacing = particle_adaptation_->getLocalSpacing(target_shape_, position);
    Real local_particle_volume_ratio = pow(lattice_spacing_ / local_particle_spacing, Dimensions);
    if ((Real)rand() / (RAND_MAX) < local_particle_volume_ratio)
    {
        ParticleGeneratorLattice::initializePositionAndVolumetricMeasure(position, volume / local_particle_volume_ratio);
        initializeSmoothingLengthRatio(local_particle_spacing);
    }
}
//=================================================================================================//
void ParticleGeneratorMultiResolution::initializeSmoothingLengthRatio(Real local_spacing)
{
    h_ratio_.push_back(particle_adaptation_->ReferenceSpacing() / local_spacing);
}
//=================================================================================================//
ParticleGeneratorSplitAndMerge::ParticleGeneratorSplitAndMerge(SPHBody &sph_body)
    : ParticleGeneratorLattice(sph_body),
      particle_adaptation_(DynamicCast<ParticleSplitAndMerge>(this, sph_body.sph_adaptation_)),
      h_ratio_(*base_particles_.getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void ParticleGeneratorSplitAndMerge::
    initializePositionAndVolumetricMeasure(const Vecd &position, Real volume)
{
    ParticleGeneratorLattice::initializePositionAndVolumetricMeasure(position, volume);
    h_ratio_.push_back(1.0);
}
//=================================================================================================//
ThickSurfaceParticleGeneratorLattice::
    ThickSurfaceParticleGeneratorLattice(SPHBody &sph_body, Real global_avg_thickness)
    : BaseParticleGeneratorLattice(sph_body), SurfaceParticleGenerator(sph_body),
      total_volume_(0), global_avg_thickness_(global_avg_thickness),
      particle_spacing_(sph_body.sph_adaptation_->ReferenceSpacing()),
      avg_particle_volume_(pow(particle_spacing_, Dimensions - 1) * global_avg_thickness_),
      all_cells_(0), planned_number_of_particles_(0)
{
    lattice_spacing_ = global_avg_thickness_ > particle_spacing_ ? 0.5 * particle_spacing_ : 0.5 * global_avg_thickness_;
}
//=================================================================================================//
} // namespace SPH
