#include "adaptation.h"

#include "base_particle_dynamics.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"
#include "level_set.h"
#include "mesh_with_data_packages.hpp"
#include "vector_functions.h"

namespace SPH
{
//=================================================================================================//
SPHAdaptation::SPHAdaptation(Real resolution_ref, Real h_spacing_ratio, Real system_refinement_ratio)
    : h_spacing_ratio_(h_spacing_ratio), system_refinement_ratio_(system_refinement_ratio),
      local_refinement_level_(0), spacing_ref_(resolution_ref / system_refinement_ratio_),
      h_ref_(h_spacing_ratio_ * spacing_ref_), kernel_ptr_(makeUnique<KernelWendlandC2>(h_ref_)),
      sigma0_ref_(computeLatticeNumberDensity(Vecd())),
      spacing_min_(this->MostRefinedSpacingRegular(spacing_ref_, local_refinement_level_)),
      Vol_min_(pow(spacing_min_, Dimensions)), h_ratio_max_(spacing_ref_ / spacing_min_){};
//=================================================================================================//
Real SPHAdaptation::MostRefinedSpacing(Real coarse_particle_spacing, int local_refinement_level)
{
    return MostRefinedSpacingRegular(coarse_particle_spacing, local_refinement_level);
}
Real SPHAdaptation::MostRefinedSpacingRegular(Real coarse_particle_spacing, int local_refinement_level)
{
    return coarse_particle_spacing / pow(2.0, local_refinement_level);
}
//=================================================================================================//
Real SPHAdaptation::computeLatticeNumberDensity(Vec2d zero)
{
    Real sigma(0);
    Real cutoff_radius = kernel_ptr_->CutOffRadius();
    Real particle_spacing = ReferenceSpacing();
    int search_depth = int(cutoff_radius / particle_spacing) + 1;
    for (int j = -search_depth; j <= search_depth; ++j)
        for (int i = -search_depth; i <= search_depth; ++i)
        {
            Vec2d particle_location(i * particle_spacing, j * particle_spacing);
            Real distance = particle_location.norm();
            if (distance < cutoff_radius)
                sigma += kernel_ptr_->W(distance, particle_location);
        }
    return sigma;
}
//=================================================================================================//
Real SPHAdaptation::computeLatticeNumberDensity(Vec3d zero)
{
    Real sigma(0);
    Real cutoff_radius = kernel_ptr_->CutOffRadius();
    Real particle_spacing = ReferenceSpacing();
    int search_depth = int(cutoff_radius / particle_spacing) + 1;
    for (int k = -search_depth; k <= search_depth; ++k)
        for (int j = -search_depth; j <= search_depth; ++j)
            for (int i = -search_depth; i <= search_depth; ++i)
            {
                Vec3d particle_location(i * particle_spacing,
                                        j * particle_spacing, k * particle_spacing);
                Real distance = particle_location.norm();
                if (distance < cutoff_radius)
                    sigma += kernel_ptr_->W(distance, particle_location);
            }
    return sigma;
}
//=================================================================================================//
Real SPHAdaptation::NumberDensityScaleFactor(Real smoothing_length_ratio)
{
    return pow(smoothing_length_ratio, Dimensions);
}
//=================================================================================================//
void SPHAdaptation::resetAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
{
    h_spacing_ratio_ = h_spacing_ratio;
    spacing_ref_ = spacing_ref_ * system_refinement_ratio_ / new_system_refinement_ratio;
    system_refinement_ratio_ = new_system_refinement_ratio;
    h_ref_ = h_spacing_ratio_ * spacing_ref_;
    getKernel()->resetSmoothingLength(h_ref_);
    sigma0_ref_ = computeLatticeNumberDensity(Vecd());
    spacing_min_ = MostRefinedSpacing(spacing_ref_, local_refinement_level_);
    Vol_min_ = pow(spacing_min_, Dimensions);
    h_ratio_max_ = spacing_ref_ / spacing_min_;
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> SPHAdaptation::createCellLinkedList(const BoundingBox &domain_bounds)
{
    return makeUnique<CellLinkedList>(domain_bounds, kernel_ptr_->CutOffRadius(), *this);
}
//=================================================================================================//
UniquePtr<BaseLevelSet> SPHAdaptation::createLevelSet(Shape &shape, Real refinement_ratio)
{
    // estimate the required mesh levels
    int total_levels = (int)log10(MinimumDimension(shape.getBounds()) / ReferenceSpacing()) + 2;
    Real coarsest_spacing = ReferenceSpacing() * pow(2.0, total_levels - 1);
    MultilevelLevelSet coarser_level_sets(shape.getBounds(), coarsest_spacing / refinement_ratio,
                                          total_levels - 1, shape, *this);
    // return the finest level set only
    return makeUnique<RefinedMesh<LevelSet>>(shape.getBounds(), *coarser_level_sets.getMeshLevels().back(), shape, *this);
}
//=================================================================================================//
ParticleWithLocalRefinement::
    ParticleWithLocalRefinement(Real resolution_ref, Real h_spacing_ratio, Real system_refinement_ratio,
                                int local_refinement_level)
    : SPHAdaptation(resolution_ref, h_spacing_ratio, system_refinement_ratio), h_ratio_(nullptr)
{
    local_refinement_level_ = local_refinement_level;
    spacing_min_ = MostRefinedSpacingRegular(spacing_ref_, local_refinement_level_);
    Vol_min_ = pow(spacing_min_, Dimensions);
    h_ratio_max_ = spacing_ref_ / spacing_min_;
    // To ensure that the adaptation strictly within all level set and mesh cell linked list levels
    finest_spacing_bound_ = spacing_min_ + Eps;
    coarsest_spacing_bound_ = spacing_ref_ - Eps;
}
//=================================================================================================//
void ParticleWithLocalRefinement::initializeAdaptationVariables(BaseParticles &base_particles)
{
    SPHAdaptation::initializeAdaptationVariables(base_particles);
    h_ratio_ = base_particles.registerSharedVariable<Real>(
        "SmoothingLengthRatio", [&](size_t i) -> Real
        { return ReferenceSpacing() / base_particles.ParticleSpacing(i); });
    base_particles.addVariableToSort<Real>("SmoothingLengthRatio");
    base_particles.addVariableToReload<Real>("SmoothingLengthRatio");
}
//=================================================================================================//
size_t ParticleWithLocalRefinement::getCellLinkedListTotalLevel()
{
    return size_t(local_refinement_level_);
}
//=================================================================================================//
size_t ParticleWithLocalRefinement::getLevelSetTotalLevel()
{
    return getCellLinkedListTotalLevel() + 1;
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> ParticleWithLocalRefinement::createCellLinkedList(const BoundingBox &domain_bounds)
{
    return makeUnique<MultilevelCellLinkedList>(domain_bounds, kernel_ptr_->CutOffRadius(),
                                                getCellLinkedListTotalLevel(), *this);
}
//=================================================================================================//
UniquePtr<BaseLevelSet> ParticleWithLocalRefinement::createLevelSet(Shape &shape, Real refinement_ratio)
{
    return makeUnique<MultilevelLevelSet>(shape.getBounds(), ReferenceSpacing() / refinement_ratio,
                                          getLevelSetTotalLevel(), shape, *this);
}
//=================================================================================================//
Real ParticleRefinementByShape::smoothedSpacing(const Real &measure, const Real &transition_thickness)
{
    Real ratio_ref = measure / (2.0 * transition_thickness);
    Real target_spacing = coarsest_spacing_bound_;
    if (ratio_ref < kernel_ptr_->KernelSize())
    {
        Real weight = kernel_ptr_->W_1D(ratio_ref) / kernel_ptr_->W_1D(0.0);
        target_spacing = weight * finest_spacing_bound_ + (1.0 - weight) * coarsest_spacing_bound_;
    }
    return target_spacing;
}
//=================================================================================================//
Real ParticleRefinementNearSurface::getLocalSpacing(Shape &shape, const Vecd &position)
{
    Real phi = fabs(shape.findSignedDistance(position));
    return smoothedSpacing(phi, spacing_ref_);
}
//=================================================================================================//
Real ParticleRefinementWithinShape::getLocalSpacing(Shape &shape, const Vecd &position)
{
    Real phi = shape.findSignedDistance(position);
    return phi < 0.0 ? finest_spacing_bound_ : smoothedSpacing(phi, 2.0 * spacing_ref_);
}
//=================================================================================================//
} // namespace SPH
