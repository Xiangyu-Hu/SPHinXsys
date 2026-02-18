#include "adaptation.hpp"

#include "base_particles.hpp"
#include "cell_linked_list.h"
#include "level_set_shape.h"

namespace SPH
{
//=================================================================================================//
SPHAdaptation::SPHAdaptation(Real global_resolution, Real h_spacing_ratio, Real refinement_to_global)
    : global_resolution_(global_resolution), h_spacing_ratio_(h_spacing_ratio),
      refinement_to_global_(refinement_to_global), local_refinement_level_(0),
      spacing_ref_(global_resolution / refinement_to_global_),
      h_ref_(h_spacing_ratio_ * spacing_ref_), kernel_ptr_(makeUnique<KernelWendlandC2>(h_ref_)),
      sigma0_ref_(computeLatticeNumberDensity(Vecd())),
      spacing_min_(MostRefinedSpacing(spacing_ref_, local_refinement_level_)),
      Vol_min_(pow(spacing_min_, Dimensions)), h_ratio_max_(spacing_ref_ / spacing_min_) {};
//=================================================================================================//
Real SPHAdaptation::MostRefinedSpacing(Real spacing_ref, int local_refinement_level)
{
    return spacing_ref / pow(2.0, local_refinement_level);
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
void SPHAdaptation::resetAdaptationRatios(Real h_spacing_ratio, Real new_refinement_to_global)
{
    h_spacing_ratio_ = h_spacing_ratio;
    spacing_ref_ = spacing_ref_ * refinement_to_global_ / new_refinement_to_global;
    refinement_to_global_ = new_refinement_to_global;
    h_ref_ = h_spacing_ratio_ * spacing_ref_;
    getKernel()->resetSmoothingLength(h_ref_);
    sigma0_ref_ = computeLatticeNumberDensity(Vecd());
    spacing_min_ = MostRefinedSpacing(spacing_ref_, local_refinement_level_);
    Vol_min_ = pow(spacing_min_, Dimensions);
    h_ratio_max_ = spacing_ref_ / spacing_min_;
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> SPHAdaptation::
    createCellLinkedList(const BoundingBoxd &domain_bounds, BaseParticles &base_particles)
{
    return makeUnique<CellLinkedList<SPHAdaptation>>(
        domain_bounds, kernel_ptr_->CutOffRadius(), base_particles, *this);
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> SPHAdaptation::createFinestCellLinkedList(
    const BoundingBoxd &domain_bounds, BaseParticles &base_particles)
{
    Real grid_spacing = kernel_ptr_->CutOffRadius() / pow(2.0, local_refinement_level_);
    return makeUnique<CellLinkedList<SPHAdaptation>>(
        domain_bounds, grid_spacing, base_particles, *this);
}
//=================================================================================================//
UniquePtr<LevelSet> SPHAdaptation::createLevelSet(Shape &shape, Real refinement) const
{
    // estimate the required mesh levels
    int total_levels = (int)log10(shape.getBounds().MinimumDimension() / ReferenceSpacing()) + 2;
    Real coarsest_spacing = ReferenceSpacing() * pow(2.0, total_levels - 1);
    LevelSet coarser_level_sets(shape.getBounds(), coarsest_spacing / refinement,
                                total_levels - 1, shape, *this, refinement);
    // return the finest level set only
    return makeUnique<LevelSet>(shape.getBounds(), &coarser_level_sets, shape, *this, refinement);
}
//=================================================================================================//
Real AdaptiveByShape::smoothedSpacing(const Real &measure, const Real &transition_thickness)
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
AdaptiveSmoothingLength::AdaptiveSmoothingLength(
    Real global_resolution, Real h_spacing_ratio, Real refinement_to_global, int local_refinement_level)
    : SPHAdaptation(global_resolution, h_spacing_ratio, refinement_to_global),
      dv_h_ratio_(nullptr), dv_h_level_(nullptr), h_ratio_(nullptr), h_level_(nullptr)
{
    local_refinement_level_ = local_refinement_level;
    spacing_min_ = MostRefinedSpacing(spacing_ref_, local_refinement_level_);
    Vol_min_ = pow(spacing_min_, Dimensions);
    h_ratio_max_ = spacing_ref_ / spacing_min_;
    // To ensure that the adaptation strictly within all level set and mesh cell linked list levels
    finest_spacing_bound_ = spacing_min_ + Eps;
    coarsest_spacing_bound_ = spacing_ref_ - Eps;
    max_cut_off_radius_ = kernel_ptr_->KernelSize() * h_ref_;
}
//=================================================================================================//
void AdaptiveSmoothingLength::initializeAdaptationVariables(BaseParticles &base_particles)
{
    SPHAdaptation::initializeAdaptationVariables(base_particles);
    dv_h_ratio_ = base_particles.registerStateVariable<Real>(
        "SmoothingLengthRatio", [&](size_t i) -> Real
        { return ReferenceSpacing() / base_particles.ParticleSpacing(i); });
    dv_h_level_ = base_particles.registerStateVariable<int>("SmoothingLengthLevel");
    h_ratio_ = dv_h_ratio_->Data();
    h_level_ = dv_h_level_->Data();
    base_particles.addEvolvingVariable<Real>("SmoothingLengthRatio");
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> AdaptiveSmoothingLength::
    createCellLinkedList(const BoundingBoxd &domain_bounds, BaseParticles &base_particles)
{
    return makeUnique<CellLinkedList<AdaptiveSmoothingLength>>(
        domain_bounds, kernel_ptr_->CutOffRadius(),
        local_refinement_level_, base_particles, *this);
}
//=================================================================================================//
UniquePtr<LevelSet> AdaptiveSmoothingLength::createLevelSet(Shape &shape, Real refinement) const
{
    // one more level for interpolation
    return makeUnique<LevelSet>(shape.getBounds(), ReferenceSpacing() / refinement,
                                local_refinement_level_ + 1, shape, *this, refinement);
}
//=================================================================================================//
AdaptiveSmoothingLength::SmoothedSpacing::SmoothedSpacing(AdaptiveSmoothingLength &encloser)
    : smoothing_kernel_(*encloser.kernel_ptr_),
      kernel_size_(smoothing_kernel_.KernelSize()), inv_w0_(1.0 / smoothing_kernel_.normalized_W(0)),
      finest_spacing_bound_(encloser.finest_spacing_bound_),
      coarsest_spacing_bound_(encloser.coarsest_spacing_bound_) {}
//=================================================================================================//
Real AdaptiveNearSurface::getLocalSpacing(Shape &shape, const Vecd &position)
{
    Real phi = fabs(shape.findSignedDistance(position));
    return smoothedSpacing(phi, spacing_ref_);
}
//=================================================================================================//
AdaptiveNearSurface::LocalSpacing::LocalSpacing(
    AdaptiveNearSurface &encloser, LevelSetShape &level_set_shape)
    : smoothed_spacing_(encloser), level_set_(level_set_shape.getLevelSet()),
      spacing_ref_(encloser.spacing_ref_) {}
//=================================================================================================//
Real AdaptiveWithinShape::getLocalSpacing(Shape &shape, const Vecd &position)
{
    Real phi = shape.findSignedDistance(position);
    return phi < 0.0 ? finest_spacing_bound_ : smoothedSpacing(phi, 2.0 * spacing_ref_);
}
//=================================================================================================//
AdaptiveWithinShape::LocalSpacing::LocalSpacing(
    AdaptiveWithinShape &encloser, LevelSetShape &level_set_shape)
    : smoothed_spacing_(encloser), level_set_(level_set_shape.getLevelSet()),
      spacing_ref_(encloser.spacing_ref_) {}
//=================================================================================================//
AnisotropicAdaptation::AnisotropicAdaptation(
    Real global_resolution, Real h_spacing_ratio_, Real refinement_to_global, int local_refinement_level)
    : AdaptiveSmoothingLength(
          global_resolution, h_spacing_ratio_, refinement_to_global, local_refinement_level),
      dv_scaling_(nullptr), dv_orientation_(nullptr),
      dv_deformation_matrix_(nullptr), dv_deformation_det_(nullptr) {}
//=================================================================================================//
void AnisotropicAdaptation::initializeAdaptationVariables(BaseParticles &particles)
{
    AdaptiveSmoothingLength::initializeAdaptationVariables(particles);
    dv_scaling_ = particles.registerStateVariable<Vecd>("AnisotropicScaling", Vecd::Ones().eval());
    dv_orientation_ = particles.registerStateVariable<Vecd>("AnisotropicOrientation", Vecd::Zero().eval());
    dv_deformation_matrix_ = particles.registerStateVariable<Matd>("AnisotropicMatrix", Matd::Identity().eval());
    dv_deformation_det_ = particles.registerStateVariable<Real>("AnisotropicDeterminate", 1.0);
    particles.addVariableToWrite<Vecd>(dv_scaling_);
    particles.addVariableToWrite<Vecd>(dv_orientation_);
}
//=================================================================================================//
PrescribedAnisotropy::PrescribedAnisotropy(
    const Vecd &scaling, const Vecd &orientation,
    Real global_resolution, Real h_spacing_ratio_, Real refinement_to_global)
    : AnisotropicAdaptation(global_resolution, h_spacing_ratio_, refinement_to_global, 0),
      scaling_ref_(scaling), orientation_ref_(orientation)
{
    deformation_matrix_ref_ =
        scaling_ref_.cwiseInverse().asDiagonal() * RotationMatrix(Vecd::UnitX(), orientation_ref_);
    max_cut_off_radius_ = kernel_ptr_->KernelSize() * h_ref_ * scaling_ref_.maxCoeff();
}
//=================================================================================================//
void PrescribedAnisotropy::initializeAdaptationVariables(BaseParticles &particles)
{
    AnisotropicAdaptation::initializeAdaptationVariables(particles);
    UnsignedInt total_real_particles = particles.TotalRealParticles();
    dv_scaling_->fill(0, total_real_particles, [&](size_t i) -> Vecd
                      { return scaling_ref_; });
    dv_orientation_->fill(0, total_real_particles, [&](size_t i) -> Vecd
                          { return orientation_ref_; });
    dv_deformation_matrix_->fill(0, total_real_particles, [&](size_t i) -> Matd
                                 { return deformation_matrix_ref_; });
    dv_deformation_det_->fill(0, total_real_particles, [&](size_t i) -> Real
                              { return deformation_matrix_ref_.determinant(); });
}
//=================================================================================================//
} // namespace SPH
