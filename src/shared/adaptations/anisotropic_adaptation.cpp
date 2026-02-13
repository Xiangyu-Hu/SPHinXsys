#include "anisotropic_adaptation.h"

namespace SPH
{
//=================================================================================================//
AnisotropicAdaptation::AnisotropicAdaptation(
    const Vecd &scaling, const Vecd &orientation,
    Real global_resolution, Real h_spacing_ratio, Real refinement_to_global)
    : SPHAdaptation(global_resolution, h_spacing_ratio, refinement_to_global),
      scaling_(scaling), orientation_(orientation), deformation_matrix_(Matd::Identity()),
      spacing_ref_min_(spacing_ref_ * scaling_.cwiseMin()), h_ref_max_(h_ref_ * scaling_.cwiseMax())
{
    // rotation from local reference normal to orientation
    deformation_matrix_ = scaling.cwiseInverse().asDiagonal() *
                          Eigen::Quaterniond::FromTwoVectors(local_n0, orientation).toRotationMatrix();
    spacing_min_ = spacing_ref_min_;
}
//=================================================================================================//
UniquePtr<BaseCellLinkedList> AnisotropicAdaptation::
    createCellLinkedList(const BoundingBoxd &domain_bounds, BaseParticles &base_particles)
{
    return makeUnique<CellLinkedList<AnisotropicAdaptation>>(
        domain_bounds, kernel_ptr_->KernelSize() * h_ref_max_, base_particles, *this);
}
//=================================================================================================//
UniquePtr<LevelSet> AnisotropicAdaptation::createLevelSet(Shape &shape, Real refinement) const
{
    // estimate the required minimum grid spacing for the level set
    int total_levels = (int)log10(shape.getBounds().MinimumDimension() / spacing_ref_min_) + 2;
    Real coarsest_spacing = spacing_ref_min_ * pow(2.0, total_levels - 1);
    LevelSet coarser_level_sets(shape.getBounds(), coarsest_spacing / refinement,
                                total_levels - 1, shape, *this, refinement);
    // return the finest level set only
    return makeUnique<LevelSet>(shape.getBounds(), &coarser_level_sets, shape, *this, refinement);
}
//=================================================================================================//
} // namespace SPH
