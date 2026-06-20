#ifndef APHI_TEAM7_NATIVE_TEAM7_SHAPES_H
#define APHI_TEAM7_NATIVE_TEAM7_SHAPES_H

#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Coil + plate STL union: air MR spacing is refined vs distance to this shape. */
class AphiTeam7NativeInnerBoundaryShape : public ComplexShape
{
  public:
    explicit AphiTeam7NativeInnerBoundaryShape(const std::string &shape_name, Real stl_scale)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./input/coil.stl", Vec3d::Zero(), stl_scale);
        add<TriangleMeshShapeSTL>("./input/plate.stl", Vec3d::Zero(), stl_scale);
    }
};

/** Tolerant mesh-level lookup for coarse air particles (rc_ref can exceed coarsest grid spacing). */
class Team7TolerantAdaptiveCellLinkedList : public CellLinkedList<AdaptiveSmoothingLength>
{
  public:
    Team7TolerantAdaptiveCellLinkedList(BoundingBoxd tentative_bounds, Real reference_grid_spacing,
                                        UnsignedInt total_levels, BaseParticles &base_particles,
                                        SPHAdaptation &sph_adaptation)
        : CellLinkedList<AdaptiveSmoothingLength>(tentative_bounds, reference_grid_spacing, total_levels,
                                                  base_particles, sph_adaptation)
    {
    }

    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override
    {
        const UnsignedInt level = getMeshLevelWithTolerance(kernel_.CutOffRadius(h_ratio_[particle_index]));
        h_level_[particle_index] = static_cast<int>(level);
        const UnsignedInt linear_index = getMesh(level).LinearCellIndexFromPosition(particle_position);
        cell_index_lists_[linear_index].emplace_back(particle_index);
    }

    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override
    {
        const UnsignedInt level = getMeshLevelWithTolerance(kernel_.CutOffRadius(h_ratio_[particle_index]));
        const UnsignedInt linear_index = getMesh(level).LinearCellIndexFromPosition(particle_position);
        cell_data_lists_[linear_index].emplace_back(std::make_pair(particle_index, particle_position));
    }

  private:
    UnsignedInt getMeshLevelWithTolerance(Real particle_cutoff_radius)
    {
        for (UnsignedInt level = resolution_levels_; level != 0; --level)
        {
            const Real grid_spacing = getMesh(level - 1).GridSpacing();
            const Real tolerance = SMAX(SqrtEps, Real(1.0e-6) * grid_spacing);
            if (particle_cutoff_radius - grid_spacing < tolerance)
            {
                return level - 1;
            }
        }
        if (resolution_levels_ > 0)
        {
            const UnsignedInt coarsest_level = resolution_levels_ - 1;
            const Real coarsest_spacing = getMesh(coarsest_level).GridSpacing();
            const Real fallback_tolerance = SMAX(SqrtEps, Real(1.0e-4) * coarsest_spacing);
            if (particle_cutoff_radius - coarsest_spacing < fallback_tolerance)
            {
                return coarsest_level;
            }
        }
        std::cout << "\n Error: Team7TolerantAdaptiveCellLinkedList level searching out of bound!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
        return 999;
    }
};

/** Uniform / near-surface adaptive spacing with tolerant CLL (reload-safe). */
class AphiTeam7NativeAdaptiveNearSurface : public AdaptiveNearSurface
{
  public:
    template <typename... Args>
    explicit AphiTeam7NativeAdaptiveNearSurface(Args &&...args)
        : AdaptiveNearSurface(std::forward<Args>(args)...)
    {
    }

    UniquePtr<BaseCellLinkedList> createCellLinkedList(const BoundingBoxd &domain_bounds,
                                                       BaseParticles &base_particles) override
    {
        return makeUnique<Team7TolerantAdaptiveCellLinkedList>(
            domain_bounds, getKernel()->CutOffRadius(), LocalRefinementLevel() + 1, base_particles, *this);
    }
};

/** Air spacing vs distance to coil/plate union. */
class AphiTeam7NativeAdaptiveNearInnerSurface : public AphiTeam7NativeAdaptiveNearSurface
{
  public:
    AphiTeam7NativeAdaptiveNearInnerSurface(Real global_resolution, Real h_spacing_ratio, Real refinement_to_global,
                                            int local_refinement_level, Shape *inner_shape)
        : AphiTeam7NativeAdaptiveNearSurface(global_resolution, h_spacing_ratio, refinement_to_global,
                                             local_refinement_level),
          inner_shape_(inner_shape)
    {
    }
    Real getLocalSpacing(Shape &shape, const Vecd &position) override
    {
        (void)shape;
        const Real phi = std::fabs(inner_shape_->findSignedDistance(position));
        return smoothedSpacing(phi, spacing_ref_);
    }

  private:
    Shape *inner_shape_;
};

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_TEAM7_SHAPES_H
