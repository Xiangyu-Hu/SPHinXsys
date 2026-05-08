#ifndef EM_ADAPTIVE_CELL_LINKED_LIST_H
#define EM_ADAPTIVE_CELL_LINKED_LIST_H

#include "sphinxsys.h"

namespace SPH
{
namespace extra_em_mesh
{
class TolerantAdaptiveCellLinkedList : public CellLinkedList<AdaptiveSmoothingLength>
{
  public:
    TolerantAdaptiveCellLinkedList(BoundingBoxd tentative_bounds,
                                   Real reference_grid_spacing,
                                   UnsignedInt total_levels,
                                   BaseParticles &base_particles,
                                   SPHAdaptation &sph_adaptation)
        : CellLinkedList<AdaptiveSmoothingLength>(tentative_bounds,
                                                  reference_grid_spacing,
                                                  total_levels,
                                                  base_particles,
                                                  sph_adaptation) {}

    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override
    {
        UnsignedInt level = getMeshLevelWithTolerance(kernel_.CutOffRadius(h_ratio_[particle_index]));
        h_level_[particle_index] = static_cast<int>(level);
        UnsignedInt linear_index = getMesh(level).LinearCellIndexFromPosition(particle_position);
        cell_index_lists_[linear_index].emplace_back(particle_index);
    }

    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override
    {
        UnsignedInt level = getMeshLevelWithTolerance(kernel_.CutOffRadius(h_ratio_[particle_index]));
        UnsignedInt linear_index = getMesh(level).LinearCellIndexFromPosition(particle_position);
        cell_data_lists_[linear_index].emplace_back(std::make_pair(particle_index, particle_position));
    }

  private:
    UnsignedInt getMeshLevelWithTolerance(Real particle_cutoff_radius)
    {
        for (UnsignedInt level = resolution_levels_; level != 0; --level)
        {
            Real grid_spacing = getMesh(level - 1).GridSpacing();
            Real tolerance = SMAX(SqrtEps, Real(1.0e-6) * grid_spacing);
            if (particle_cutoff_radius - grid_spacing < tolerance)
            {
                return level - 1;
            }
        }

        if (resolution_levels_ > 0)
        {
            UnsignedInt coarsest_level = resolution_levels_ - 1;
            Real coarsest_spacing = getMesh(coarsest_level).GridSpacing();
            Real fallback_tolerance = SMAX(SqrtEps, Real(1.0e-4) * coarsest_spacing);
            if (particle_cutoff_radius - coarsest_spacing < fallback_tolerance)
            {
                return coarsest_level;
            }
        }

        std::cout << "\n Error: TolerantAdaptiveCellLinkedList level searching out of bound!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
        return 999;
    }
};

template <class BaseAdaptation>
class AdaptiveWithTolerantCellLinkedList : public BaseAdaptation
{
  public:
    typedef typename BaseAdaptation::CellLinkedListIdentifier CellLinkedListIdentifier;
    typedef typename BaseAdaptation::SmoothingLengthRatioType SmoothingLengthRatioType;

    template <typename... Args>
    explicit AdaptiveWithTolerantCellLinkedList(Args &&...args)
        : BaseAdaptation(std::forward<Args>(args)...) {}

    UniquePtr<BaseCellLinkedList> createCellLinkedList(const BoundingBoxd &domain_bounds,
                                                       BaseParticles &base_particles) override
    {
        return makeUnique<TolerantAdaptiveCellLinkedList>(
            domain_bounds,
            this->getKernel()->CutOffRadius(),
            this->LocalRefinementLevel() + 1,
            base_particles,
            *this);
    }
};
} // namespace extra_em_mesh
} // namespace SPH

#endif // EM_ADAPTIVE_CELL_LINKED_LIST_H
