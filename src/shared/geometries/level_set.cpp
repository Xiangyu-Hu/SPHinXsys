#include "level_set.h"

#include "adaptation.h"
#include "base_kernel.h"

namespace SPH
{
//=================================================================================================//
BaseLevelSet ::BaseLevelSet(Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), sph_adaptation_(sph_adaptation)
{
    if (!shape_.isValid())
    {
        std::cout << "\n BaseLevelSet Error: shape_ is invalid." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
}
//=================================================================================================//
LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size,
                   Shape &shape, SPHAdaptation &sph_adaptation)
    : MeshWithGridDataPackages<4>(tentative_bounds, data_spacing, buffer_size),
      BaseLevelSet(shape, sph_adaptation),
      global_h_ratio_(sph_adaptation.ReferenceSpacing() / data_spacing),
      phi_(*registerMeshVariable<Real>("Levelset")),
      near_interface_id_(*registerMeshVariable<int>("NearInterfaceID")),
      phi_gradient_(*registerMeshVariable<Vecd>("LevelsetGradient")),
      kernel_weight_(*registerMeshVariable<Real>("KernelWeight")),
      kernel_gradient_(*registerMeshVariable<Vecd>("KernelGradient")),
      kernel_(*sph_adaptation.getKernel()) {}
//=================================================================================================//
LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing,
                   Shape &shape, SPHAdaptation &sph_adaptation)
    : LevelSet(tentative_bounds, data_spacing, 4, shape, sph_adaptation)
{
    initialize_data_in_a_cell.exec();
    finish_data_packages.exec();
}
//=================================================================================================//
Vecd LevelSet::probeNormalDirection(const Vecd &position)
{
    return probe_normal_direction.exec(position);
}
//=================================================================================================//
Vecd LevelSet::probeLevelSetGradient(const Vecd &position)
{
    return probe_level_set_gradient.exec(position);
}
//=================================================================================================//
Real LevelSet::probeSignedDistance(const Vecd &position)
{
    return probe_signed_distance.exec(position);
}
//=================================================================================================//
Real LevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    return probe_kernel_integral.exec(position);
}
//=================================================================================================//
Vecd LevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    return probe_kernel_gradient_integral.exec(position);
}
//=================================================================================================//
void LevelSet::cleanInterface(Real small_shift_factor)
{
    mark_near_interface.setSmallShiftFactor(small_shift_factor);
    mark_near_interface.exec();
    redistance_interface.exec();
    reinitialize_level_set.exec();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();
}
//=============================================================================================//
void LevelSet::correctTopology(Real small_shift_factor)
{
    mark_near_interface.setSmallShiftFactor(small_shift_factor);
    mark_near_interface.exec();
    for (size_t i = 0; i != 10; ++i)
        diffuse_level_set_sign.exec();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();
}
//=================================================================================================//
bool LevelSet::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    Arrayi cell_pos = mesh_data_.CellIndexFromPosition(position);
    for (int i = 0; i != position.size(); ++i)
    {
        if (cell_pos[i] < 2)
            is_bounded = false;
        if (cell_pos[i] > (all_cells_[i] - 2))
            is_bounded = false;
    }
    return is_bounded;
}
//=============================================================================================//
bool LevelSet::isWithinCorePackage(Vecd position)
{
    Arrayi cell_index = mesh_data_.CellIndexFromPosition(position);
    return mesh_data_.isCoreDataPackage(cell_index);
}
//=============================================================================================//
RefinedLevelSet::RefinedLevelSet(BoundingBox tentative_bounds, LevelSet &coarse_level_set,
                                 Shape &shape, SPHAdaptation &sph_adaptation)
    : RefinedMesh(tentative_bounds, coarse_level_set, 4, shape, sph_adaptation)
{
    MeshAllDynamics<InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(mesh_data_, coarse_mesh_, shape_);
    initialize_data_in_a_cell_from_coarse.exec();

    finish_data_packages.exec();
}
//=============================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : MultilevelMesh<BaseLevelSet, LevelSet, RefinedLevelSet>(
          tentative_bounds, reference_data_spacing, total_levels, shape, sph_adaptation) {}
//=================================================================================================//
size_t MultilevelLevelSet::getCoarseLevel(Real h_ratio)
{
    for (size_t level = total_levels_; level != 0; --level)
        if (h_ratio > mesh_levels_[level - 1]->global_h_ratio_)
            return level - 1; // jump out the loop!

    std::cout << "\n Error: LevelSet level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void MultilevelLevelSet::cleanInterface(Real small_shift_factor)
{
    mesh_levels_.back()->cleanInterface(small_shift_factor);
}
//=============================================================================================//
void MultilevelLevelSet::correctTopology(Real small_shift_factor)
{
    mesh_levels_.back()->correctTopology(small_shift_factor);
}
//=============================================================================================//
Real MultilevelLevelSet::probeSignedDistance(const Vecd &position)
{
    return mesh_levels_[getProbeLevel(position)]->probeSignedDistance(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeNormalDirection(const Vecd &position)
{
    return mesh_levels_[getProbeLevel(position)]->probeNormalDirection(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeLevelSetGradient(const Vecd &position)
{
    return mesh_levels_[getProbeLevel(position)]->probeLevelSetGradient(position);
}
//=============================================================================================//
size_t MultilevelLevelSet::getProbeLevel(const Vecd &position)
{
    for (size_t level = total_levels_; level != 0; --level)
        if (mesh_levels_[level - 1]->isWithinCorePackage(position))
            return level - 1; // jump out of the loop!
    return 0;
}
//=================================================================================================//
Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (mesh_levels_[coarse_level + 1]->global_h_ratio_ - h_ratio) /
                 (mesh_levels_[coarse_level + 1]->global_h_ratio_ - mesh_levels_[coarse_level]->global_h_ratio_);
    Real coarse_level_value = mesh_levels_[coarse_level]->probeKernelIntegral(position);
    Real fine_level_value = mesh_levels_[coarse_level + 1]->probeKernelIntegral(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (mesh_levels_[coarse_level + 1]->global_h_ratio_ - h_ratio) /
                 (mesh_levels_[coarse_level + 1]->global_h_ratio_ - mesh_levels_[coarse_level]->global_h_ratio_);
    Vecd coarse_level_value = mesh_levels_[coarse_level]->probeKernelGradientIntegral(position);
    Vecd fine_level_value = mesh_levels_[coarse_level + 1]->probeKernelGradientIntegral(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
bool MultilevelLevelSet::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    for (size_t l = 0; l != total_levels_; ++l)
    {
        if (!mesh_levels_[l]->probeIsWithinMeshBound(position))
        {
            is_bounded = false;
            break;
        };
    }
    return is_bounded;
}
//=============================================================================================//
} // namespace SPH
