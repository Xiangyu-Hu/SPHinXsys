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
Real BaseLevelSet::CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing)
{
    Real squared_norm_inv = 1.0 / (phi_gradient.squaredNorm() + TinyReal);
    Real volume_fraction(0);
    for (size_t i = 0; i != Dimensions; ++i)
    {
        volume_fraction += phi_gradient[i] * phi_gradient[i] * squared_norm_inv *
                           Heaviside(phi / (ABS(phi_gradient[i]) + TinyReal), 0.5 * data_spacing);
    }
    return volume_fraction;
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
    finishDataPackages();
}
//=================================================================================================//
void LevelSet::finishDataPackages()
{
    tag_a_cell_is_inner_package.exec();

    mesh_data_.organizeOccupiedPackages();
    initialize_index_mesh.exec();
    initialize_cell_neighborhood.exec();
    resizeMeshVariableData();


    Real far_field_distance = grid_spacing_ * (Real)buffer_width_;
    initializeDataForSingularPackage(0, -far_field_distance);
    initializeDataForSingularPackage(1, far_field_distance);

    initialize_basic_data_for_a_package.exec();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();
}
//=================================================================================================//
Vecd LevelSet::probeNormalDirection(const Vecd &position)
{
    Vecd probed_value = probeLevelSetGradient(position);

    Real threshold = 1.0e-2 * data_spacing_;
    while (probed_value.norm() < threshold)
    {
        Vecd jittered = position; // jittering
        for (int l = 0; l != position.size(); ++l)
            jittered[l] += rand_uniform(-0.5, 0.5) * 0.5 * data_spacing_;
        probed_value = probeLevelSetGradient(jittered);
    }
    return probed_value.normalized();
}
//=================================================================================================//
Vecd LevelSet::probeLevelSetGradient(const Vecd &position)
{
    return probeMesh(phi_gradient_, position);
}
//=================================================================================================//
Real LevelSet::probeSignedDistance(const Vecd &position)
{
    return probeMesh(phi_, position);
}
//=================================================================================================//
Real LevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    return probeMesh(kernel_weight_, position);
}
//=================================================================================================//
Vecd LevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    return probeMesh(kernel_gradient_, position);
}
//=================================================================================================//
void LevelSet::redistanceInterface()
{
    // redistance_interface.exec();
    package_parallel_for(
        [&](size_t package_index) {
            std::pair<Arrayi, int> &metadata = meta_data_cell_[package_index];
            if (metadata.second == 1)
            {
                redistanceInterfaceForAPackage(PackageIndexFromCellIndex(metadata.first));
            }
        });
}
//=================================================================================================//
void LevelSet::cleanInterface(Real small_shift_factor)
{
    markNearInterface(small_shift_factor);
    redistanceInterface();
    reinitializeLevelSet();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();
}
//=============================================================================================//
void LevelSet::correctTopology(Real small_shift_factor)
{
    markNearInterface(small_shift_factor);
    for (size_t i = 0; i != 10; ++i)
        diffuseLevelSetSign();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();
}
//=================================================================================================//
bool LevelSet::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    Arrayi cell_pos = CellIndexFromPosition(position);
    for (int i = 0; i != position.size(); ++i)
    {
        if (cell_pos[i] < 2)
            is_bounded = false;
        if (cell_pos[i] > (all_cells_[i] - 2))
            is_bounded = false;
    }
    return is_bounded;
}
//=================================================================================================//
Real LevelSet::upwindDifference(Real sign, Real df_p, Real df_n)
{
    if (sign * df_p >= 0.0 && sign * df_n >= 0.0)
        return df_n;
    if (sign * df_p <= 0.0 && sign * df_n <= 0.0)
        return df_p;
    if (sign * df_p > 0.0 && sign * df_n < 0.0)
        return 0.0;

    Real df = df_p;
    if (sign * df_p < 0.0 && sign * df_n > 0.0)
    {
        Real ss = sign * (fabs(df_p) - fabs(df_n)) / (df_p - df_n);
        if (ss > 0.0)
            df = df_n;
    }

    return df;
}
//=============================================================================================//
void RefinedLevelSet::initializeDataInACellFromCoarse(const Arrayi &cell_index)
{
    Vecd cell_position = CellPositionFromIndex(cell_index);
    size_t package_index = coarse_mesh_.probeSignedDistance(cell_position) < 0.0 ? 0 : 1;
    assignDataPackageIndex(cell_index, package_index);
    if (coarse_mesh_.isWithinCorePackage(cell_position))
    {
        Real signed_distance = shape_.findSignedDistance(cell_position);
        Vecd normal_direction = shape_.findNormalDirection(cell_position);
        Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
        if (measure < grid_spacing_)
        {
            std::pair<size_t, int> occupied;
            occupied.first = cell_index[0] * all_cells_[1] + cell_index[1];
            occupied.second = 1;

            mesh_data_.assignDataPackageIndex(cell_index, 2);
            mesh_data_.registerOccupied(occupied);
        }
    }
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
