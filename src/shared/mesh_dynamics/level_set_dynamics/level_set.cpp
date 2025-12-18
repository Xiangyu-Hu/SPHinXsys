#include "level_set.hpp"
#include "adaptation.h"
#include "base_kernel.h"
#include <type_traits>

namespace SPH
{
//=================================================================================================//
LevelSet::LevelSet(
    BoundingBoxd tentative_bounds, MeshWithGridDataPackagesType *coarse_data,
    Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : MeshWithGridDataPackages<4>(
          "LevelSet_" + shape.getName(), 1, tentative_bounds,
          coarse_data->getResolutionLevel(0).GridSpacing() * 0.5, 4, 2),
      shape_(shape), refinement_ratio_(refinement_ratio)
{
    Real data_spacing = coarse_data->getResolutionLevel(0).DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / data_spacing / refinement_ratio;
    Real smoothing_length = sph_adaptation.ReferenceSmoothingLength() / global_h_ratio;
    global_h_ratio_vec_.push_back(global_h_ratio);
    neighbor_method_set_.push_back(
        neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
            *sph_adaptation.getKernel(), smoothing_length, data_spacing));

    initializeLevel(0, coarse_data, coarse_data->getMeshes().size() - 1);
}
//=================================================================================================//
LevelSet::LevelSet(
    BoundingBoxd tentative_bounds, Real data_spacing,
    size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : MeshWithGridDataPackages<4>(
          "LevelSet_" + shape.getName(), total_levels, tentative_bounds,
          data_spacing * Real(4), 4, 2),
      shape_(shape), refinement_ratio_(refinement_ratio)
{
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / data_spacing / refinement_ratio;
    Real smoothing_length = sph_adaptation.ReferenceSmoothingLength() / global_h_ratio;
    global_h_ratio_vec_.push_back(global_h_ratio);
    neighbor_method_set_.push_back(
        neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
            *sph_adaptation.getKernel(), smoothing_length, data_spacing));

    initializeLevel(0, this, 0);
    for (size_t level = 1; level < resolution_levels_; ++level)
    {
        data_spacing *= 0.5;     // Halve the data spacing
        global_h_ratio *= 2;     // Double the ratio
        smoothing_length *= 0.5; // Halve the smoothing length
        global_h_ratio_vec_.push_back(global_h_ratio);
        neighbor_method_set_.push_back(
            neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
                *sph_adaptation.getKernel(), smoothing_length, data_spacing));

        initializeLevel(0, this, level - 1);
    }
}
//=================================================================================================//
void LevelSet::initializeLevel(
    UnsignedInt level, MeshWithGridDataPackagesType *coarse_data, UnsignedInt coarse_level)
{
    if (level == 0)
    {
        MeshAllDynamics<execution::ParallelPolicy, InitialCellTagging>
            initial_cell_tagging(*this, shape_);
        initial_cell_tagging.exec();
    }
    else
    {
        MeshAllDynamics<execution::ParallelPolicy, InitialCellTaggingFromCoarse>
            initial_cell_tagging_from_coarse(*this, level, *coarse_data, coarse_level, shape_);
        initial_cell_tagging_from_coarse.exec();
    }

    MeshAllDynamics<execution::ParallelPolicy, InnerCellTagging> tag_a_cell_is_inner_package(*this, level);
    tag_a_cell_is_inner_package.exec();
    organizeOccupiedPackages(level);
    PackageSort<execution::ParallelPolicy> pkg_sort(*this, level);
    pkg_sort.exec();

    /* All initializations in `FinishDataPackages` are achieved on CPU. */
    FinishDataPackages finish_data_packages(*this, level, shape_);
    finish_data_packages.exec();
}
//=================================================================================================//
size_t LevelSet::getCoarseLevel(Real h_ratio)
{
    for (size_t level = resolution_levels_; level != 0; --level)
        if (h_ratio > global_h_ratio_vec_[level - 1])
            return level - 1; // jump out the loop!

    std::cout << "\n Error: LevelSet level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void LevelSet::cleanInterface(UnsignedInt repeat_times)
{
    DynamicCast<RepeatTimes>(this, *clean_interface_keeper_.get())(repeat_times);
    clean_interface_keeper_->exec();
    sync_mesh_variables_to_probe_();
}
//=============================================================================================//
void LevelSet::correctTopology()
{
    correct_topology_keeper_->exec();
    sync_mesh_variables_to_probe_();
}
//=============================================================================================//
Real LevelSet::probeSignedDistance(const Vecd &position)
{
    return (*probe_signed_distance_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
Vecd LevelSet::probeNormalDirection(const Vecd &position)
{
    return (*probe_normal_direction_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
Vecd LevelSet::probeLevelSetGradient(const Vecd &position)
{
    return (*probe_level_set_gradient_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
size_t LevelSet::getProbeLevel(const Vecd &position)
{
    for (size_t level = resolution_levels_; level != 0; --level)
    {
        if (mesh_index_handler_set_[level - 1]->isWithinCorePackage(cell_pkg_index_, pkg_type_, position))
            return level - 1; // jump out of the loop!
    }
    return 0;
}
//=================================================================================================//
Real LevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    if (resolution_levels_ == 1)
    {
        return (*probe_kernel_integral_set_[0])(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Real coarse_level_value = (*probe_kernel_integral_set_[coarse_level])(position);
    Real fine_level_value = (*probe_kernel_integral_set_[coarse_level + 1])(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Vecd LevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    // std::cout << "probe kernel gradient integral" << std::endl;
    if (resolution_levels_ == 1)
    {
        return (*probe_kernel_gradient_integral_set_[0])(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Vecd coarse_level_value = (*probe_kernel_gradient_integral_set_[coarse_level])(position);
    Vecd fine_level_value = (*probe_kernel_gradient_integral_set_[coarse_level + 1])(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Matd LevelSet::probeKernelSecondGradientIntegral(const Vecd &position, Real h_ratio)
{
    if (resolution_levels_ == 1)
    {
        return (*probe_kernel_second_gradient_integral_set_[0])(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Matd coarse_level_value = (*probe_kernel_second_gradient_integral_set_[coarse_level])(position);
    Matd fine_level_value = (*probe_kernel_second_gradient_integral_set_[coarse_level + 1])(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
void LevelSet::writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence)
{
    sync_mesh_variables_to_write_();
    writeMeshFieldToPlt(partial_file_name, sequence);
}
//=============================================================================================//
} // namespace SPH
