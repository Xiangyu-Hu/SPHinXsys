#include "level_set.hpp"
#include "adaptation.h"
#include "base_kernel.h"
#include <type_traits>

namespace SPH
{
//=================================================================================================//
LevelSet::LevelSet(
    BoundingBoxd tentative_bounds, SparseMeshField<4> *coarse_data,
    Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : SparseMeshField<4>(
          "LevelSet_" + shape.getName(), 1, tentative_bounds,
          coarse_data->getFinestMesh().GridSpacing() * 0.5, 4, 2),
      shape_(shape), refinement_ratio_(refinement_ratio), ca_global_h_ratio_(nullptr)
{
    StdVec<Real> global_h_ratio_vec;
    Real data_spacing = coarse_data->getFinestMesh().DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / data_spacing / refinement_ratio;
    Real smoothing_length = sph_adaptation.ReferenceSmoothingLength() / global_h_ratio;
    global_h_ratio_vec.push_back(global_h_ratio);
    neighbor_method_set_.push_back(
        neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
            *sph_adaptation.getKernel(), smoothing_length, data_spacing));

    initializeLevel(0, coarse_data, coarse_data->ResolutionLevels() - 1);
    ca_global_h_ratio_ = createUniqueEnity<Real, ConstantArray>("GlobalHRatio", global_h_ratio_vec);
}
//=================================================================================================//
LevelSet::LevelSet(
    BoundingBoxd tentative_bounds, Real data_spacing,
    size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : SparseMeshField<4>(
          "LevelSet_" + shape.getName(), total_levels, tentative_bounds,
          data_spacing * Real(4), 4, 2),
      shape_(shape), refinement_ratio_(refinement_ratio), ca_global_h_ratio_(nullptr)
{
    StdVec<Real> global_h_ratio_vec;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / data_spacing / refinement_ratio;
    Real smoothing_length = sph_adaptation.ReferenceSmoothingLength() / global_h_ratio;
    global_h_ratio_vec.push_back(global_h_ratio);
    neighbor_method_set_.push_back(
        neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
            *sph_adaptation.getKernel(), smoothing_length, data_spacing));

    initializeLevel(0);
    for (size_t level = 1; level < resolution_levels_; ++level)
    {
        data_spacing *= 0.5;     // Halve the data spacing
        global_h_ratio *= 2;     // Double the ratio
        smoothing_length *= 0.5; // Halve the smoothing length
        global_h_ratio_vec.push_back(global_h_ratio);
        neighbor_method_set_.push_back(
            neighbor_method_keeper_.template createPtr<NeighborMethod<SPHAdaptation, SPHAdaptation>>(
                *sph_adaptation.getKernel(), smoothing_length, data_spacing));

        initializeLevel(level, this, level - 1);
    }
    ca_global_h_ratio_ = createUniqueEnity<Real, ConstantArray>("GlobalHRatioi", global_h_ratio_vec);
}
//=================================================================================================//
void LevelSet::initializeLevel(
    UnsignedInt level, SparseMeshField<4> *coarse_data, UnsignedInt coarse_level)
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

    /* All initializations in `FinishPackageDatas` are achieved on CPU. */
    FinishPackageDatas finish_data_packages(*this, level, shape_);
    finish_data_packages.exec();
}
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
    return (*probe_signed_distance_)(position);
}
//=============================================================================================//
Vecd LevelSet::probeNormalDirection(const Vecd &position)
{
    return (*probe_level_set_gradient_)(position).normalized();
}
//=============================================================================================//
Vecd LevelSet::probeLevelSetGradient(const Vecd &position)
{
    return (*probe_level_set_gradient_)(position);
}
//=================================================================================================//
Real LevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    return (*probe_kernel_integral_)(position, h_ratio);
}
//=================================================================================================//
Vecd LevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    return (*probe_kernel_gradient_integral_)(position, h_ratio);
}
//=================================================================================================//
Matd LevelSet::probeKernelSecondGradientIntegral(const Vecd &position, Real h_ratio)
{
    return (*probe_kernel_second_gradient_integral_)(position, h_ratio);
}
//=================================================================================================//
void LevelSet::writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence)
{
    sync_mesh_variables_to_write_();
    SparseMeshField<4>::writeMeshFieldToPlt(partial_file_name, sequence);
}
//=============================================================================================//
} // namespace SPH
