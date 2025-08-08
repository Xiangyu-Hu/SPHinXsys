#include "level_set.hpp"
#include "adaptation.h"
#include "base_kernel.h"
#include <type_traits>

namespace SPH
{
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, MeshWithGridDataPackagesType *coarse_data,
    Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(1)
{
    Real reference_data_spacing = coarse_data->DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing / refinement_ratio;
    kernel_ = makeUnique<SingularVariable<KernelTabulatedCK>>(
        "levelset_kernel", KernelTabulatedCK(*sph_adaptation.getKernel()));
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(reference_data_spacing, tentative_bounds, coarse_data);
}
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, Real reference_data_spacing,
    size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio)
    : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(total_levels)
{
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing / refinement_ratio;
    global_h_ratio_vec_.push_back(global_h_ratio);
    kernel_ = makeUnique<SingularVariable<KernelTabulatedCK>>(
        "levelset_kernel", KernelTabulatedCK(*sph_adaptation.getKernel()));

    initializeLevel(reference_data_spacing, tentative_bounds);
    for (size_t level = 1; level < total_levels_; ++level)
    {
        reference_data_spacing *= 0.5; // Halve the data spacing
        global_h_ratio *= 2;           // Double the ratio
        global_h_ratio_vec_.push_back(global_h_ratio);

        initializeLevel(reference_data_spacing, tentative_bounds, mesh_data_set_[level - 1]);
    }
}
//=================================================================================================//
void MultilevelLevelSet::initializeLevel(
    Real reference_data_spacing, BoundingBox tentative_bounds, MeshWithGridDataPackagesType *coarse_data)
{
    MeshWithGridDataPackagesType *mesh_data =
        mesh_data_ptr_vector_keeper_
            .template createPtr<MeshWithGridDataPackagesType>(
                tentative_bounds, reference_data_spacing, 4);
    mesh_data_set_.push_back(mesh_data);

    if (coarse_data == nullptr)
    {
        MeshAllDynamics<execution::ParallelPolicy, InitialCellTagging> initial_cell_tagging(*mesh_data, shape_);
        initial_cell_tagging.exec();
    }
    else
    {
        MeshAllDynamics<execution::ParallelPolicy, InitialCellTaggingFromCoarse>
            initial_cell_tagging_from_coarse(*mesh_data, *coarse_data, shape_);
        initial_cell_tagging_from_coarse.exec();
    }

    MeshAllDynamics<execution::ParallelPolicy, InnerCellTagging> tag_a_cell_is_inner_package(*mesh_data);
    tag_a_cell_is_inner_package.exec();
    mesh_data->organizeOccupiedPackages();

    mesh_data->registerMeshVariable<Real>("LevelSet");
    mesh_data->registerMeshVariable<int>("NearInterfaceID");
    mesh_data->registerMeshVariable<Vecd>("LevelSetGradient");
    mesh_data->addMeshVariableToWrite<Real>("LevelSet");
    mesh_data->addMeshVariableToWrite<int>("NearInterfaceID");
    mesh_data->addMeshVariableToWrite<Vecd>("LevelSetGradient");

    /* All initializations in `FinishDataPackages` are achieved on CPU. */
    FinishDataPackages finish_data_packages(*mesh_data, shape_);
    finish_data_packages.exec();
}
//=================================================================================================//
size_t MultilevelLevelSet::getCoarseLevel(Real h_ratio)
{
    for (size_t level = total_levels_; level != 0; --level)
        if (h_ratio > global_h_ratio_vec_[level - 1])
            return level - 1; // jump out the loop!

    std::cout << "\n Error: LevelSet level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void MultilevelLevelSet::cleanInterface(Real small_shift_factor)
{
    clean_interface_keeper_->exec(small_shift_factor);
}
//=============================================================================================//
void MultilevelLevelSet::correctTopology(Real small_shift_factor)
{
    correct_topology_keeper_->exec(small_shift_factor);
}
//=============================================================================================//
Real MultilevelLevelSet::probeSignedDistance(const Vecd &position)
{
    return (*probe_signed_distance_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeNormalDirection(const Vecd &position)
{
    return (*probe_normal_direction_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeLevelSetGradient(const Vecd &position)
{
    return (*probe_level_set_gradient_set_[getProbeLevel(position)])(position);
}
//=============================================================================================//
size_t MultilevelLevelSet::getProbeLevel(const Vecd &position)
{
    for (size_t level = total_levels_; level != 0; --level)
    {
        if (mesh_data_set_[level - 1]->isWithinCorePackage(cell_package_index_set_[level - 1],
                                                           meta_data_cell_set_[level - 1],
                                                           position))
            return level - 1; // jump out of the loop!
    }
    return 0;
}
//=================================================================================================//
Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    // std::cout << "probe kernel integral" << std::endl;
    if (mesh_data_set_.size() == 1)
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
Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    // std::cout << "probe kernel gradient integral" << std::endl;
    if (mesh_data_set_.size() == 1)
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
Matd MultilevelLevelSet::probeKernelSecondGradientIntegral(const Vecd &position, Real h_ratio)
{
    if (mesh_data_set_.size() == 1)
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
bool MultilevelLevelSet::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    for (size_t l = 0; l != total_levels_; ++l)
    {
        ProbeIsWithinMeshBound probe_is_within_mesh_bound{*mesh_data_set_[l]};
        if (!probe_is_within_mesh_bound.update(position))
        {
            is_bounded = false;
            break;
        };
    }
    return is_bounded;
}
//=================================================================================================//
void MultilevelLevelSet::writeMeshFieldToPlt(const std::string &partial_file_name)
{
    sync_mesh_variable_data_();
    resetProbes();
    for (size_t l = 0; l != total_levels_; ++l)
    {
        std::string full_file_name = partial_file_name + "_" + std::to_string(l) + ".dat";
        std::ofstream out_file(full_file_name.c_str(), std::ios::app);
        mesh_data_set_[l]->writeMeshVariableToPlt(out_file);
        out_file.close();
    }
}
//=================================================================================================//
void MultilevelLevelSet::writeBackgroundMeshToPlt(const std::string &partial_file_name)
{
    for (size_t l = 0; l != total_levels_; ++l)
    {
        std::string full_file_name = partial_file_name + "_" + std::to_string(l) + ".dat";
        std::ofstream out_file(full_file_name.c_str(), std::ios::app);
        mesh_data_set_[l]->writeDiscreteVariableToPlt(out_file);
        out_file.close();
    }
}
//=============================================================================================//
} // namespace SPH
