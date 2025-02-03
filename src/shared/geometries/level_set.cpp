#include "level_set.h"

#include "adaptation.h"
#include "base_kernel.h"

namespace SPH
{
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : MultilevelLevelSet(par, tentative_bounds, reference_data_spacing, total_levels, shape, sph_adaptation){};
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data, Shape &shape, SPHAdaptation &sph_adaptation)
    : MultilevelLevelSet(par, tentative_bounds, coarse_data, shape, sph_adaptation){};
//=================================================================================================//
template <class ExecutionPolicy>
MultilevelLevelSet::MultilevelLevelSet(
    const ExecutionPolicy &ex_policy, BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : MultilevelLevelSet(ex_policy, tentative_bounds, coarse_data->DataSpacing() * 0.5, 1, shape, sph_adaptation){};
//=================================================================================================//
template <class ExecutionPolicy>
MultilevelLevelSet::MultilevelLevelSet(
    const ExecutionPolicy &ex_policy, BoundingBox tentative_bounds, Real reference_data_spacing,
    size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), kernel_(*sph_adaptation.getKernel()), shape_(shape), total_levels_(total_levels)
{
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(ex_policy, 0, reference_data_spacing, global_h_ratio, tentative_bounds);
    cell_package_index_set_.push_back(mesh_data_set_[0]->cell_package_index_.DelegatedDataField(ex_policy));
    meta_data_cell_set_.push_back(mesh_data_set_[0]->meta_data_cell_.DelegatedDataField(ex_policy));

    for (size_t level = 1; level < total_levels_; ++level) {
        reference_data_spacing *= 0.5;  // Halve the data spacing
        global_h_ratio *= 2;            // Double the ratio
        global_h_ratio_vec_.push_back(global_h_ratio);

        initializeLevel(ex_policy, level, reference_data_spacing, global_h_ratio, tentative_bounds);
        cell_package_index_set_.push_back(mesh_data_set_[level]->cell_package_index_.DelegatedDataField(ex_policy));
        meta_data_cell_set_.push_back(mesh_data_set_[level]->meta_data_cell_.DelegatedDataField(ex_policy));
    }

    clean_interface = makeUnique<CleanInterface>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
    correct_topology = makeUnique<CorrectTopology>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::initializeLevel(const ExecutionPolicy &ex_policy, size_t level,
                                         Real reference_data_spacing, Real global_h_ratio,
                                         BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data)
{
    mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));

    RegisterMeshVariable().exec(mesh_data_set_[level]);

    if (coarse_data == nullptr) {
        MeshAllDynamicsCK<execution::ParallelPolicy, InitializeDataInACell> initialize_data_in_a_cell(*mesh_data_set_[level], shape_);
        initialize_data_in_a_cell.exec();
    } else {
        MeshAllDynamicsCK<execution::ParallelPolicy, InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(*mesh_data_set_[level], *coarse_data, shape_);
        initialize_data_in_a_cell_from_coarse.exec();
    }

    /* All initializations in `FinishDataPackages` are achieved on CPU. */
    FinishDataPackages finish_data_packages(*mesh_data_set_[level], shape_);
    MeshInnerDynamicsCK<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{*mesh_data_set_[level]};
    MeshInnerDynamicsCK<execution::ParallelPolicy, UpdateKernelIntegrals> update_kernel_integrals{*mesh_data_set_[level], kernel_, global_h_ratio};
    finish_data_packages.exec();
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();

    registerProbes(ex_policy, level);
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::registerProbes(const ExecutionPolicy &ex_policy, size_t level)
{
    probe_signed_distance_set_.push_back(
        probe_signed_distance_vector_keeper_
            .template createPtr<ProbeSignedDistance>(ex_policy, mesh_data_set_[level]));
    probe_normal_direction_set_.push_back(
        probe_normal_direction_vector_keeper_
            .template createPtr<ProbeNormalDirection>(ex_policy, mesh_data_set_[level]));
    probe_level_set_gradient_set_.push_back(
        probe_level_set_gradient_vector_keeper_
            .template createPtr<ProbeLevelSetGradient>(ex_policy, mesh_data_set_[level]));
    probe_kernel_integral_set_.push_back(
        probe_kernel_integral_vector_keeper_
            .template createPtr<ProbeKernelIntegral>(ex_policy, mesh_data_set_[level]));
    probe_kernel_gradient_integral_set_.push_back(
        probe_kernel_gradient_integral_vector_keeper_
            .template createPtr<ProbeKernelGradientIntegral>(ex_policy, mesh_data_set_[level]));
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
    clean_interface->exec(small_shift_factor);
}
//=============================================================================================//
void MultilevelLevelSet::correctTopology(Real small_shift_factor)
{
    correct_topology->exec(small_shift_factor);
}
//=============================================================================================//
Real MultilevelLevelSet::probeSignedDistance(const Vecd &position)
{
    return probe_signed_distance_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeNormalDirection(const Vecd &position)
{
    return probe_normal_direction_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeLevelSetGradient(const Vecd &position)
{
    return probe_level_set_gradient_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
size_t MultilevelLevelSet::getProbeLevel(const Vecd &position)
{
    for (size_t level = total_levels_; level != 0; --level){
        if(mesh_data_set_[level - 1]->isWithinCorePackage(cell_package_index_set_[level - 1],
                                                          meta_data_cell_set_[level - 1],
                                                          position))
            return level - 1; // jump out of the loop!
    }
    return 0;
}
//=================================================================================================//
Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Real coarse_level_value = probe_kernel_integral_set_[coarse_level]->update(position);
    Real fine_level_value = probe_kernel_integral_set_[coarse_level + 1]->update(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_gradient_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Vecd coarse_level_value = probe_kernel_gradient_integral_set_[coarse_level]->update(position);
    Vecd fine_level_value = probe_kernel_gradient_integral_set_[coarse_level + 1]->update(position);

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
//=============================================================================================//
} // namespace SPH
