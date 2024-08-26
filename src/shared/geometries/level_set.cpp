#include "level_set.h"

#include "adaptation.h"
#include "base_kernel.h"

namespace SPH
{
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), kernel_(*sph_adaptation.getKernel()), shape_(shape), total_levels_(total_levels)
{
    mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));

    mesh_data_set_[0]->registerMeshVariable<Real>("Levelset");
    mesh_data_set_[0]->registerMeshVariable<int>("NearInterfaceID");
    mesh_data_set_[0]->registerMeshVariable<Vecd>("LevelsetGradient");
    mesh_data_set_[0]->registerMeshVariable<Real>("KernelWeight");
    mesh_data_set_[0]->registerMeshVariable<Vecd>("KernelGradient");

    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);
    MeshAllDynamics<InitializeDataInACell> initialize_data_in_a_cell{*mesh_data_set_[0], shape_};
    FinishDataPackages finish_data_packages{*mesh_data_set_[0], shape_, kernel_, global_h_ratio};
    initialize_data_in_a_cell.exec();
    finish_data_packages.exec();

    for (size_t level = 1; level != total_levels_; ++level)
    {
        reference_data_spacing *= 0.5;
        global_h_ratio *= 2;
        global_h_ratio_vec_.push_back(global_h_ratio);
        /** all mesh levels aligned at the lower bound of tentative_bounds */
        mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));
        mesh_data_set_[level]->registerMeshVariable<Real>("Levelset");
        mesh_data_set_[level]->registerMeshVariable<int>("NearInterfaceID");
        mesh_data_set_[level]->registerMeshVariable<Vecd>("LevelsetGradient");
        mesh_data_set_[level]->registerMeshVariable<Real>("KernelWeight");
        mesh_data_set_[level]->registerMeshVariable<Vecd>("KernelGradient");
        MeshAllDynamics<InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(*mesh_data_set_[level], *mesh_data_set_[level - 1], shape_);
        FinishDataPackages finish_data_packages{*mesh_data_set_[level], shape_, kernel_, global_h_ratio};
        initialize_data_in_a_cell_from_coarse.exec();
        finish_data_packages.exec();
    }
}
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
    BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data, Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), kernel_(*sph_adaptation.getKernel()), shape_(shape)
{
    total_levels_ = 1;
    Real reference_data_spacing = coarse_data->DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));
    mesh_data_set_[0]->registerMeshVariable<Real>("Levelset");
    mesh_data_set_[0]->registerMeshVariable<int>("NearInterfaceID");
    mesh_data_set_[0]->registerMeshVariable<Vecd>("LevelsetGradient");
    mesh_data_set_[0]->registerMeshVariable<Real>("KernelWeight");
    mesh_data_set_[0]->registerMeshVariable<Vecd>("KernelGradient");

    MeshAllDynamics<InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(*mesh_data_set_[0], *coarse_data, shape_);
    FinishDataPackages finish_data_packages{*mesh_data_set_[0], shape_, kernel_, global_h_ratio};
    initialize_data_in_a_cell_from_coarse.exec();
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
    CleanInterface clean_interface{*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back()};
    clean_interface.exec(small_shift_factor);
}
//=============================================================================================//
void MultilevelLevelSet::correctTopology(Real small_shift_factor)
{
    CorrectTopology correct_topology{*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back()};
    correct_topology.exec(small_shift_factor);
}
//=============================================================================================//
Real MultilevelLevelSet::probeSignedDistance(const Vecd &position)
{
    MeshCalculateDynamics<Real, ProbeSignedDistance> probe_signed_distance{*mesh_data_set_[getProbeLevel(position)]};
    return probe_signed_distance.exec(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeNormalDirection(const Vecd &position)
{
    MeshCalculateDynamics<Vecd, ProbeNormalDirection> probe_normal_direction{*mesh_data_set_[getProbeLevel(position)]};
    return probe_normal_direction.exec(position);
}
//=============================================================================================//
Vecd MultilevelLevelSet::probeLevelSetGradient(const Vecd &position)
{
    MeshCalculateDynamics<Vecd, ProbeLevelSetGradient> probe_levelset_gradient{*mesh_data_set_[getProbeLevel(position)]};
    return probe_levelset_gradient.exec(position);
}
//=============================================================================================//
size_t MultilevelLevelSet::getProbeLevel(const Vecd &position)
{
    for (size_t level = total_levels_; level != 0; --level){
        MeshCalculateDynamics<bool, IsWithinCorePackage> is_within_core_package{*mesh_data_set_[level - 1]};
        if(is_within_core_package.exec(position))
            return level - 1; // jump out of the loop!
    }
    return 0;
}
//=================================================================================================//
Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position)
{
    MeshCalculateDynamics<Real, ProbeKernelIntegral> refine_probe{*mesh_data_set_[0]};
    return refine_probe.exec(position);
}
//=================================================================================================//
Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        MeshCalculateDynamics<Real, ProbeKernelIntegral> refine_probe{*mesh_data_set_[0]};
        return refine_probe.exec(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    MeshCalculateDynamics<Real, ProbeKernelIntegral> coarse_probe{*mesh_data_set_[coarse_level]};
    MeshCalculateDynamics<Real, ProbeKernelIntegral> fine_probe{*mesh_data_set_[coarse_level + 1]};
    Real coarse_level_value = coarse_probe.exec(position);
    Real fine_level_value = fine_probe.exec(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position)
{
    MeshCalculateDynamics<Vecd, ProbeKernelGradientIntegral> refine_probe{*mesh_data_set_[0]};
    return refine_probe.exec(position);
}
//=================================================================================================//
Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        MeshCalculateDynamics<Vecd, ProbeKernelGradientIntegral> refine_probe{*mesh_data_set_[0]};
        return refine_probe.exec(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    MeshCalculateDynamics<Vecd, ProbeKernelGradientIntegral> coarse_probe{*mesh_data_set_[coarse_level]};
    MeshCalculateDynamics<Vecd, ProbeKernelGradientIntegral> fine_probe{*mesh_data_set_[coarse_level + 1]};
    Vecd coarse_level_value = coarse_probe.exec(position);
    Vecd fine_level_value = fine_probe.exec(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
bool MultilevelLevelSet::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    for (size_t l = 0; l != total_levels_; ++l)
    {
        MeshCalculateDynamics<bool, ProbeIsWithinMeshBound> probe_is_within_mesh_bound{*mesh_data_set_[l]};
        if (!probe_is_within_mesh_bound.exec(position))
        {
            is_bounded = false;
            break;
        };
    }
    return is_bounded;
}
//=============================================================================================//
} // namespace SPH
