#include "mesh_local_dynamics.hpp"

namespace SPH
{
//=================================================================================================//
BaseMeshLocalDynamics::BaseMeshLocalDynamics(MeshWithGridDataPackagesType &data_mesh)
    : data_mesh_(data_mesh),
      all_cells_(data_mesh.AllCells()),
      grid_spacing_(data_mesh.GridSpacing()),
      data_spacing_(data_mesh.DataSpacing()),
      num_singular_pkgs_(data_mesh.NumSingularPackages()) {}
//=================================================================================================//
InitialCellTagging::InitialCellTagging(MeshWithGridDataPackagesType &data_mesh, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh), shape_(shape),
      bmv_cell_contain_id_(*data_mesh.registerBKGMeshVariable<int>(
          "CellContainID", // default value is 2, indicating not near interface
          [&](UnsignedInt index)
          { return 2; })) // the value is not touched in this class
{
    data_mesh.addBKGMeshVariableToWrite<int>("CellContainID");
}
//=================================================================================================//
void InitialCellTagging::UpdateKernel::update(const Arrayi &cell_index)
{
    data_mesh_->assignDataPackageIndex(cell_index, 1); // outside far field by default
    UnsignedInt index_1d = data_mesh_->LinearCellIndex(cell_index);
    cell_contain_id_[index_1d] = 2; // default value is 2, indicating not near interface
    Vecd cell_position = data_mesh_->CellPositionFromIndex(cell_index);
    Real signed_distance = shape_->findSignedDistance(cell_position);
    Vecd normal_direction = shape_->findNormalDirection(cell_position);
    Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
    if (measure < grid_spacing_)
    {
        UnsignedInt sort_index = base_dynamics->SortIndexFromCellIndex(cell_index);
        data_mesh_->assignDataPackageIndex(cell_index, 2);
        data_mesh_->registerOccupied(sort_index, 1);
    }
}
//=================================================================================================//
InitialCellTaggingFromCoarse::InitialCellTaggingFromCoarse(
    MeshWithGridDataPackagesType &data_mesh,
    MeshWithGridDataPackagesType &coarse_mesh, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh),
      coarse_mesh_(coarse_mesh), shape_(shape),
      far_field_distance_(data_mesh.GridSpacing() * (Real)data_mesh.BufferWidth()),
      bmv_cell_contain_id_(*data_mesh.registerBKGMeshVariable<int>(
          "CellContainID", // default value is 2, indicating not near interface
          [&](UnsignedInt index)
          { return 2; })),
      bmv_cell_pkg_index__coarse_(coarse_mesh.getCellPackageIndex()),
      dv_pkg_cell_info_coarse_(coarse_mesh.dvPkgCellInfo())
{
    data_mesh.addBKGMeshVariableToWrite<int>("CellContainID");
}
//=================================================================================================//
void InitialCellTaggingFromCoarse::UpdateKernel::update(const Arrayi &cell_index)
{
    Vecd cell_position = data_mesh_->CellPositionFromIndex(cell_index);
    Real phi = probe_coarse_phi_(cell_position);
    UnsignedInt package_index = phi < 0.0 ? 0 : 1;
    data_mesh_->assignDataPackageIndex(cell_index, package_index);

    UnsignedInt index_1d = data_mesh_->LinearCellIndex(cell_index);
    cell_contain_id_[index_1d] = 2;
    if (ABS(phi) > far_field_distance_)
    {
        cell_contain_id_[index_1d] = phi < 0.0 ? -1 : 1;
    }

    if (coarse_mesh_->isWithinCorePackage(
            cell_pkg_index_coarse_, pkg_cell_info_coarse_, cell_position))
    {
        Real signed_distance = shape_->findSignedDistance(cell_position);
        Vecd normal_direction = shape_->findNormalDirection(cell_position);
        Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
        if (measure < grid_spacing_)
        {
            UnsignedInt sort_index = base_dynamics_->SortIndexFromCellIndex(cell_index);
            data_mesh_->assignDataPackageIndex(cell_index, 2);
            data_mesh_->registerOccupied(sort_index, 1);
        }
    }
}
//=================================================================================================//
void InnerCellTagging::UpdateKernel::update(const Arrayi &cell_index)
{
    if (isInnerPackage(cell_index) && !data_mesh_->isInnerDataPackage(cell_index))
    {
        data_mesh_->registerOccupied(base_dynamics->SortIndexFromCellIndex(cell_index), 0);
    }
}
//=================================================================================================//
InitializeCellPackageInfo::InitializeCellPackageInfo(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      dv_pkg_cell_info_(data_mesh.dvPkgCellInfo()),
      bmv_cell_pkg_index_(data_mesh.getCellPackageIndex()) {}
//=================================================================================================//
void InitializeCellPackageInfo::UpdateKernel::update(const UnsignedInt &package_index)
{
    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs = data_mesh_->getOccupiedDataPackages();
    UnsignedInt sort_index = occupied_data_pkgs[package_index - num_singular_pkgs_].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    UnsignedInt linear_index = data_mesh_->LinearCellIndex(cell_index);
    cell_pkg_index_[linear_index] = package_index;
    std::pair<Arrayi, int> &metadata = pkg_cell_info_[package_index];
    metadata.first = cell_index;
    metadata.second = occupied_data_pkgs[package_index - num_singular_pkgs_].second;
}
//=================================================================================================//
InitializeCellNeighborhood::InitializeCellNeighborhood(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh), dv_pkg_cell_info_(data_mesh.dvPkgCellInfo()),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      bmv_cell_pkg_index_(data_mesh.getCellPackageIndex())
{
    CellNeighborhood *neighbor = dv_cell_neighborhood_.Data();
    for (UnsignedInt i = 0; i != num_singular_pkgs_; i++)
    {
        mesh_for_each(
            -Arrayi::Ones(), Arrayi::Ones() * 2,
            [&](const Arrayi &index)
            {
                neighbor[i](index + Arrayi::Ones()) = i;
            });
    }
}
//=============================================================================================//
void InitializeCellNeighborhood::UpdateKernel::update(const UnsignedInt &package_index)
{
    CellNeighborhood &current = cell_neighborhood_[package_index];
    Arrayi cell_index = pkg_cell_info_[package_index].first;
    mesh_for_each(
        -Arrayi::Ones(), Arrayi::Ones() * 2,
        [&](const Arrayi &index)
        {
            current(index + Arrayi::Ones()) =
                data_mesh_->PackageIndexFromCellIndex(cell_pkg_index_, cell_index + index);
        });
}
//=============================================================================================//
InitializeBasicPackageData::InitializeBasicPackageData(
    MeshWithGridDataPackagesType &data_mesh, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh), shape_(shape),
      far_field_distance(data_mesh.GridSpacing() * (Real)data_mesh.BufferWidth()),
      dv_pkg_cell_info_(data_mesh.dvPkgCellInfo()),
      mv_phi_(*data_mesh.registerMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.registerMeshVariable<int>("NearInterfaceID"))
{
    data_mesh.addMeshVariableToWrite<Real>("LevelSet");
    initializeSingularPackages(0, -far_field_distance);
    initializeSingularPackages(1, far_field_distance);
}
//=============================================================================================//
NearInterfaceCellTagging::NearInterfaceCellTagging(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      bmv_cell_contain_id_(*data_mesh.getBKGMeshVariable<int>("CellContainID")),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")) {}
//=============================================================================================//
CellContainDiffusion::CellContainDiffusion(
    MeshWithGridDataPackagesType &data_mesh, SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh),
      bmv_cell_contain_id_(*data_mesh.getBKGMeshVariable<int>("CellContainID")),
      bmv_cell_package_index_(data_mesh.getCellPackageIndex()),
      sv_count_modified_(sv_count_modified) {}
//=============================================================================================//
UpdateLevelSetGradient::UpdateLevelSetGradient(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
UpdateKernelIntegrals::UpdateKernelIntegrals(
    MeshWithGridDataPackagesType &data_mesh, NeighborMethod<SingleValued> &neighbor_method)
    : BaseMeshLocalDynamics(data_mesh), neighbor_method_(neighbor_method),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      dv_pkg_cell_info_(data_mesh.dvPkgCellInfo()),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      bmv_cell_pkg_index_(data_mesh.getCellPackageIndex()),
      mv_kernel_weight_(*data_mesh.registerMeshVariable<Real>("KernelWeight")),
      mv_kernel_gradient_(*data_mesh.registerMeshVariable<Vecd>("KernelGradient")),
      mv_kernel_second_gradient_(*data_mesh.registerMeshVariable<Matd>("KernelSecondGradient")),
      far_field_distance(data_mesh.GridSpacing() * (Real)data_mesh.BufferWidth())
{
    initializeSingularPackages(0, -far_field_distance);
    initializeSingularPackages(1, far_field_distance);
}
//=============================================================================================//
MarkCutInterfaces::MarkCutInterfaces(MeshWithGridDataPackagesType &data_mesh, Real perturbation_ratio)
    : BaseMeshLocalDynamics(data_mesh),
      threshold_(data_spacing_ * sqrt(Dimensions)), perturbation_ratio_(perturbation_ratio),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
MarkNearInterface::MarkNearInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      threshold_(data_spacing_ * sqrt(Dimensions)),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
ReinitializeLevelSet::ReinitializeLevelSet(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
RedistanceInterface::RedistanceInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
DiffuseLevelSetSign::DiffuseLevelSetSign(
    MeshWithGridDataPackagesType &data_mesh, SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      sv_count_modified_(sv_count_modified) {}
//=================================================================================================//
} // namespace SPH
