#include "mesh_local_dynamics.hpp"

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

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
void InitialCellTagging::UpdateKernel::update(const Arrayi &cell_index)
{
    Vecd cell_position = data_mesh_->CellPositionFromIndex(cell_index);
    Real signed_distance = shape_->findSignedDistance(cell_position);
    Vecd normal_direction = shape_->findNormalDirection(cell_position);
    Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
    if (measure < grid_spacing_)
    {
        size_t sort_index = base_dynamics->SortIndexFromCellIndex(cell_index);
        data_mesh_->assignDataPackageIndex(cell_index, 2);
        data_mesh_->registerOccupied(sort_index, 1);
    }
    else
    { // this is not reliable for non-water tightening shape, to be corrected later
        size_t package_index = shape_->checkContain(cell_position) ? 0 : 1;
        data_mesh_->assignDataPackageIndex(cell_index, package_index);
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
void InitializeIndexMesh::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = data_mesh_->occupied_data_pkgs_[package_index - num_singular_pkgs_].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    data_mesh_->assignDataPackageIndex(cell_index, package_index);
}
//=================================================================================================//
void InitialCellTaggingFromCoarse::UpdateKernel::update(const Arrayi &cell_index)
{
    Vecd cell_position = data_mesh_->CellPositionFromIndex(cell_index);
    size_t package_index = probe_coarse_phi_(cell_position) < 0.0 ? 0 : 1;
    data_mesh_->assignDataPackageIndex(cell_index, package_index);
    if (coarse_mesh_->isWithinCorePackage(coarse_mesh_->cell_pkg_index_.Data(),
                                          coarse_mesh_->pkg_cell_info_.Data(),
                                          cell_position))
    {
        Real signed_distance = shape_->findSignedDistance(cell_position);
        Vecd normal_direction = shape_->findNormalDirection(cell_position);
        Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
        if (measure < grid_spacing_)
        {
            size_t sort_index = base_dynamics_->SortIndexFromCellIndex(cell_index);
            data_mesh_->assignDataPackageIndex(cell_index, 2);
            data_mesh_->registerOccupied(sort_index, 1);
        }
    }
}
//=================================================================================================//
InitializeCellNeighborhood::InitializeCellNeighborhood(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh), dv_pkg_cell_info_(data_mesh.getMetaDataCell()),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      bmv_cell_pkg_index_(data_mesh.getCellPackageIndex()) {}
//=============================================================================================//
void InitializeCellNeighborhood::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = data_mesh_->occupied_data_pkgs_[package_index - num_singular_pkgs_].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    CellNeighborhood &current = cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = pkg_cell_info_[package_index];
    metadata.first = cell_index;
    metadata.second = data_mesh_->occupied_data_pkgs_[package_index - num_singular_pkgs_].second;

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
      dv_pkg_cell_info_(data_mesh.getMetaDataCell()),
      mv_phi_(*data_mesh.registerMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.registerMeshVariable<int>("CellNearInterfaceID"))
{
    data_mesh.addMeshVariableToWrite<Real>("LevelSet");
    initializeSingularPackages(0, -far_field_distance);
    initializeSingularPackages(1, far_field_distance);
}
//=============================================================================================//
NearInterfaceCellTagging::NearInterfaceCellTagging(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      dv_cell_near_interface_id_(
          data_mesh.registerBKGMeshVariable<int>(
              "CellNearInterfaceID", [&](UnsignedInt index)
              { return MaxInt; })),
      dv_phi_(data_mesh.getMeshVariable<Real>("LevelSet")) {}
//=============================================================================================//
SingularPackageCorrection::SingularPackageCorrection(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      dv_cell_near_interface_id_(data_mesh.getBKGMeshVariable<int>("CellNearInterfaceID")),
      sv_count_modified_("IsModified", 0) {}
//=============================================================================================//
UpdateLevelSetGradient::UpdateLevelSetGradient(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
UpdateKernelIntegrals::UpdateKernelIntegrals(
    MeshWithGridDataPackagesType &data_mesh, KernelTabulatedCK *kernel, Real global_h_ratio)
    : BaseMeshLocalDynamics(data_mesh), kernel_(kernel), global_h_ratio_(global_h_ratio),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      dv_pkg_cell_info_(data_mesh.getMetaDataCell()),
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
MarkNearInterface::MarkNearInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("CellNearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
ReinitializeLevelSet::ReinitializeLevelSet(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("CellNearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
RedistanceInterface::RedistanceInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("CellNearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
DiffuseLevelSetSign::DiffuseLevelSetSign(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("CellNearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=================================================================================================//
} // namespace SPH
