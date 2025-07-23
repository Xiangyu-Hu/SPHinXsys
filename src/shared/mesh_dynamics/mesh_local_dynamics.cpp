#include "mesh_local_dynamics.hpp"

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
void InitializeDataInACell::UpdateKernel::update(const Arrayi &cell_index)
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
    {
        size_t package_index = shape_->checkContain(cell_position) ? 0 : 1;
        data_mesh_->assignDataPackageIndex(cell_index, package_index);
    }
}
//=================================================================================================//
void TagACellIsInnerPackage::UpdateKernel::update(const Arrayi &cell_index)
{
    if (isInnerPackage(cell_index) && !data_mesh_->isInnerDataPackage(cell_index))
    {
        data_mesh_->registerOccupied(base_dynamics->SortIndexFromCellIndex(cell_index), 0);
    }
}
//=================================================================================================//
void InitializeIndexMesh::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = data_mesh_->occupied_data_pkgs_[package_index - 2].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    data_mesh_->assignDataPackageIndex(cell_index, package_index);
}
//=================================================================================================//
void InitializeDataInACellFromCoarse::UpdateKernel::update(const Arrayi &cell_index)
{
    Vecd cell_position = data_mesh_->CellPositionFromIndex(cell_index);
    size_t package_index = probe_coarse_phi_(cell_position) < 0.0 ? 0 : 1;
    data_mesh_->assignDataPackageIndex(cell_index, package_index);
    if(coarse_mesh_->isWithinCorePackage(coarse_mesh_->cell_package_index_.Data(),
                                         coarse_mesh_->meta_data_cell_.Data(),
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
} // namespace SPH
//=================================================================================================//
