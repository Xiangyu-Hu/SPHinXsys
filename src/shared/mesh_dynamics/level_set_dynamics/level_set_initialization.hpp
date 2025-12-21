#include "level_set_initialization.h"

#ifndef LEVEL_SET_INITIALIZATION_HPP
#define LEVEL_SET_INITIALIZATION_HPP

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitialCellTagging::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      cell_pkg_index_(encloser.mcv_cell_pkg_index_.DelegatedData(ex_policy)),
      index_handler_(encloser.index_handler_),
      grid_spacing_(index_handler_.GridSpacing()),
      shape_(&encloser.shape_),
      cell_contain_id_(encloser.mcv_cell_contain_id_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitialCellTaggingFromCoarse::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : shape_(&encloser.shape_),
      occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      boundary_pkg_index_offset_(encloser.boundary_pkg_index_offset_),
      cell_pkg_index_(encloser.mcv_cell_pkg_index_.DelegatedData(ex_policy)),
      index_handler_(encloser.index_handler_),
      coarse_index_handler_(
          encloser.coarse_mesh_.getMesh(encloser.coarse_resolution_level_)),
      grid_spacing_(index_handler_.GridSpacing()),
      far_field_distance_(grid_spacing_ * (Real)index_handler_.BufferWidth()),
      probe_coarse_phi_(ex_policy, encloser.coarse_mesh_, "LevelSet"),
      cell_contain_id_(encloser.mcv_cell_contain_id_.DelegatedData(ex_policy)),
      cell_pkg_index_coarse_(encloser.mcv_cell_pkg_index_coarse_.DelegatedData(ex_policy)),
      pkg_type_coarse_(encloser.dv_pkg_type_coarse_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InnerCellTagging::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      index_handler_(encloser.index_handler_),
      cell_pkg_index_(encloser.mcv_cell_pkg_index_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitializeCellNeighborhood::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      cell_pkg_index_(encloser.mcv_cell_pkg_index_.DelegatedData(ex_policy)),
      num_boundary_pkgs_(encloser.data_mesh_.NumBoundaryPackages()) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitializeBasicPackageData::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_.DelegatedData(ex_policy)),
      shape_(&encloser.shape_),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
NearInterfaceCellTagging::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_.DelegatedData(ex_policy)),
      cell_contain_id_(encloser.mcv_cell_contain_id_.DelegatedData(ex_policy)),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void NearInterfaceCellTagging::UpdateKernel::update(const UnsignedInt &package_index)
{
    UnsignedInt index_1d = pkg_1d_cell_index_[package_index];
    PackageVariableData<Real> &pkg_phi = phi_[package_index];
    Real phi0 = pkg_phi(Arrayi::Zero());
    cell_contain_id_[index_1d] = phi0 > 0.0 ? 1 : -1;
    bool is_sign_changed = mesh_any_of(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                                       [&](const Arrayi &data_index) -> bool
                                       {
                                           return pkg_phi(data_index) * phi0 < 0.0;
                                       });
    if (is_sign_changed)
        cell_contain_id_[index_1d] = 0;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
CellContainDiffusion::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      boundary_pkg_index_offset_(encloser.boundary_pkg_index_offset_),
      cell_contain_id_(encloser.mcv_cell_contain_id_.DelegatedData(ex_policy)),
      cell_pkg_index_(encloser.mcv_cell_package_index_.DelegatedData(ex_policy)),
      count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void CellContainDiffusion::UpdateKernel::update(const Arrayi &cell_index)
{
    UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
    if (cell_contain_id_[index_1d] == 2)
    {
        if (mesh_any_of(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                index_handler_.AllCells().min(cell_index + 2 * Arrayi::Ones()),
                [&](const Arrayi &index)
                {
                    UnsignedInt neighbor_1d = index_handler_.LinearCellIndex(index);
                    return cell_contain_id_[neighbor_1d] == -1;
                }))
        {
            cell_contain_id_[index_1d] = -1;
            cell_pkg_index_[index_1d] = boundary_pkg_index_offset_; // inside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
        else if (mesh_any_of(
                     Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                     index_handler_.AllCells().min(cell_index + 2 * Arrayi::Ones()),
                     [&](const Arrayi &index)
                     {
                         UnsignedInt neighbor_1d = index_handler_.LinearCellIndex(index);
                         return cell_contain_id_[neighbor_1d] == 1;
                     }))
        {
            cell_contain_id_[index_1d] = 1;
            cell_pkg_index_[index_1d] = boundary_pkg_index_offset_ + 1; // outside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
    }
}
//=================================================================================================//
} // namespace SPH
#endif // LEVEL_SET_INITIALIZATION_HPP