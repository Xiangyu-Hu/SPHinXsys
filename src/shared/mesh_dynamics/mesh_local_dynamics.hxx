#include "mesh_local_dynamics.h"

#ifndef MESH_LOCAL_DYNAMICS_HXX
#define MESH_LOCAL_DYNAMICS_HXX

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitialCellTagging::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
      index_handler_(encloser.index_handler_),
      grid_spacing_(index_handler_.GridSpacing()),
      shape_(&encloser.shape_),
      cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitialCellTaggingFromCoarse::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : shape_(&encloser.shape_),
      occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
      index_handler_(encloser.index_handler_),
      coarse_index_handler_(encloser.coarse_mesh_.getIndexHandler()),
      grid_spacing_(index_handler_.GridSpacing()),
      far_field_distance_(grid_spacing_ * (Real)index_handler_.BufferWidth()),
      probe_coarse_phi_(ex_policy, &encloser.coarse_mesh_),
      cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
      cell_pkg_index_coarse_(encloser.bmv_cell_pkg_index_coarse_.DelegatedData(ex_policy)),
      pkg_type_coarse_(encloser.dv_pkg_type_coarse_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InnerCellTagging::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : occupied_data_pkgs_(&encloser.occupied_data_pkgs_),
      index_handler_(encloser.index_handler_),
      cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
InitializeCellNeighborhood::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
      num_singular_pkgs_(encloser.data_mesh_.NumSingularPackages()) {}
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
      cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
CellContainDiffusion::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
      cell_package_index_(encloser.bmv_cell_package_index_.DelegatedData(ex_policy)),
      count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateLevelSetGradient::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateKernelIntegrals::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      kernel_weight_(encloser.mv_kernel_weight_.DelegatedData(ex_policy)),
      kernel_gradient_(encloser.mv_kernel_gradient_.DelegatedData(ex_policy)),
      kernel_second_gradient_(encloser.mv_kernel_second_gradient_.DelegatedData(ex_policy)),
      pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_.DelegatedData(ex_policy)),
      kernel_(ex_policy, encloser.neighbor_method_),
      index_handler_(encloser.index_handler_),
      data_spacing_(index_handler_.DataSpacing()),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
      probe_signed_distance_(ex_policy, &encloser.data_mesh_),
      cutoff_radius_(encloser.neighbor_method_.CutOffRadius()),
      depth_(static_cast<int>(std::ceil((cutoff_radius_ - Eps) / data_spacing_))) {}
//=================================================================================================//
inline void UpdateKernelIntegrals::UpdateKernel::update(const UnsignedInt &package_index)
{
    assignByDataIndex(
        kernel_weight_[package_index], [&](const Arrayi &data_index) -> Real
        { return computeKernelIntegral(package_index, data_index); });
    assignByDataIndex(
        kernel_gradient_[package_index], [&](const Arrayi &data_index) -> Vecd
        { return computeKernelGradientIntegral(package_index, data_index); });
    assignByDataIndex(
        kernel_second_gradient_[package_index], [&](const Arrayi &data_index) -> Matd
        { return computeKernelSecondGradientIntegral(package_index, data_index); });
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
ReinitializeLevelSet::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Real ReinitializeLevelSet::UpdateKernel::upwindDifference(Real sign, Real df_p, Real df_n)
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
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
MarkCutInterfaces::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      threshold_(index_handler_.DataSpacing() * sqrt(Dimensions)),
      perturbation_(threshold_ * encloser.perturbation_ratio_),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
MarkNearInterface::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      threshold_(index_handler_.GridSpacing() * sqrt(Dimensions)),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
RedistanceInterface::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
DiffuseLevelSetSign::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)) {}
//=================================================================================================//
//=================================================================================================//
} // namespace SPH

#endif // MESH_LOCAL_DYNAMICS_HXX