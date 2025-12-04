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
inline void NearInterfaceCellTagging::UpdateKernel::update(const UnsignedInt &package_index)
{
    UnsignedInt index_1d = pkg_1d_cell_index_[package_index];
    MeshVariableData<Real> &pkg_phi = phi_[package_index];
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
      cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
      cell_package_index_(encloser.bmv_cell_package_index_.DelegatedData(ex_policy)),
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
            cell_package_index_[index_1d] = 0; // inside far field package updated
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
            cell_package_index_[index_1d] = 1; // outside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateLevelSetGradient::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &index)
                  {
                      phi_gradient_[package_index](index) =
                          regularizedCentralDifference(
                              phi_, cell_neighborhood_[package_index], index,
                              [](Real dp, Real dm)
                              {
                                  return 0.5 * (dp + dm); // no upwinding
                              }) /
                          data_spacing_;
                  });
}
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
      depth_(static_cast<int>(std::ceil((cutoff_radius_ - Eps) / data_spacing_))),
      bounding_box_(BoundingBoxi(Arrayi::Constant(depth_))) {}
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
template <typename DataType, typename FunctionByDataIndex>
DataType UpdateKernelIntegrals::UpdateKernel::
    computeIntegral(Real phi, const UnsignedInt &package_index, const Arrayi &grid_index,
                    const DataType &initial_value, const FunctionByDataIndex &function_by_index)
{
    DataType integral = initial_value;
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each(
            bounding_box_.lower_, bounding_box_.upper_ + Arrayi::Ones(),
            [&](const Arrayi &search_index)
            {
                DataPackagePair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + search_index);
                Real phi_neighbor = phi_[neighbor_meta.first](neighbor_meta.second);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first](neighbor_meta.second);
                    Vecd displacement = -search_index.cast<Real>().matrix() * data_spacing_;
                    if (displacement.norm() < cutoff_radius_)
                        integral += function_by_index(displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return integral;
}
//=============================================================================================//
inline Real UpdateKernelIntegrals::UpdateKernel::
    computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    Real integral = computeIntegral(phi, package_index, data_index, 0.0,
                                    [&](const Vecd &displacement) -> Real
                                    { return kernel_.W(displacement); });
    return phi > cutoff_radius_ ? 1.0 : integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Vecd UpdateKernelIntegrals::UpdateKernel::
    computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    Vecd integral = Vecd::Zero();
    integral = computeIntegral(phi, package_index, data_index, integral,
                               [&](const Vecd &displacement) -> Vecd
                               { return kernel_.dW(displacement) * displacement /
                                        (displacement.norm() + TinyReal); });
    return integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Matd UpdateKernelIntegrals::UpdateKernel::
    computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    Matd integral = Matd::Zero();
    integral = computeIntegral(phi, package_index, data_index, integral,
                               [&](const Vecd &displacement) -> Matd
                               { return kernel_.d2W(displacement) *
                                        displacement * displacement.transpose() /
                                        (displacement.squaredNorm() + TinyReal); });
    return integral * data_spacing_ * data_spacing_;
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
//=============================================================================================//
inline void ReinitializeLevelSet::UpdateKernel::update(const UnsignedInt &package_index)
{
    auto &phi_pkg = phi_[package_index];
    auto &near_interface_id_pkg = near_interface_id_[package_index];
    auto &neighborhood = cell_neighborhood_[package_index];

    mesh_for_each(
        Arrayi::Zero(), Arrayi::Constant(pkg_size),
        [&](const Arrayi &index)
        {
            if (near_interface_id_pkg(index) != 0)
            {
                Real phi_0 = phi_pkg(index);
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);

                Vecd difference = regularizedCentralDifference(
                    phi_, neighborhood, index,
                    [&](Real dp, Real dm)
                    {
                        return upwindDifference(sign, dp, dm);
                    });

                phi_pkg(index) -= sign * (difference.norm() - data_spacing_) / Real(Dimensions);
            }
        });
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