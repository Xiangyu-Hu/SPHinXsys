#include "mesh_local_dynamics.h"

#ifndef MESH_LOCAL_DYNAMICS_HPP
#define MESH_LOCAL_DYNAMICS_HPP

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
      kernel_(ex_policy, encloser.neighbor_method_),
      data_spacing_(encloser.index_handler_.DataSpacing()),
      data_cell_volume_(math::pow(data_spacing_, Dimensions)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      cutoff_radius_(encloser.neighbor_method_.CutOffRadius()),
      bounding_box_(BoundingBoxi(Arrayi::Constant(
          static_cast<int>(std::ceil((cutoff_radius_ - Eps) / data_spacing_))))) {}
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
inline Real UpdateKernelIntegrals::UpdateKernel::
    CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing)
{
    Real squared_norm_inv = 1.0 / (phi_gradient.squaredNorm() + TinyReal);
    Real volume_fraction(0);
    for (UnsignedInt i = 0; i != Dimensions; ++i)
    {
        volume_fraction += phi_gradient[i] * phi_gradient[i] * squared_norm_inv *
                           Heaviside(phi / (ABS(phi_gradient[i]) + TinyReal), 0.5 * data_spacing);
    }
    return volume_fraction;
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
    return phi > cutoff_radius_ ? 1.0 : integral * data_cell_volume_;
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
    return integral * data_cell_volume_;
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
    return integral * data_cell_volume_;
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
//=============================================================================================//
inline void MarkCutInterfaces::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageData<Real, 5> corner_averages;
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Constant(5),
        [&](const Arrayi &index)
        {
            corner_averages(index) =
                CornerAverage(phi_, index, Arrayi::Constant(-1),
                              cell_neighborhood_[package_index], (Real)0);
        });

    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &index)
                  {
                      // first assume far cells
                      Real phi_0 = phi_addrs(index);
                      int near_interface_id = phi_0 > 0.0 ? 2 : -2;
                      if (fabs(phi_0) < perturbation_)
                      {
                          near_interface_id = 0;
                          Real phi_average_0 = corner_averages(index);
                          // find inner and outer cut cells
                          mesh_for_each(
                              Arrayi::Zero(), Arrayi::Constant(2),
                              [&](const Arrayi &shift)
                              {
                                  Real phi_average = corner_averages(index + shift);
                                  if ((phi_average_0 - perturbation_) * (phi_average - perturbation_) < 0.0)
                                      near_interface_id = 1;
                                  if ((phi_average_0 + perturbation_) * (phi_average + perturbation_) < 0.0)
                                      near_interface_id = -1;
                              });
                          // find zero cut cells
                          mesh_for_each(
                              Arrayi::Zero(), Arrayi::Constant(2),
                              [&](const Arrayi &shift)
                              {
                                  Real phi_average = corner_averages(index + shift);
                                  if (phi_average_0 * phi_average < 0.0)
                                      near_interface_id = 0;
                              });
                      }
                      // assign this is to package
                      near_interface_id_addrs(index) = near_interface_id;
                  });
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
MarkNearInterface::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : index_handler_(encloser.index_handler_),
      threshold_(index_handler_.GridSpacing() * sqrt(Dimensions)),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=============================================================================================//
inline void MarkNearInterface::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Constant(pkg_size),
        [&](const Arrayi &index)
        {
            near_interface_id_[package_index](index) = 3; // undetermined
            Real phi0 = phi_[package_index](index);
            if (ABS(phi0) < 2.0 * threshold_) // only consider data close to the interface
            {
                bool is_sign_changed = mesh_any_of(
                    Arrayi::Constant(-1), Arrayi::Constant(2), // check in the 3x3x3 neighborhood
                    [&](const Arrayi &shift) -> bool
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            index + shift, cell_neighborhood_[package_index]);

                        return phi0 * phi_[neighbour_index.first](neighbour_index.second) < 0.0;
                    });

                if (is_sign_changed)
                {
                    if (ABS(phi0) < threshold_)
                    {
                        near_interface_id_[package_index](index) = 0; // cut cell
                    }
                }
                else
                {
                    near_interface_id_[package_index](index) = phi0 > 0.0 ? 1 : -1; // in the band
                }
            }
        });
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
RedistanceInterface::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=============================================================================================//
inline void RedistanceInterface::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Constant(pkg_size),
        [&](const Arrayi &index)
        {
            int near_interface_id = near_interface_id_[package_index](index);
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each(
                    Arrayi::Constant(-1), Arrayi::Constant(2), // check in the 3x3x3 neighborhood
                    [&](const Arrayi &shift)
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            index + shift, cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id =
                            near_interface_id_[neighbour_index.first](neighbour_index.second);
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each(
                        Arrayi::Constant(-4), Arrayi::Constant(5), // check in the 4x4x4 neighborhood
                        [&](const Arrayi &shift)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                index + shift, cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id(neighbour_index.second) >= 1)
                            {
                                Real phi_p_ = neighbor_phi(neighbour_index.second);
                                Vecd norm_to_face = neighbor_phi_gradient(neighbour_index.second);
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(
                                    min_distance_p,
                                    (Vecd(shift.cast<Real>()) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index](index) = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index](index) = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each(
                        Arrayi::Constant(-4), Arrayi::Constant(5), // check in the 4x4x4 neighborhood
                        [&](const Arrayi &shift)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                index + shift, cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id(neighbour_index.second) <= -1)
                            {
                                Real phi_n_ = neighbor_phi(neighbour_index.second);
                                Vecd norm_to_face = neighbor_phi_gradient(neighbour_index.second);
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(
                                    min_distance_n,
                                    (Vecd(shift.cast<Real>()) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index](index) = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index](index) = 1;
                }
            }
        });
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
DiffuseLevelSetSign::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void DiffuseLevelSetSign::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Constant(pkg_size),
        [&](const Arrayi &index)
        {
            if (near_interface_id_[package_index](index) == 3) // check undetermined only
            {
                Real phi = phi_[package_index](index);
                if (mesh_any_of(
                        Arrayi::Constant(-1), Arrayi::Constant(2), // check in the 3x3x3 neighborhood
                        [&](const Arrayi &shift) -> bool
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                index + shift, cell_neighborhood_[package_index]);
                            return near_interface_id_[neighbour_index.first](neighbour_index.second) == 1;
                        }))
                {
                    phi_[package_index](index) = ABS(phi);
                    near_interface_id_[package_index](index) = 1; // mark as positive band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
                else if (mesh_any_of(
                             Arrayi::Constant(-1), Arrayi::Constant(2), // check in the 3x3x3 neighborhood
                             [&](const Arrayi &shift) -> bool
                             {
                                 DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                     index + shift, cell_neighborhood_[package_index]);
                                 return near_interface_id_[neighbour_index.first](neighbour_index.second) == -1;
                             }))
                {
                    phi_[package_index](index) = -ABS(phi);
                    near_interface_id_[package_index](index) = -1; // mark as negative band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
            }
        });
}
//=================================================================================================//
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_HPP