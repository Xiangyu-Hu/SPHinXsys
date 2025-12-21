#include "level_set_correction.h"

#ifndef LEVEL_SET_CORRECTION_HPP
#define LEVEL_SET_CORRECTION_HPP

namespace SPH
{
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
                        PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                        PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                            PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                            PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                            PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                                 PackageDataPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
//=============================================================================================//
template <class ExecutionPolicy>
CleanInterface<ExecutionPolicy>::CleanInterface(
    SparseMeshField<4> &mesh_data, UnsignedInt resolution_level,
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method, Real refinement_ratio)
    : RepeatTimes(), BaseDynamics<void>(),
      neighbor_method_(neighbor_method),
      update_level_set_gradient{mesh_data, resolution_level},
      update_kernel_integrals{mesh_data, resolution_level, neighbor_method_},
      mark_cut_interfaces{mesh_data, resolution_level, 0.5 * refinement_ratio},
      redistance_interface{mesh_data, resolution_level},
      reinitialize_level_set{mesh_data, resolution_level} {}
//=============================================================================================//
template <class ExecutionPolicy>
CorrectTopology<ExecutionPolicy>::CorrectTopology(
    SparseMeshField<4> &mesh_data, UnsignedInt resolution_level,
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method)
    : BaseDynamics<void>(), neighbor_method_(neighbor_method),
      update_level_set_gradient(mesh_data, resolution_level),
      update_kernel_integrals(mesh_data, resolution_level, neighbor_method),
      mark_near_interface(mesh_data, resolution_level),
      diffuse_level_set_sign(mesh_data, resolution_level, sv_count_modified_) {}
//=================================================================================================//
} // namespace SPH
#endif // LEVEL_SET_CORRECTION_HPP