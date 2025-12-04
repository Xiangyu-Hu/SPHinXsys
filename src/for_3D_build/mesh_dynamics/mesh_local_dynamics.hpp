#ifndef MESH_LOCAL_DYNAMICS_3D_HPP
#define MESH_LOCAL_DYNAMICS_3D_HPP

#include "data_package_function.hpp"
#include "mesh_local_dynamics.hxx"

namespace SPH
{
//=============================================================================================//
inline void MarkCutInterfaces::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageData<Real, 5> corner_averages;
    mesh_for_each3d<0, 5>(
        [&](int i, int j, int k)
        {
            corner_averages[i][j][k] = CornerAverage(phi_, Arrayi(i, j, k), Arrayi(-1, -1, -1),
                                                     cell_neighborhood_[package_index], (Real)0);
        });

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j][k];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < perturbation_)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j][k];
                // find inner and outer cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if ((phi_average_0 - perturbation_) * (phi_average - perturbation_) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + perturbation_) * (phi_average + perturbation_) < 0.0)
                            near_interface_id = -1;
                    });
                // find zero cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if (phi_average_0 * phi_average < 0.0)
                            near_interface_id = 0;
                    });
            }
            // assign this is to package
            near_interface_id_addrs[i][j][k] = near_interface_id;
        });
}
//=============================================================================================//
inline void MarkNearInterface::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            near_interface_id_[package_index][i][j][k] = 3; // undetermined
            Real phi0 = phi_[package_index][i][j][k];
            if (ABS(phi0) < 2.0 * threshold_) // only consider data close to the interface
            {
                bool is_sign_changed = mesh_any_of3d<-1, 2>( // check in the 3x3x3 neighborhood
                    [&](int l, int m, int n) -> bool
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);

                        return phi0 * phi_[neighbour_index.first]
                                          [neighbour_index.second[0]]
                                          [neighbour_index.second[1]]
                                          [neighbour_index.second[2]] <
                               0.0;
                    });

                if (is_sign_changed)
                {
                    if (ABS(phi0) < threshold_)
                    {
                        near_interface_id_[package_index][i][j][k] = 0; // cut cell
                    }
                }
                else
                {
                    near_interface_id_[package_index][i][j][k] = phi0 > 0.0 ? 1 : -1; // in the band
                }
            }
        });
}
//=============================================================================================//
inline void RedistanceInterface::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            int near_interface_id = near_interface_id_[package_index][i][j][k];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each3d<-1, 2>(
                    [&](int r, int s, int t)
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + r, j + s, k + t), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id = near_interface_id_[neighbour_index.first]
                                                                           [neighbour_index.second[0]]
                                                                           [neighbour_index.second[1]]
                                                                           [neighbour_index.second[2]];
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]]
                                                                         [neighbour_index.second[1]]
                                                                         [neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(
                                    min_distance_p,
                                    (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j][k] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j][k] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]]
                                                                         [neighbour_index.second[1]]
                                                                         [neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(
                                    min_distance_n,
                                    (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j][k] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j][k] = 1;
                }
            }
        });
}
//=============================================================================================//
inline void DiffuseLevelSetSign::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            if (near_interface_id_[package_index][i][j][k] == 3) // check undetermined only
            {
                Real phi = phi_[package_index][i][j][k];
                if (mesh_any_of3d<-1, 2>(
                        [&](int l, int m, int n) -> bool
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                            return near_interface_id_[neighbour_index.first]
                                                     [neighbour_index.second[0]]
                                                     [neighbour_index.second[1]]
                                                     [neighbour_index.second[2]] == 1;
                        }))
                {
                    phi_[package_index][i][j][k] = ABS(phi);
                    near_interface_id_[package_index][i][j][k] = 1; // mark as positive band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
                else if (mesh_any_of3d<-1, 2>(
                             [&](int l, int m, int n) -> bool
                             {
                                 DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                     Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                                 return near_interface_id_[neighbour_index.first]
                                                          [neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] == -1;
                             }))
                {
                    phi_[package_index][i][j][k] = -ABS(phi);
                    near_interface_id_[package_index][i][j][k] = -1; // mark as negative band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif // MESH_LOCAL_DYNAMICS_3D_HPP