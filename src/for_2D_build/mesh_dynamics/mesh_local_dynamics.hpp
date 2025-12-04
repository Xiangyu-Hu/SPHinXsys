#ifndef MESH_LOCAL_DYNAMICS_2D_HPP
#define MESH_LOCAL_DYNAMICS_2D_HPP

#include "data_package_function.hpp"
#include "mesh_local_dynamics.hxx"

namespace SPH
{
//=============================================================================================//
inline void MarkNearInterface::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            near_interface_id_[package_index][i][j] = 3; // undetermined
            Real phi0 = phi_[package_index][i][j];
            if (ABS(phi0) < 2.0 * threshold_) // only consider data close to the interface
            {
                bool is_sign_changed = mesh_any_of2d<-1, 2>( // check in the 3x3x3 neighborhood
                    [&](int l, int m) -> bool
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + l, j + m), cell_neighborhood_[package_index]);

                        return phi0 * phi_[neighbour_index.first]
                                          [neighbour_index.second[0]]
                                          [neighbour_index.second[1]] <
                               0.0;
                    });

                if (is_sign_changed)
                {
                    if (ABS(phi0) < threshold_)
                    {
                        near_interface_id_[package_index][i][j] = 0; // cut cell
                    }
                }
                else
                {
                    near_interface_id_[package_index][i][j] = phi0 > 0.0 ? 1 : -1; // in the band
                }
            }
        });
}
//=============================================================================================//
inline void RedistanceInterface::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            int near_interface_id = near_interface_id_[package_index][i][j];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each2d<-1, 2>(
                    [&](int r, int s)
                    {
                        DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + r, j + s), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id =
                            near_interface_id_[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            DataPackagePair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j] = 1;
                }
            }
        });
}
//=============================================================================================//
inline void DiffuseLevelSetSign::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_[package_index][i][j]) > 1)
            {
                mesh_find_if2d<-1, 2>(
                    [&](int l, int m) -> bool
                    {
                        DataPackagePair neighbour_index =
                            NeighbourIndexShift<pkg_size>(
                                Arrayi(i + l, j + m), cell_neighborhood_[package_index]);
                        int near_interface_id =
                            near_interface_id_[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_[package_index][i][j];
                            near_interface_id_[package_index][i][j] = near_interface_id;
                            phi_[package_index][i][j] =
                                near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif // MESH_LOCAL_DYNAMICS_2D_HPP