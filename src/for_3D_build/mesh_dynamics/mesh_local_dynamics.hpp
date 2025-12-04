#ifndef MESH_LOCAL_DYNAMICS_3D_HPP
#define MESH_LOCAL_DYNAMICS_3D_HPP

#include "data_package_function.hpp"
#include "mesh_local_dynamics.hxx"

namespace SPH
{
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