#ifndef MESH_LOCAL_DYNAMICS_2D_HPP
#define MESH_LOCAL_DYNAMICS_2D_HPP

#include "data_package_function.hpp"
#include "mesh_local_dynamics.hxx"

namespace SPH
{
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