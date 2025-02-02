#include "mesh_local_dynamics.hpp"
#include "mesh_iterators.hpp"

namespace SPH
{
//=============================================================================================//
size_t BaseMeshLocalDynamics::SortIndexFromCellIndex(const Arrayi &cell_index)
{
    return cell_index[0] * all_cells_[1] * all_cells_[2] + cell_index[1] * all_cells_[2] + cell_index[2];
}
//=============================================================================================//
Arrayi BaseMeshLocalDynamics::CellIndexFromSortIndex(const size_t &sort_index)
{
    Array3i cell_index;
    cell_index[0] = sort_index / (all_cells_[1] * all_cells_[2]);
    cell_index[1] = (sort_index / all_cells_[2]) % all_cells_[1];
    cell_index[2] = sort_index % all_cells_[2];

    return cell_index;
}
//=============================================================================================//
std::pair<size_t, Arrayi> BaseMeshLocalDynamics::
    NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour)
{
    std::pair<size_t, Arrayi> result;
    Arrayi neighbour_index = (shift_index + pkg_size * Arrayi::Ones()) / pkg_size;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]][neighbour_index[2]];
    result.second = (shift_index + pkg_size * Arrayi::Ones()) - neighbour_index * pkg_size;

    return result;
}
//=============================================================================================//
void InitializeDataForSingularPackage::update(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.DataField()[package_index];
    auto &near_interface_id = near_interface_id_.DataField()[package_index];
    auto &phi_gradient = phi_gradient_.DataField()[package_index];
    auto &kernel_weight = kernel_weight_.DataField()[package_index];
    auto &kernel_gradient = kernel_gradient_.DataField()[package_index];

    for_each_cell_data(
        [&](int i, int j, int k)
        {
            phi[i][j][k] = far_field_level_set;
            near_interface_id[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j][k] = Vecd::Ones();
            kernel_weight[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j][k] = Vec3d::Zero();
        });
}
//=============================================================================================//
bool TagACellIsInnerPackage::UpdateKernel::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array3i::Zero().max(cell_index - Array3i::Ones()),
        all_cells_.min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return mesh_data_->isInnerDataPackage(Arrayi(l, m, n));    //actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_->occupied_data_pkgs_[package_index-2].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    CellNeighborhood &current = cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = meta_data_cell_[package_index];
    metadata.first = cell_index;
    metadata.second = mesh_data_->occupied_data_pkgs_[package_index-2].second;
    for (int l = -1; l < 2; l++)
        for (int m = -1; m < 2; m++)
            for (int n = -1; n < 2; n++)
            {
                current[l + 1][m + 1][n + 1] = mesh_data_->PackageIndexFromCellIndex(cell_package_index_, cell_index + Arrayi(l, m, n));
            }
}
//=============================================================================================//
void InitializeBasicDataForAPackage::UpdateKernel::update(const size_t &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j, int k)
        {
            Vec3d position = index_handler_->DataPositionFromIndex(cell_index, Array3i(i, j, k));
            phi[i][j][k] = shape_->findSignedDistance(position);
            near_interface_id[i][j][k] = phi[i][j][k] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
SYCL_EXTERNAL void UpdateLevelSetGradient::UpdateKernel::update(const size_t &package_index)
{
    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = phi_gradient_[package_index];

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j, int k)
        {
            NeighbourIndex x1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + 1, j, k), neighborhood);
            NeighbourIndex x2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i - 1, j, k), neighborhood);
            NeighbourIndex y1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j + 1, k), neighborhood);
            NeighbourIndex y2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j - 1, k), neighborhood);
            NeighbourIndex z1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j, k + 1), neighborhood);
            NeighbourIndex z2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j, k - 1), neighborhood);
            Real dphidx = GET_NEIGHBOR_VAL(phi_, x1) - GET_NEIGHBOR_VAL(phi_, x2);
            Real dphidy = GET_NEIGHBOR_VAL(phi_, y1) - GET_NEIGHBOR_VAL(phi_, y2);
            Real dphidz = GET_NEIGHBOR_VAL(phi_, z1) - GET_NEIGHBOR_VAL(phi_, z2);
            pkg_data[i][j][k] = 0.5 * Vecd(dphidx, dphidy, dphidz) / data_spacing_;
        });
}
//=============================================================================================//
void UpdateKernelIntegrals::UpdateKernel::update(const size_t &package_index)
{
    auto &pkg_data_kernel_weight = kernel_weight_[package_index];
    auto &pkg_data_kernel_gradient = kernel_gradient_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;

    Real cut_cell_volume_fraction[7 * 7 * 7];
    Arrayi local_index[7 * 7 * 7];
    int idx = 0;
    mesh_for_each3d<-3, 4>(
        [&](int i, int j, int k)
        {
            NeighbourIndex neighbor_meta = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j, k), cell_neighborhood_[package_index]);
            Real phi_neighbor = GET_NEIGHBOR_VAL(phi_, neighbor_meta);
            if(phi_neighbor > -data_spacing_){
                Vecd phi_gradient = GET_NEIGHBOR_VAL(phi_gradient_, neighbor_meta);
                cut_cell_volume_fraction[idx] = CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                local_index[idx] = Arrayi(i, j, k);
                idx++;
            }
        });
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k){
            Vec3d position = index_handler_->DataPositionFromIndex(cell_index, Arrayi(i, j, k));
            std::pair<Real, Vecd> ret = computeKernelIntegral(position, package_index, Arrayi(i, j, k), cut_cell_volume_fraction, local_index, idx);
            pkg_data_kernel_weight[i][j][k] = ret.first;
            pkg_data_kernel_gradient[i][j][k] = ret.second;
        });
}
//=============================================================================================//
std::pair<Real, Vecd> UpdateKernelIntegrals::UpdateKernel::computeKernelIntegral(const Vecd &position,
                                                                                 const size_t &package_index,
                                                                                 const Arrayi &grid_index,
                                                                                 Real *cut_cell_volume_fraction,
                                                                                 Arrayi *local_index,
                                                                                 int n)
{
    Real phi = probe_signed_distance_.update(position);
    Real integral_kernel_weight(0);
    Vecd integral_kernel_gradient = Vecd::Zero();
    if (fabs(phi) < threshold)
    {
        for(int i = 0; i < n; i++){
            Vecd displacement = (grid_index - local_index[i]).cast<Real>().matrix() * data_spacing_;
            Real distance = displacement.norm();
            integral_kernel_weight += kernel_->W(global_h_ratio_, distance, displacement) * cut_cell_volume_fraction[i];
            integral_kernel_gradient += kernel_->dW(global_h_ratio_, distance, displacement) * cut_cell_volume_fraction[i] * displacement / (distance + TinyReal);
        }
    }

    std::pair<Real, Vecd> ret;
    ret.first = phi > threshold ? 1.0 : integral_kernel_weight * data_spacing_ * data_spacing_ * data_spacing_;
    ret.second = integral_kernel_gradient * data_spacing_ * data_spacing_ * data_spacing_;
    return ret;
}
//=============================================================================================//
void ReinitializeLevelSet::UpdateKernel::update(const size_t &package_index)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];
    auto &neighborhood = cell_neighborhood_[package_index];

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j, int k)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j][k] != 0)
            {
                Real phi_0 = phi_addrs[i][j][k];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                NeighbourIndex x1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + 1, j, k), neighborhood);
                NeighbourIndex x2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i - 1, j, k), neighborhood);
                NeighbourIndex y1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j + 1, k), neighborhood);
                NeighbourIndex y2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j - 1, k), neighborhood);
                NeighbourIndex z1 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j, k + 1), neighborhood);
                NeighbourIndex z2 = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i, j, k - 1), neighborhood);
                Real dv_x = upwindDifference(sign, GET_NEIGHBOR_VAL(phi_, x1) - phi_0, phi_0 - GET_NEIGHBOR_VAL(phi_, x2));
                Real dv_y = upwindDifference(sign, GET_NEIGHBOR_VAL(phi_, y1) - phi_0, phi_0 - GET_NEIGHBOR_VAL(phi_, y2));
                Real dv_z = upwindDifference(sign, GET_NEIGHBOR_VAL(phi_, z1) - phi_0, phi_0 - GET_NEIGHBOR_VAL(phi_, z2));
                phi_addrs[i][j][k] -= 0.3 * sign * (Vec3d(dv_x, dv_y, dv_z).norm() - data_spacing_);
            }
        });
}
//=============================================================================================//
void MarkNearInterface::UpdateKernel::update(const size_t &package_index, Real small_shift_factor)
{
    Real small_shift = data_spacing_ * small_shift_factor; 
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageDataMatrix<Real, 5> corner_averages;
    mesh_for_each3d<0, 5>(
        [&](int i, int j, int k)
        {
            corner_averages[i][j][k] = BaseMeshLocalDynamics::CornerAverage(phi_, Arrayi(i, j, k), Arrayi(-1, -1, -1), cell_neighborhood_[package_index]);
        });

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j, int k)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j][k];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < small_shift)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j][k];
                // find inner and outer cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
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
void RedistanceInterface::UpdateKernel::update(const size_t &package_index)
{
    BaseMeshLocalDynamics::for_each_cell_data(
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
                        NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + r, j + s, k + t), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id = GET_NEIGHBOR_VAL(near_interface_id_, neighbour_index);
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
                            NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
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
                            NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
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
void DiffuseLevelSetSign::UpdateKernel::update(const size_t &package_index)
{
    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j, int k)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_[package_index][i][j][k]) > 1)
            {
                mesh_find_if3d<-1, 2>(
                    [&](int l, int m, int n) -> bool
                    {
                        NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                        int near_interface_id = GET_NEIGHBOR_VAL(near_interface_id_, neighbour_index);
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_[package_index][i][j][k];
                            near_interface_id_[package_index][i][j][k] = near_interface_id;
                            phi_[package_index][i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
void WriteMeshFieldToPlt::update(std::ofstream &output_file)
{
    Arrayi number_of_operation = mesh_data_.global_mesh_.AllGridPoints();

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "z, "
                << "phi, "
                << "n_x, "
                << "n_y, "
                << "n_z "
                << "\n";
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
                << "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[0] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[1] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[2] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_.DataField(), Arrayi(i, j, k), &mesh_data_, cell_package_index_.DataField())
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_.DataField(), Arrayi(i, j, k), &mesh_data_, cell_package_index_.DataField())[0]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_.DataField(), Arrayi(i, j, k), &mesh_data_, cell_package_index_.DataField())[1]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_.DataField(), Arrayi(i, j, k), &mesh_data_, cell_package_index_.DataField())[2]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(near_interface_id_.DataField(), Arrayi(i, j, k), &mesh_data_, cell_package_index_.DataField())
                            << " ";
            }
            output_file << " \n";
        }
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//