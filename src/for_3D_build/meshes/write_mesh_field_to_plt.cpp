#include "write_mesh_field_to_plt.h"

namespace SPH
{
//=============================================================================================//
std::vector<Block3D> clusterActiveCells3D(const std::vector<Coord3D> &activeCells)
{
    std::unordered_set<Coord3D, Coord3DHash> activeSet(activeCells.begin(), activeCells.end());
    std::vector<Block3D> blocks;

    while (!activeSet.empty())
    {
        Coord3D start = *activeSet.begin();
        int x0 = start[0], y0 = start[1], z0 = start[2];
        int x1 = x0;

        // Extend in x
        while (activeSet.count({x1 + 1, y0, z0}))
            ++x1;

        // Extend in y
        int y1 = y0;
        bool extendY = true;
        while (extendY)
        {
            for (int x = x0; x <= x1; ++x)
            {
                if (!activeSet.count({x, y1 + 1, z0}))
                {
                    extendY = false;
                    break;
                }
            }
            if (extendY)
                ++y1;
        }

        // Extend in z
        int z1 = z0;
        bool extendZ = true;
        while (extendZ)
        {
            for (int y = y0; y <= y1; ++y)
            {
                for (int x = x0; x <= x1; ++x)
                {
                    if (!activeSet.count({x, y, z1 + 1}))
                    {
                        extendZ = false;
                        break;
                    }
                }
                if (!extendZ)
                    break;
            }
            if (extendZ)
                ++z1;
        }

        // Remove covered cells
        for (int z = z0; z <= z1; ++z)
        {
            for (int y = y0; y <= y1; ++y)
            {
                for (int x = x0; x <= x1; ++x)
                {
                    activeSet.erase({x, y, z});
                }
            }
        }

        blocks.push_back({{x0, y0, z0}, {x1, y1, z1}});
    }

    return blocks;
}
//=============================================================================================//
void WriteMeshFieldToPlt::update(std::ofstream &output_file)
{
    StdVec<Coord3D> active_cells;
    mesh_for_each(Array3i::Zero(), mesh_data_.AllCells(), [&](const Arrayi &cell_index)
                  { 
      if (mesh_data_.isCoreDataPackage(cell_index))
      {
        active_cells.push_back({cell_index[0], cell_index[1], cell_index[2]});
      } });

    StdVec<Block3D> clustered_blocks = clusterActiveCells3D(active_cells);

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
                << "n_z, "
                << "near_interface_id "
                << "\n";
    for (size_t l = 0; l != clustered_blocks.size(); ++l)
    {
        auto &block = clustered_blocks[l];
        Array3i lower_bound_cell_index = Array3i(block.first[0], block.first[1], block.first[2]);
        Array3i block_size = Array3i(block.second[0] - block.first[0] + 1,
                                     block.second[1] - block.first[1] + 1,
                                     block.second[2] - block.first[2] + 1);
        Array3i global_lower_bound = mesh_data_.DataPackageSize() * lower_bound_cell_index - Array3i::Ones();
        Array3i global_upper_bound = mesh_data_.DataPackageSize() * (lower_bound_cell_index + block_size) + Array3i::Ones();
        output_file << "zone i=" << block_size[0] * mesh_data_.DataPackageSize() + 2 << "  j="
                    << block_size[1] * mesh_data_.DataPackageSize() + 2 << "  k="
                    << block_size[2] * mesh_data_.DataPackageSize() + 2
                    << "  DATAPACKING=POINT \n";

        mesh_for_column_major(
            global_lower_bound, global_upper_bound,
            [&](const Array3i &global_index)
            {
                Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(global_index);
                output_file << data_position[0] << " ";
                output_file << data_position[1] << " ";
                output_file << data_position[2] << " ";
                output_file << mesh_data_.DataValueFromGlobalIndex(phi_, global_index) << " ";
                output_file << mesh_data_.DataValueFromGlobalIndex(phi_gradient_, global_index)[0] << " ";
                output_file << mesh_data_.DataValueFromGlobalIndex(phi_gradient_, global_index)[1] << " ";
                output_file << mesh_data_.DataValueFromGlobalIndex(phi_gradient_, global_index)[2] << " ";
                output_file << mesh_data_.DataValueFromGlobalIndex(near_interface_id_, global_index) << " ";
                output_file << " \n";
            });
        output_file << " \n";
    }
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//