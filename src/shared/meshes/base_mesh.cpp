#include "base_mesh.hpp"

namespace SPH
{
//=================================================================================================//
Mesh::Mesh(BoundingBoxd tentative_bounds, Real grid_spacing,
           UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset)
    : grid_spacing_(grid_spacing), buffer_width_(buffer_width),
      linear_cell_index_offset_(linear_cell_index_offset)
{
    Vecd mesh_buffer = Real(buffer_width) * grid_spacing * Vecd::Ones();
    mesh_lower_bound_ = tentative_bounds.lower_ - mesh_buffer;
    Vecd tentative_dimension = tentative_bounds.upper_ + mesh_buffer - mesh_lower_bound_;
    all_grid_points_ = ceil(tentative_dimension.array() / grid_spacing).cast<int>() + Arrayi::Ones();
    all_cells_ = all_grid_points_ - Arrayi::Ones();
}
//=================================================================================================//
Mesh::Mesh(Vecd mesh_lower_bound, Real grid_spacing, Arrayi all_grid_points)
    : mesh_lower_bound_{mesh_lower_bound}, grid_spacing_{grid_spacing},
      buffer_width_(0), all_grid_points_{all_grid_points},
      linear_cell_index_offset_(0)
{
    all_cells_ = all_grid_points_ - Arrayi::Ones();
}
//=================================================================================================//
Vecd Mesh::CellPositionFromIndex(const Arrayi &cell_index) const
{
    return CellLowerCornerPosition(cell_index) + 0.5 * Vecd::Ones() * grid_spacing_;
}
//=================================================================================================//
Vecd Mesh::GridPositionFromIndex(const Arrayi &grid_index) const
{
    return mesh_lower_bound_ + grid_index.cast<Real>().matrix() * grid_spacing_;
}
//=================================================================================================//
MultiLevelMeshField::MultiLevelMeshField(
    const std::string &name, BoundingBoxd tentative_bounds,
    Real Reference_grid_spacing, UnsignedInt buffer_width, size_t total_levels)
    : BaseMeshField(name), total_levels_(total_levels), total_number_of_cells_(0)
{
    for (size_t level = 0; level < total_levels; ++level)
    {
        meshes_.push_back(mesh_ptrs_keeper_.template createPtr<Mesh>(
            tentative_bounds, Reference_grid_spacing / math::pow(2.0, level),
            buffer_width, total_number_of_cells_));
        total_number_of_cells_ += meshes_.back()->NumberOfCells();
    };
}
//=================================================================================================//
void MultiLevelMeshField::writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence)
{
    for (UnsignedInt l = 0; l != meshes_.size(); ++l)
    {
        std::string full_file_name = partial_file_name + "_" + std::to_string(l) +
                                     "_" + std::to_string(sequence) + ".dat";
        std::ofstream out_file(full_file_name.c_str(), std::ios::app);
        writeCellVariableToPltByMesh(*meshes_[l], out_file);
        out_file.close();
    }
}
//=================================================================================================//
} // namespace SPH