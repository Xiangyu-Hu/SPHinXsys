#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
Mesh::Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
    : grid_spacing_(grid_spacing), buffer_width_(buffer_width)
{
    Vecd mesh_buffer = Real(buffer_width) * grid_spacing * Vecd::Ones();
    mesh_lower_bound_ = tentative_bounds.first_ - mesh_buffer;
    Vecd tentative_dimension = tentative_bounds.second_ + mesh_buffer - mesh_lower_bound_;
    all_grid_points_ = ceil(tentative_dimension.array() / grid_spacing).cast<int>() + Arrayi::Ones();
    all_cells_ = all_grid_points_ - Arrayi::Ones();
}
//=================================================================================================//
Mesh::Mesh(Vecd mesh_lower_bound, Real grid_spacing, Arrayi all_grid_points)
    : mesh_lower_bound_{mesh_lower_bound}, grid_spacing_{grid_spacing},
      buffer_width_(0), all_grid_points_{all_grid_points}
{
    all_cells_ = all_grid_points_ - Arrayi::Ones();
}
//=================================================================================================//
Vecd Mesh::CellLowerCornerPosition(const Arrayi &cell_index) const
{
    return mesh_lower_bound_ + cell_index.cast<Real>().matrix() * grid_spacing_;
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
} // namespace SPH