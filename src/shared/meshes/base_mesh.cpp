#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
BaseMesh::BaseMesh(Arrayi all_grid_points)
    : all_grid_points_{all_grid_points} {};
//=================================================================================================//
BaseMesh::BaseMesh(Vecd mesh_lower_bound, Real grid_spacing, Arrayi all_grid_points)
    : mesh_lower_bound_{mesh_lower_bound}, grid_spacing_{grid_spacing},
      all_grid_points_{all_grid_points} {}
//=================================================================================================//
BaseMesh::BaseMesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
    : BaseMesh()
{
    grid_spacing_ = grid_spacing;
    Vecd mesh_buffer = Real(buffer_width) * grid_spacing * Vecd::Ones();
    mesh_lower_bound_ = tentative_bounds.first_ - mesh_buffer;
    Vecd tentative_dimension = tentative_bounds.second_ + mesh_buffer - mesh_lower_bound_;
    all_grid_points_ = ceil(tentative_dimension.array() / grid_spacing).cast<int>() + Arrayi::Ones();
}
//=================================================================================================//
Arrayi BaseMesh::CellIndexFromPosition(const Vecd &position)
{
    return floor((position - mesh_lower_bound_).array() / grid_spacing_)
        .cast<int>()
        .max(Arrayi::Zero())
        .min(all_grid_points_ - 2 * Arrayi::Ones());
}
//=================================================================================================//
Vecd BaseMesh::CellPositionFromIndex(const Arrayi &cell_index)
{
    return mesh_lower_bound_ + (cell_index.cast<Real>().matrix() + 0.5 * Vecd::Ones()) * grid_spacing_;
}
//=================================================================================================//
Vecd BaseMesh::GridPositionFromIndex(const Arrayi &grid_index)
{
    return mesh_lower_bound_ + grid_index.cast<Real>().matrix() * grid_spacing_;
}
//=================================================================================================//
size_t BaseMesh::MortonCode(const size_t &i)
{
    size_t x = i;
    x &= 0x3ff;
    x = (x | x << 16) & 0x30000ff;
    x = (x | x << 8) & 0x300f00f;
    x = (x | x << 4) & 0x30c30c3;
    x = (x | x << 2) & 0x9249249;
    return x;
}
//=================================================================================================//
Mesh::Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
    : BaseMesh(tentative_bounds, grid_spacing, buffer_width),
      all_cells_{this->AllCellsFromAllGridPoints(this->AllGridPoints())},
      buffer_width_{buffer_width} {}
//=================================================================================================//
Mesh::Mesh(Vecd mesh_lower_bound, Arrayi all_cells, Real grid_spacing)
    : BaseMesh(), all_cells_{all_cells}
{
    mesh_lower_bound_ = mesh_lower_bound;
    grid_spacing_ = grid_spacing;
    all_grid_points_ = AllGridPointsFromAllCells(all_cells_);
}
//=================================================================================================//
void Mesh::copyMeshProperties(Mesh *another_mesh)
{
    mesh_lower_bound_ = another_mesh->mesh_lower_bound_;
    grid_spacing_ = another_mesh->grid_spacing_;
    all_grid_points_ = another_mesh->all_grid_points_;
    all_cells_ = another_mesh->all_cells_;
    buffer_width_ = another_mesh->buffer_width_;
}
//=================================================================================================//
} // namespace SPH
