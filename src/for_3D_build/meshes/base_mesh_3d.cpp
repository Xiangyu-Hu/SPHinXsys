#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
Arrayi Mesh::transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i) const
{
    size_t row_times_column_size = mesh_size[1] * mesh_size[2];
    size_t page = i / row_times_column_size;
    size_t left_over = (i - page * row_times_column_size);
    size_t row_size = mesh_size[2];
    size_t column = left_over / row_size;
    return Arrayi(page, column, left_over - column * row_size);
}
//=================================================================================================//
} // namespace SPH
