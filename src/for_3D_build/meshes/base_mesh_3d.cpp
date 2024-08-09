#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
Arrayi Mesh::transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i)
{
    size_t row_times_column_size = mesh_size[1] * mesh_size[2];
    size_t page = i / row_times_column_size;
    size_t left_over = (i - page * row_times_column_size);
    size_t row_size = mesh_size[2];
    size_t column = left_over / row_size;
    return Arrayi(page, column, left_over - column * row_size);
}
//=================================================================================================//
size_t Mesh::transferMeshIndexTo1D(const Arrayi &mesh_size, const Arrayi &mesh_index)
{
    return mesh_index[0] * mesh_size[1] * mesh_size[2] +
           mesh_index[1] * mesh_size[2] +
           mesh_index[2];
}
//=================================================================================================//
size_t Mesh::transferMeshIndexToMortonOrder(const Arrayi &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1) | (MortonCode(mesh_index[2]) << 2);
}
} // namespace SPH
