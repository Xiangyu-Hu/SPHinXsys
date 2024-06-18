#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
Array3i BaseMesh::transfer1DtoMeshIndex(const Array3i &mesh_size, size_t i)
{
    size_t row_times_column_size = mesh_size[1] * mesh_size[2];
    size_t page = i / row_times_column_size;
    size_t left_over = (i - page * row_times_column_size);
    size_t row_size = mesh_size[2];
    size_t column = left_over / row_size;
    return Array3i(page, column, left_over - column * row_size);
}
//=================================================================================================//
size_t BaseMesh::transferMeshIndexTo1D(const Array3i &mesh_size, const Array3i &mesh_index)
{
    return mesh_index[0] * mesh_size[1] * mesh_size[2] +
           mesh_index[1] * mesh_size[2] +
           mesh_index[2];
}
//=================================================================================================//
size_t BaseMesh::transferMeshIndexToMortonOrder(const Array3i &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1) | (MortonCode(mesh_index[2]) << 2);
}
} // namespace SPH
