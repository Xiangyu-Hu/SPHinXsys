#include "base_mesh.h"

namespace SPH
{
//=============================================================================================//
Array2i BaseMesh::transfer1DtoMeshIndex(const Array2i &mesh_size, size_t i)
{
    size_t row_size = mesh_size[1];
    size_t column = i / row_size;
    return Array2i(column, i - column * row_size);
}
//=============================================================================================//
size_t BaseMesh::transferMeshIndexTo1D(const Array2i &mesh_size, const Array2i &mesh_index)
{
    return mesh_index[0] * mesh_size[1] + mesh_index[1];
}
//=============================================================================================//
size_t BaseMesh::transferMeshIndexToMortonOrder(const Array2i &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
}
} // namespace SPH
//=============================================================================================//
