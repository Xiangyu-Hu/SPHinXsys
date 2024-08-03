#include "base_mesh.h"

namespace SPH
{
//=============================================================================================//
Arrayi Mesh::transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i)
{
    size_t row_size = mesh_size[1];
    size_t column = i / row_size;
    return Arrayi(column, i - column * row_size);
}
//=============================================================================================//
size_t Mesh::transferMeshIndexTo1D(const Arrayi &mesh_size, const Arrayi &mesh_index)
{
    return mesh_index[0] * mesh_size[1] + mesh_index[1];
}
//=============================================================================================//
size_t Mesh::transferMeshIndexToMortonOrder(const Arrayi &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
}
} // namespace SPH
//=============================================================================================//
