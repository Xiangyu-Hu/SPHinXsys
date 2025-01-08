#include "base_mesh.h"

namespace SPH
{
//=============================================================================================//
Arrayi Mesh::transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i) const
{
    size_t row_size = mesh_size[1];
    size_t column = i / row_size;
    return Arrayi(column, i - column * row_size);
}
} // namespace SPH
//=============================================================================================//
