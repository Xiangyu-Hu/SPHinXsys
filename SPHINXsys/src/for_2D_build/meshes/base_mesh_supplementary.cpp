#include "base_mesh.h"

namespace SPH {
	//=============================================================================================//
	Arrayi BaseMesh::transfer1DtoMeshIndex(const Arrayi &number_of_mesh_indexes, size_t i)
	{
		size_t row_size = number_of_mesh_indexes[1];
		size_t column = i / row_size;
		return Arrayi(column, i - column * row_size);
	}
	//=============================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(const Arrayi &number_of_mesh_indexes, const Arrayi &mesh_index)
	{
		return mesh_index[0] * number_of_mesh_indexes[1] + mesh_index[1];
	}
    //=============================================================================================//
    size_t BaseMesh::transferMeshIndexToMortonOrder(const Arrayi &mesh_index)
    {
        return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
    }
}
//=============================================================================================//
