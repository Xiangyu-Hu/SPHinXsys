#include "base_mesh.h"

namespace SPH {
	//=============================================================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(const Vecu &number_of_mesh_indexes, size_t i)
	{
		size_t row_size = number_of_mesh_indexes[1];
		size_t column = i / row_size;
		return Vecu(column, i - column * row_size);
	}
	//=============================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(const Vecu &number_of_mesh_indexes, const Vecu &mesh_index)
	{
		return mesh_index[0] * number_of_mesh_indexes[1] + mesh_index[1];
	}
    //=============================================================================================//
    size_t BaseMesh::transferMeshIndexToMortonOrder(const Vecu &mesh_index)
    {
        return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
    }
}
//=============================================================================================//
