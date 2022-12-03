#include "base_mesh.h"

namespace SPH
{
	//=================================================================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(const Vecu &number_of_mesh_indexes, size_t i)
	{
		size_t row_times_column_size = number_of_mesh_indexes[1] * number_of_mesh_indexes[2];
		size_t page = i / row_times_column_size;
		size_t left_over = (i - page * row_times_column_size);
		size_t row_size = number_of_mesh_indexes[2];
		size_t column = left_over / row_size;
		return Vecu(page, column, left_over - column * row_size);
	}
	//=================================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(const Vecu &number_of_mesh_indexes, const Vecu &mesh_index)
	{
		return mesh_index[0] * number_of_mesh_indexes[1] * number_of_mesh_indexes[2] + mesh_index[1] * number_of_mesh_indexes[2] + mesh_index[2];
	}
	//=================================================================================================//
	size_t BaseMesh::transferMeshIndexToMortonOrder(const Vecu &mesh_index)
	{
		return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1) | (MortonCode(mesh_index[2]) << 2);
	}
}
