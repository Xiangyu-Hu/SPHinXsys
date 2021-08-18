/**
 * @file 	base_mesh_supplementary.cpp
 * @author	Luhui Han, Chi ZHang Yongchuan YU and Xiangyu Hu
 */

#include "base_mesh.h"

namespace SPH {
	//=================================================================================================//
	void MeshIterator(const Vec3u &index_begin, const Vec3u &index_end, MeshFunctor& mesh_functor, Real dt)
	{
		for (size_t i = index_begin[0]; i != index_end[0]; ++i)
			for (size_t j = index_begin[1]; j != index_end[1]; ++j)
				for (size_t k = index_begin[2]; k != index_end[2]; ++k)
				{
					mesh_functor(Vecu(i, j, k), dt);
				}
	}
	//=================================================================================================//
	void MeshIterator_parallel(const Vecu &index_begin, const Vecu &index_end, MeshFunctor& mesh_functor, Real dt)
	{
		parallel_for(blocked_range3d<size_t>
			(index_begin[0], index_end[0], index_begin[1], index_end[1], index_begin[2], index_end[2]),
			[&](const blocked_range3d<size_t>& r) {
				for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
						for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
						{
							mesh_functor(Vecu(i, j, k), dt);
						}
			}, ap);
	}
	//=================================================================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(const Vecu &number_of_grid_points, size_t i)
	{
		size_t row_times_column_size = number_of_grid_points[1] * number_of_grid_points[2];
		size_t page = i / row_times_column_size;
		size_t left_over = (i - page * row_times_column_size);
		size_t row_size = number_of_grid_points[2];
		size_t column = left_over / row_size;
		return Vecu(page, column, left_over - column * row_size);
	}
	//=================================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(const Vecu &number_of_grid_points, const Vecu &grid_index)
	{
		return grid_index[0] * number_of_grid_points[1] * number_of_grid_points[2] 
			 + grid_index[1] * number_of_grid_points[2] 
			 + grid_index[2];
	}
    //=================================================================================================//
    size_t BaseMesh::transferMeshIndexToMortonOrder(const Vecu &grid_index)
    {
        return MortonCode(grid_index[0]) | (MortonCode(grid_index[1]) << 1)
                | (MortonCode(grid_index[2]) << 2);
    }
}
