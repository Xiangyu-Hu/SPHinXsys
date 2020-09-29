/**
 * @file 	base_mesh_supplementary.cpp
 * @author	Luhui Han, Chi ZHang Yongchuan YU and Xiangyu Hu
 * @version	0.1
 */

#include "base_mesh.h"
#include "base_mesh.hpp"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "neighbor_relation.h"
#include "base_data_package.h"

namespace SPH {
	//=================================================================================================//
	void MeshIterator(Vec3u index_begin, Vec3u index_end, MeshFunctor& mesh_functor, Real dt)
	{
		for (size_t i = index_begin[0]; i != index_end[0]; ++i)
			for (size_t j = index_begin[1]; j != index_end[1]; ++j)
				for (size_t k = index_begin[2]; k != index_end[2]; ++k)
				{
					mesh_functor(Vecu(i, j, k), dt);
				}
	}
	//=================================================================================================//
	void MeshIterator_parallel(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt)
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
	Vecu BaseMesh::transfer1DtoMeshIndex(Vecu number_of_grid_points, size_t i)
	{
		size_t row_times_column_size = number_of_grid_points[1] * number_of_grid_points[2];
		size_t page = i / row_times_column_size;
		size_t left_over = (i - page * row_times_column_size);
		size_t row_size = number_of_grid_points[2];
		size_t column = left_over / row_size;
		return Vecu(page, column, left_over - column * row_size);
	}
	//=================================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(Vecu number_of_grid_points, Vecu grid_index)
	{
		return grid_index[0] * number_of_grid_points[1] * number_of_grid_points[2] 
			 + grid_index[1] * number_of_grid_points[2] 
			 + grid_index[2];
	}
	//=================================================================================================//
}
