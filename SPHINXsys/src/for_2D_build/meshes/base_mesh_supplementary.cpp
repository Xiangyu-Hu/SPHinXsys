/**
 * @file 	base_mesh_supplementary.cpp
 * @author	Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
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

#include "math.h"
//=================================================================================================//
namespace SPH {
	//=============================================================================================//
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt)
	{
		for (size_t i = index_begin[0]; i != index_end[0]; ++i)
			for (size_t j = index_begin[1]; j != index_end[1]; ++j) {
				mesh_functor(Vecu(i, j), dt);
			}
	}
	//=============================================================================================//
	void MeshIterator_parallel(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt)
	{
		parallel_for(blocked_range2d<size_t>
			(index_begin[0], index_end[0], index_begin[1], index_end[1]),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
					{
						mesh_functor(Vecu(i, j), dt);
					}
			}, ap);
	}
	//=============================================================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(Vecu number_of_grid_points, size_t i)
	{
		size_t row_size = number_of_grid_points[1];
		size_t column = i / row_size;
		return Vec2u(column, i - column * row_size);
	}
	//=============================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(Vecu number_of_grid_points, Vecu grid_index)
	{
		return grid_index[0] * number_of_grid_points[1] + grid_index[1];
	}
	//=============================================================================================//
}
//=============================================================================================//
