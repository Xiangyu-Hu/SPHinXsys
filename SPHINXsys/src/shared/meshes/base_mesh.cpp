/**
 * @file 	base_mesh.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_mesh.h"

namespace SPH {

	Mesh::Mesh(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: grid_spacing_(grid_spacing), buffer_size_(buffer_size)
	{
		Vecd mesh_buffer = Vecd(Real(buffer_size_)*grid_spacing_);
		mesh_lower_bound_ = lower_bound - mesh_buffer;
		mesh_upper_bound_ = upper_bound + mesh_buffer;

		number_of_grid_points_ 	= CalcNumberOfGridPoints();
	}

	Vecu Mesh::CalcNumberOfGridPoints()
	{
		Vecu number_of_grid_points;
		Vecd zero(0);
		for (int i = 0; i < zero.size(); ++i) {
			number_of_grid_points[i] = static_cast<int>(ceil((mesh_upper_bound_[i]
				- mesh_lower_bound_[i]) / grid_spacing_));
		}

		return number_of_grid_points;
	}
//===================================================================//
	Vecu Mesh::GridIndexesFromPosition(Vecd position)
	{
		Vecd rltpos = position - mesh_lower_bound_;
		Vecu grid_pos(0);
		for (int n = 0; n < rltpos.size(); n++)
		{
			grid_pos[n] = clamp((int)floor(rltpos[n] / grid_spacing_),
				0, int(number_of_grid_points_[n]) - 1);
		}
		return grid_pos;
	}
//===================================================================//
	Vecd Mesh::GridPositionFromIndexes(Vecu grid_indexes)
	{
		Vecd grid_position;
		for (int n = 0; n < grid_position.size(); n++)
		{
			grid_position[n] = mesh_lower_bound_[n] 
				+ Real(grid_indexes[n])*grid_spacing_;
		}
		return grid_position;
	}

	Vecd Mesh::CellPositionFromIndexes(Vecu cell_indexes)
	{
		Vecd cell_position;
		for (int n = 0; n < cell_position.size(); n++)
		{
			cell_position[n] = mesh_lower_bound_[n]
				+ (Real(cell_indexes[n]) + 0.5)*grid_spacing_;
		}
		return cell_position;
	}

//===================================================================//
	BackgroundData
		::BackgroundData(Real level_set, Vecd normal_direction)
		: phi_(level_set), n_(normal_direction), kappa_(0.0)
	{

	}
//===================================================================//
	MeshBackground
		::MeshBackground(Vecd lower_bound, Vecd upper_bound, 
			Real grid_spacing, size_t buffer_size)
		: Mesh(lower_bound, upper_bound, grid_spacing, buffer_size)
	{

	}
}