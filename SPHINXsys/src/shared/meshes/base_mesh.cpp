/**
 * @file 	base_mesh.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "base_mesh.h"

namespace SPH {
	//=================================================================================================//
	BaseMesh::BaseMesh() : mesh_lower_bound_(0), grid_spacing_(1.0),
		number_of_grid_points_(0) {};
	//=================================================================================================//
	BaseMesh::BaseMesh(Vecu number_of_grid_points) : 
		mesh_lower_bound_(0), grid_spacing_(1.0),
		number_of_grid_points_(number_of_grid_points) {};
	//=================================================================================================//
	Vecu BaseMesh::GridIndexFromPosition(const Vecd& position)
	{
		Vecd rltpos = position - mesh_lower_bound_;
		Vecu gird_pos(0);
		for (int n = 0; n < rltpos.size(); n++)
		{
			gird_pos[n] = clamp((int)floor(rltpos[n] / grid_spacing_),
				0, int(number_of_grid_points_[n]) - 1);
		}
		return gird_pos;
	}
	//=================================================================================================//
	Vecd BaseMesh::GridPositionFromIndex(Vecu grid_index)
	{
		Vecd grid_position;
		for (int n = 0; n < grid_position.size(); n++)
		{
			grid_position[n] = mesh_lower_bound_[n]
				+ Real(grid_index[n]) * grid_spacing_;
		}
		return grid_position;
	}
	//=================================================================================================//
	size_t BaseMesh::MortonCode(const size_t& i)
	{
		size_t x = i;
		x &= 0x3ff;
		x = (x | x << 16) & 0x30000ff;
		x = (x | x << 8) & 0x300f00f;
		x = (x | x << 4) & 0x30c30c3;
		x = (x | x << 2) & 0x9249249;
		return x;
	}
	//=================================================================================================//
	Mesh::Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width) : 
		BaseMesh(), name_("Mesh"), buffer_width_(buffer_width), number_of_cells_(0)
	{
		initializeWithBoundingBox(tentative_bounds, grid_spacing, buffer_width);
	}
	//=================================================================================================//
	Mesh::Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real grid_spacing) : 
		BaseMesh(), name_("Mesh"), buffer_width_(0), number_of_cells_(number_of_cells)
	{
		mesh_lower_bound_ = mesh_lower_bound;
		grid_spacing_ = grid_spacing;
		number_of_grid_points_ = NumberOfGridPointsFromNumberOfCells(number_of_cells_);
	}
	//=================================================================================================//
	void Mesh::
		initializeWithBoundingBox(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
	{
		grid_spacing_ = grid_spacing;
		Vecd mesh_buffer = Vecd(Real(buffer_width) * grid_spacing);
		mesh_lower_bound_ = tentative_bounds.first - mesh_buffer;
		Vecd tentative_upper_bound = tentative_bounds.second + mesh_buffer;

		for (int i = 0; i != Dimensions; ++i) {
			number_of_cells_[i] = 
				static_cast<int>(ceil((tentative_upper_bound[i] - mesh_lower_bound_[i]) / grid_spacing));
		}
		number_of_grid_points_ = number_of_cells_ + Vecu(1);
	}
	//=================================================================================================//
	void Mesh::copyMeshProperties(Mesh* another_mesh)
	{
		mesh_lower_bound_ = another_mesh->mesh_lower_bound_;
		grid_spacing_ = another_mesh->grid_spacing_;
		number_of_grid_points_ = another_mesh->number_of_grid_points_;
		number_of_cells_ = another_mesh->number_of_cells_;
		buffer_width_ = another_mesh->buffer_width_;
	}
	//=================================================================================================//
	Vecu Mesh::CellIndexFromPosition(Vecd& position)
	{
		Vecd rltpos(0);
		Vecu cell_index(0);
		for (int n = 0; n < rltpos.size(); n++)
		{
			rltpos[n] = position[n] - mesh_lower_bound_[n] - 0.5 * grid_spacing_;
			cell_index[n] = clamp((int)floor(rltpos[n] / grid_spacing_),
				0, int(number_of_cells_[n]) - 1);
		}
		return cell_index;
	}
	//=================================================================================================//
	Vecd Mesh::CellPositionFromIndex(Vecu cell_index)
	{
		Vecd cell_position;
		for (int n = 0; n < cell_position.size(); n++)
		{
			cell_position[n] = mesh_lower_bound_[n] + (Real(cell_index[n]) + 0.5)* grid_spacing_;
		}
		return cell_position;
	}
	//=================================================================================================//
	bool Mesh::isWithinMeshBound(Vecd position)
	{
		bool is_bounded = true;
		Vecu cell_pos = CellIndexFromPosition(position);
		for (int i = 0; i != position.size(); ++i) {
			if (cell_pos[i] < 2) is_bounded = false;
			if (cell_pos[i] > (number_of_cells_[i] - 2)) is_bounded = false;
		}
		return is_bounded;
	}
	//=================================================================================================//
}
