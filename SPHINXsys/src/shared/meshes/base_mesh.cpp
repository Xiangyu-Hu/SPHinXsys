#include "base_mesh.h"

namespace SPH
{
	//=================================================================================================//
	BaseMesh::BaseMesh(Vecu number_of_grid_points)
		: number_of_grid_points_{number_of_grid_points}
	{};
	//=================================================================================================//
	BaseMesh::BaseMesh(Vecd mesh_lower_bound, Real grid_spacing, Vecu number_of_grid_points)
		: mesh_lower_bound_{mesh_lower_bound}
		, grid_spacing_{grid_spacing}
		, number_of_grid_points_{number_of_grid_points}
	{}
	//=================================================================================================//
	BaseMesh::BaseMesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width) 
		: BaseMesh()
	{
		grid_spacing_ = grid_spacing;
		Vecd mesh_buffer = Real(buffer_width) * grid_spacing * Vecd::Ones();
		mesh_lower_bound_ = tentative_bounds.first_ - mesh_buffer;
		Vecd tentative_upper_bound = tentative_bounds.second_ + mesh_buffer;
		for (int i = 0; i != Dimensions; ++i)
		{
			number_of_grid_points_[i] = 
				1 + static_cast<int>(ceil((tentative_upper_bound[i] - mesh_lower_bound_[i]) / grid_spacing));
		}
	}
	//=================================================================================================//
	Vecu BaseMesh::CellIndexFromPosition(const Vecd &position)
	{
		Vecd rltpos = position - mesh_lower_bound_;
		Vecu cell_index = Vecu::Zero();
		for (int n = 0; n < rltpos.size(); n++)
		{
			cell_index[n] =
				clamp((int)floor(rltpos[n] / grid_spacing_), 0, int(number_of_grid_points_[n]) - 2);
		}
		return cell_index;
	}
	//=================================================================================================//
	Vecd BaseMesh::CellPositionFromIndex(const Vecu &cell_index)
	{
		Vecd cell_position;
		for (int n = 0; n < cell_position.size(); n++)
		{
			cell_position[n] = mesh_lower_bound_[n] + (Real(cell_index[n]) + 0.5) * grid_spacing_;
		}
		return cell_position;
	}
	//=================================================================================================//
	Vecd BaseMesh::GridPositionFromIndex(const Vecu &grid_index)
	{
		Vecd grid_position;
		for (int n = 0; n < grid_position.size(); n++)
		{
			grid_position[n] = mesh_lower_bound_[n] + Real(grid_index[n]) * grid_spacing_;
		}
		return grid_position;
	}
	//=================================================================================================//
	size_t BaseMesh::MortonCode(const size_t &i)
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
	Mesh::Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
		: BaseMesh(tentative_bounds, grid_spacing, buffer_width)
		, number_of_cells_{this->NumberOfCellsFromNumberOfGridPoints(this->NumberOfGridPoints())}
		, buffer_width_{buffer_width}
	{}
	//=================================================================================================//
	Mesh::Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real grid_spacing)
		: BaseMesh()
		, number_of_cells_{number_of_cells}
	{
		mesh_lower_bound_ = mesh_lower_bound;
		grid_spacing_ = grid_spacing;
		number_of_grid_points_ = NumberOfGridPointsFromNumberOfCells(number_of_cells_);
	}
	//=================================================================================================//
	void Mesh::copyMeshProperties(Mesh *another_mesh)
	{
		mesh_lower_bound_ = another_mesh->mesh_lower_bound_;
		grid_spacing_ = another_mesh->grid_spacing_;
		number_of_grid_points_ = another_mesh->number_of_grid_points_;
		number_of_cells_ = another_mesh->number_of_cells_;
		buffer_width_ = another_mesh->buffer_width_;
	}
	//=================================================================================================//
}
