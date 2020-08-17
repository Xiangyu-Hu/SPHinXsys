/**
 * @file 	base_mesh.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_mesh.h"
#include "base_mesh.hpp"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	Vecu BaseMesh::GridIndexFromPosition(Vecd& position)
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
	Mesh::Mesh(Vecd lower_bound, Vecd upper_bound, Real grid_spacing,
		size_t buffer_size) : BaseMesh (), buffer_size_(buffer_size)

	{
		grid_spacing_ = grid_spacing;
		cell_spacing_ = grid_spacing;
		setMeshLowerBound(lower_bound, grid_spacing, buffer_size);
		number_of_cells_ = calcNumberOfCells(lower_bound, upper_bound, grid_spacing, buffer_size);
		number_of_grid_points_ = NumberOfGridPointsFromNumberOfCells(number_of_cells_);
	}
	//=================================================================================================//
	Mesh::Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: BaseMesh(), number_of_cells_(number_of_cells), cell_spacing_(cell_spacing), buffer_size_(0)
	{
		mesh_lower_bound_ = mesh_lower_bound;
		grid_spacing_ = cell_spacing;
		number_of_grid_points_ = NumberOfGridPointsFromNumberOfCells(number_of_cells_);
	}
	//=================================================================================================//
	void Mesh::setMeshLowerBound(Vecd lower_bound, Real grid_spacing, size_t buffer_size)
	{
		Vecd mesh_buffer = Vecd(Real(buffer_size) * grid_spacing);
		mesh_lower_bound_ = lower_bound - mesh_buffer;
	}
	//=================================================================================================//
	Vecu Mesh::calcNumberOfCells(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size)
	{
		Vecd mesh_buffer = Vecd(Real(buffer_size) * grid_spacing);
		Vecd mesh_lower_bound = lower_bound - mesh_buffer;
		Vecd tentative_upper_bound = upper_bound + mesh_buffer;

		Vecu number_of_cells(0);
		Vecd zero(0);
		for (int i = 0; i < zero.size(); ++i) {
			number_of_cells[i] = static_cast<int>(ceil((tentative_upper_bound[i]
				- mesh_lower_bound[i]) / grid_spacing));
		}

		return number_of_cells;
	}
	//=================================================================================================//
	void Mesh::copyMeshProperties(Mesh* another_mesh)
	{
		mesh_lower_bound_ = another_mesh->mesh_lower_bound_;
		grid_spacing_ = another_mesh->grid_spacing_;
		cell_spacing_ = another_mesh->cell_spacing_;
		number_of_cells_ = another_mesh->number_of_cells_;
		number_of_grid_points_ = another_mesh->number_of_grid_points_;
		buffer_size_ = another_mesh->buffer_size_;
	}
	//=================================================================================================//
	Vecu Mesh::CellIndexesFromPosition(Vecd& position)
	{
		Vecd rltpos = position - mesh_lower_bound_;
		Vecu cell_pos(0);
		for (int n = 0; n < rltpos.size(); n++)
		{
			cell_pos[n] = clamp((int)floor(rltpos[n] / cell_spacing_),
				0, int(number_of_cells_[n]) - 1);
		}
		return cell_pos;
	}
	//=================================================================================================//
	Vecd Mesh::CellPositionFromIndexes(Vecu cell_indexes)
	{
		Vecd cell_position;
		for (int n = 0; n < cell_position.size(); n++)
		{
			cell_position[n] = mesh_lower_bound_[n]
				+ (Real(cell_indexes[n]) + 0.5)* cell_spacing_;
		}
		return cell_position;
	}
	//=================================================================================================//
	bool Mesh::isWithinMeshBound(Vecd position)
	{
		bool is_bounded = true;
		Vecu cell_pos = CellIndexesFromPosition(position);
		for (int i = 0; i != position.size(); ++i) {
			if (cell_pos[i] < 2) is_bounded = false;
			if (cell_pos[i] > (number_of_cells_[i] - 2)) is_bounded = false;
		}
		return is_bounded;
	}
	//=================================================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd lower_bound,
			Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: BaseMeshWithDataPackages(lower_bound, upper_bound,
			grid_spacing, buffer_size), sph_body_(NULL) {}
	//=================================================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: BaseMeshWithDataPackages(mesh_lower_bound, number_of_cells,
			cell_spacing), sph_body_(NULL) {}
	//=================================================================================================//
	LevelSet
		::LevelSet(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>(lower_bound,
			upper_bound, grid_spacing, buffer_size)
	{
		Real far_field_distance = cell_spacing_ * (Real)buffer_size;
		LevelSetDataPackage* negative_far_field = new LevelSetDataPackage();
		negative_far_field->initializeWithUniformData(-far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(negative_far_field);
		LevelSetDataPackage* positive_far_field = new LevelSetDataPackage();
		positive_far_field->initializeWithUniformData(far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(positive_far_field);
	}
	//=================================================================================================//
	LevelSet
		::LevelSet(SPHBody* sph_body, Vecd lower_bound,
			Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: LevelSet(lower_bound,	upper_bound, grid_spacing, buffer_size)
	{
		sph_body_ = sph_body;
		initializeDataPackages();
	}
	//=================================================================================================//
	LevelSet
		::LevelSet(Vecd mesh_lower_bound, 
			Vecu number_of_cells, Real cell_spacing)
		: MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>(mesh_lower_bound,
			number_of_cells, cell_spacing) {}
	//=================================================================================================//
	void LevelSet::initializeDataPackages()
	{
		MeshFunctor initialize_data_in_a_cell = std::bind(&LevelSet::initializeDataInACell, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, initialize_data_in_a_cell);
		MeshFunctor tag_a_cell_inner_pkg = std::bind(&LevelSet::tagACellIsInnerPackage, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, tag_a_cell_inner_pkg);
		MeshFunctor initial_address_in_a_cell = std::bind(&LevelSet::initializeAddressesInACell, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, initial_address_in_a_cell);
		updateNormalDirection();
	}
	//=================================================================================================//
	void LevelSet::initializeAddressesInACell(Vecu cell_index, Real dt)
	{
		initializePackageAddressesInACell<LevelSetData>(cell_index);
	}
	//=================================================================================================//
	void LevelSet::updateNormalDirection()
	{
		PackageFunctor<void, LevelSetDataPackage> compute_normalized_gradient
			= std::bind(&LevelSet::template computeNormalizedGradient<LevelSetData, &LevelSetData::phi_, &LevelSetData::n_>, this, _1, _2);
		PackageIterator_parallel<LevelSetDataPackage>(inner_data_pkgs_, compute_normalized_gradient);
	}
	//=================================================================================================//
	Vecd LevelSet::probeNormalDirection(Vecd position)
	{
		return probeMesh<Vecd, LevelSetData, &LevelSetData::n_>(position);
	}
	//=================================================================================================//
	Real LevelSet::probeLevelSet(Vecd position)
	{
		return probeMesh<Real, LevelSetData, &LevelSetData::phi_>(position);
	}
	//=================================================================================================//
	MultiresolutionLevelSet
		::MultiresolutionLevelSet(SPHBody* sph_body, Vecd lower_bound, Vecd upper_bound,
		Real reference_cell_spacing, size_t total_levels, size_t buffer_size)
		: MultilevelMesh<BaseLevelSet, LevelSet>(lower_bound, upper_bound,
			reference_cell_spacing, total_levels, buffer_size)
	{
		sph_body_ = sph_body;
		for (size_t level = 0; level != total_levels_; ++level) {
			mesh_levels_[level]->setSPHBody(sph_body);
			mesh_levels_[level]->initializeDataPackages();
		}
	}
	//=================================================================================================//
}
