/**
 * @file 	base_mesh.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_mesh.h"
#include "base_mesh.hpp"
#include "base_body.h"

namespace SPH {
	//===================================================================//
	Vecu BaseMesh::GridIndexesFromPosition(Vecd& position)
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
	//===================================================================//
	Vecd BaseMesh::GridPositionFromIndexes(Vecu grid_indexes)
	{
		Vecd grid_position;
		for (int n = 0; n < grid_position.size(); n++)
		{
			grid_position[n] = mesh_lower_bound_[n]
				+ Real(grid_indexes[n]) * grid_spacing_;
		}
		return grid_position;
	}
	//===================================================================//
	Mesh::Mesh(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, 
		size_t buffer_size) : BaseMesh () 

	{
		buffer_size_ = buffer_size;
		grid_spacing_ = grid_spacing;
		cell_spacing_ = grid_spacing;
		setMeshLowerBound(lower_bound, grid_spacing, buffer_size);
		number_of_cells_ = calcNumberOfCells(lower_bound, upper_bound, grid_spacing, buffer_size);
		number_of_grid_points_ = getNumberOfGridPoints(number_of_cells_);
	}
	//===================================================================//
	Mesh::Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: BaseMesh(), number_of_cells_(number_of_cells), cell_spacing_(cell_spacing)
	{
		mesh_lower_bound_ = mesh_lower_bound;
		grid_spacing_ = cell_spacing;
		buffer_size_ = 0;
		number_of_grid_points_ = getNumberOfGridPoints(number_of_cells_);
	}
	void Mesh::setMeshLowerBound(Vecd lower_bound, Real grid_spacing, size_t buffer_size)
	{
		Vecd mesh_buffer = Vecd(Real(buffer_size) * grid_spacing);
		mesh_lower_bound_ = lower_bound - mesh_buffer;
	}
	//===================================================================//
	Vecu Mesh::calcNumberOfCells(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size)
	{
		Vecd mesh_buffer = Vecd(Real(buffer_size) * grid_spacing);
		Vecd mesh_lower_bound = lower_bound - mesh_buffer;
		Vecd tentive_upper_bound = upper_bound + mesh_buffer;

		Vecu number_of_cells(0);
		Vecd zero(0);
		for (int i = 0; i < zero.size(); ++i) {
			number_of_cells[i] = static_cast<int>(ceil((tentive_upper_bound[i]
				- mesh_lower_bound[i]) / grid_spacing));
		}

		return number_of_cells;
	}
	void Mesh::copyMeshProperties(Mesh* another_mesh)
	{
		mesh_lower_bound_ = another_mesh->mesh_lower_bound_;
		grid_spacing_ = another_mesh->grid_spacing_;
		cell_spacing_ = another_mesh->cell_spacing_;
		number_of_cells_ = another_mesh->number_of_cells_;
		number_of_grid_points_ = another_mesh->number_of_grid_points_;
	}
	//===================================================================//
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
	//===================================================================//
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
	//===================================================================//
	LevelSetData
		::LevelSetData(Real level_set, Vecd normal_direction)
		: phi_(level_set), n_(normal_direction), kappa_(0.0)
	{

	}
	//===================================================================//
	MeshBackground
		::MeshBackground(Vecd lower_bound, Vecd upper_bound, 
			Real grid_spacing, size_t buffer_size)
		: Mesh(lower_bound, upper_bound, grid_spacing, buffer_size)
	{
		number_of_grid_points_ = getNumberOfGridPoints(number_of_cells_);
	}
	//===================================================================//
	void BaseDataPackage
		::initializePackageGoemetry(Vecd& pkg_lower_bound, Real data_spacing)
	{
		mesh_lower_bound_ = pkg_lower_bound - Vecd(data_spacing * 0.5);;
		grid_spacing_ = data_spacing;
		data_lower_bound_ = pkg_lower_bound + Vecd(data_spacing * 0.5);
	}
	//===================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd lower_bound,
			Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: BaseMeshWithDataPackages(lower_bound, upper_bound,
			grid_spacing, buffer_size), sph_body_(NULL) {}
	//===================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: BaseMeshWithDataPackages(mesh_lower_bound, number_of_cells,
			cell_spacing), sph_body_(NULL) {}
	//===================================================================//
	LevelSet
		::LevelSet(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>(lower_bound,
			upper_bound, grid_spacing, buffer_size)
	{
		Real far_field_distance = cell_spacing_ * (Real)buffer_size;
		LevelSetDataPackage* positve_far_field = new LevelSetDataPackage();
		positve_far_field->initializeWithUniformData(far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(positve_far_field);
		LevelSetDataPackage* negatve_far_field = new LevelSetDataPackage();
		negatve_far_field->initializeWithUniformData(-far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(negatve_far_field);
	}
	//===================================================================//
	LevelSet
		::LevelSet(SPHBody* sph_body, Vecd lower_bound,
			Vecd upper_bound, Real grid_spacing, size_t buffer_size)
		: LevelSet(lower_bound,	upper_bound, grid_spacing, buffer_size)
	{
		sph_body_ = sph_body;
		InitializeDataPackages();
	}
	//===================================================================//
	LevelSet
		::LevelSet(Vecd mesh_lower_bound, 
			Vecu number_of_cells, Real cell_spacing)
		: MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>(mesh_lower_bound,
			number_of_cells, cell_spacing) {}
	//===================================================================//
	void LevelSet::InitializeDataPackages()
	{
		MeshFunctor initialize_data_in_a_cell = std::bind(&LevelSet::initializeDataInACell, this, _1, _2);
		MeshIterator(Vecu(0), number_of_cells_, initialize_data_in_a_cell);
		MeshFunctor tag_a_cell_inner_pkg = std::bind(&LevelSet::tagACellIsInnerPackage, this, _1, _2);
		MeshIterator(Vecu(0), number_of_cells_, tag_a_cell_inner_pkg);
		MeshFunctor initial_address_in_a_cell = std::bind(&LevelSet::initializeAdressesInACell, this, _1, _2);
		MeshIterator(Vecu(0), number_of_cells_, initial_address_in_a_cell);
	}
	//=================================================================================================//
	void LevelSet::initializeAdressesInACell(Vecu cell_index, Real dt)
	{
		initializeOneVariableAdressesInACell<Real, &LevelSetDataPackage::phi_addrs_, &LevelSetDataPackage::phi_>(cell_index);
		initializeOneVariableAdressesInACell<Vecd, &LevelSetDataPackage::n_addrs_, &LevelSetDataPackage::n_>(cell_index);
	}
	//=================================================================================================//
	Vecd LevelSet::probeNormalDirection(Vecd position)
	{
		return probeMesh<Vecd, &LevelSetDataPackage::n_addrs_, &LevelSetDataPackage::n_>(position);
	}
	//=================================================================================================//
	Real LevelSet::probeLevelSet(Vecd position)
	{
		return probeMesh<Real, &LevelSetDataPackage::phi_addrs_, &LevelSetDataPackage::phi_>(position);
	}
	//===================================================================//
	MultiresolutionLevelSet
		::MultiresolutionLevelSet(SPHBody* sph_body, Vecd lower_bound, Vecd upper_bound,
		Real reference_cell_spacing, size_t total_levels, size_t buffer_size)
		: MultilevelMesh<BaseLevelSet, LevelSet>(lower_bound, upper_bound,
			reference_cell_spacing, total_levels, buffer_size)
	{
		sph_body_ = sph_body;
		for (size_t level = 0; level != total_levels_; ++level) {
			mesh_levels_[level]->setSPHBody(sph_body);
			mesh_levels_[level]->InitializeDataPackages();
		}
	}
	//===================================================================//

}