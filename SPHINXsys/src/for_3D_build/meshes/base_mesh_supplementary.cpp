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
	Vecu BaseMesh::transfer1DtoMeshIndex(Vecu mesh_size, size_t i)
	{
		size_t row_times_column_size = mesh_size[1] * mesh_size[2];
		size_t page = i / row_times_column_size;
		size_t left_over = (i - page * row_times_column_size);
		size_t row_size = mesh_size[2];
		size_t column = left_over / row_size;
		return Vecu(page, column, left_over - column * row_size);
	}
	//=================================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index)
	{
		return mesh_index[0] * mesh_size[1] * mesh_size[2] 
			 + mesh_index[1] * mesh_size[2] 
			 + mesh_index[2];
	}
	//=================================================================================================//
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (size_t i = 0; i != pkg_size_; ++i)
			for (size_t j = 0; j != pkg_size_; ++j)
				for (size_t k = 0; k != pkg_size_; ++k)
				{
					pkg_data_[i][j][k].phi_ = level_set;
					pkg_data_[i][j][k].n_ = normal_direction;
				}
		for (size_t l = 0; l != number_of_grid_points_[0]; ++l)
			for (size_t m = 0; m != number_of_grid_points_[1]; ++m)
				for (size_t n = 0; n != number_of_grid_points_[2]; ++n) {
					pkg_data_addrs_[l][m][n] = &(pkg_data_[0][0][0]);
				}

	}
	//=================================================================================================//
	void  LevelSetDataPackage
		::initializeDataPackage(SPHBody* sph_body)
	{
		for (size_t i = 0; i != pkg_size_; ++i)
			for (size_t j = 0; j != pkg_size_; ++j)
				for (size_t k = 0; k != pkg_size_; ++k)
				{
					Vec3d position = data_lower_bound_ + Vec3d((Real)i * grid_spacing_, (Real)j * grid_spacing_, (Real)k * grid_spacing_);
					Real phi_from_surface = (sph_body->ClosestPointOnBodySurface(position) - position).norm();
					pkg_data_[i][j][k].phi_
						= sph_body->checkBodyShapeContain(position) ? - phi_from_surface : phi_from_surface;
				}
	}
	//=================================================================================================//
	void LevelSet::initializeDataInACell(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

		Vecd cell_position = CellPositionFromIndexes(cell_index);
		Vecd closet_pnt_on_face = sph_body_->ClosestPointOnBodySurface(cell_position);
		Real measure = getMinAbsoluteElement(closet_pnt_on_face - cell_position);
		if (measure < cell_spacing_) {
			LevelSetDataPackage* new_data_pkg = data_pkg_pool_.malloc();
			Vecd pkg_lower_bound = GridPositionFromCellPosition(cell_position);
			new_data_pkg->initializePackageGeometry(pkg_lower_bound, data_spacing_);
			new_data_pkg->initializeDataPackage(sph_body_);
			core_data_pkgs_.push_back(new_data_pkg);
			new_data_pkg->is_core_pkg_ = true;
			data_pkg_addrs_[i][j][k] = new_data_pkg;
		}
		else {
			data_pkg_addrs_[i][j][k] = sph_body_->checkBodyShapeContain(cell_position) ?
				singular_data_pkgs_addrs[0] : singular_data_pkgs_addrs[1];
		}
	}	
	//=================================================================================================//
	void LevelSet::tagACellIsInnerPackage(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

		bool is_inner_pkg = false;

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++n)
					if (data_pkg_addrs_[l][m][n]->is_core_pkg_) is_inner_pkg = true;
		
		if (is_inner_pkg) 
		{
			LevelSetDataPackage* current_data_pkg = data_pkg_addrs_[i][j][k];
			if (current_data_pkg->is_core_pkg_) 
			{
				current_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(current_data_pkg);
			}
			else {
				LevelSetDataPackage* new_data_pkg = data_pkg_pool_.malloc();
				Vecd cell_position = CellPositionFromIndexes(cell_index);
				Vecd pkg_lower_bound = GridPositionFromCellPosition(cell_position);
				new_data_pkg->initializePackageGeometry(pkg_lower_bound, data_spacing_);
				new_data_pkg->initializeDataPackage(sph_body_);
				new_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(new_data_pkg);
				data_pkg_addrs_[i][j][k] = new_data_pkg;
			}
		}
	}
	//=================================================================================================//
	void LevelSet::writeMeshToPltFile(ofstream& output_file)
	{
		Vecu number_of_operation = total_number_of_data_points_;

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "z, " << "phi, " << "n_x, "<< "n_y, " << "n_z " << "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = DataPositionFromGlobalIndex(Vecu(i, j, k));
					output_file << data_position[0] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = DataPositionFromGlobalIndex(Vecu(i, j, k));
					output_file << data_position[1] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = DataPositionFromGlobalIndex(Vecu(i, j, k));
					output_file << data_position[2] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<LevelSetData>(Vecu(i, j, k)).phi_ << " ";

				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<LevelSetData>(Vecu(i, j, k)).n_[0] << " ";

				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<LevelSetData>(Vecu(i, j, k)).n_[1] << " ";

				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<LevelSetData>(Vecu(i, j, k)).n_[2] << " ";

				}
				output_file << " \n";
			}
	}
	//=================================================================================================//
}
