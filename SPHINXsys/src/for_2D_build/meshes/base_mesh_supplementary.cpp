/**
 * @file 	base_mesh_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
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

namespace SPH {
	//===========================================================//
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt)
	{
		for (size_t i = index_begin[0]; i != index_end[0]; ++i)
			for (size_t j = index_begin[1]; j != index_end[1]; ++j) {
				mesh_functor(Vecu(i, j), dt);
			}
	}
	//===================================================================//
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
	//===========================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(Vecu mesh_size, size_t i)
	{
		int row_size = mesh_size[1];
		int column = i / row_size;
		return Vecu(column, i - column * row_size);
	}
	//===================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index)
	{
		return mesh_index[0] * mesh_size[1] + mesh_index[1];
	}
	//===================================================================//
	void MeshBackground
		::AllocateMeshDataMatrix()
	{
		Allocate2dArray(mesh_background_data_, number_of_grid_points_);
	}
	//===================================================================//
	void MeshBackground
		::DeleteMeshDataMatrix()
	{
		Delete2dArray(mesh_background_data_, number_of_grid_points_);
	}
	//===================================================================//
	void MeshBackground::InitializeLevelSetData(SPHBody &body)
	{
		//intialise the corresponding level set .
		Vecu number_of_operation = number_of_grid_points_;
		parallel_for(blocked_range2d<size_t>
			(0, number_of_operation[0], 0, number_of_operation[1]),
			[&](const blocked_range2d<size_t>& r) {
			for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
				for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
				{
					Vec2d grid_position = GridPositionFromIndexes(Vecu(i, j));
					Vec2d closet_pnt_on_face(0, 0);
					Real phi_from_surface = 0.0;

					body.ClosestPointOnBodySurface(grid_position, closet_pnt_on_face, phi_from_surface);
					mesh_background_data_[i][j].phi_ = phi_from_surface;
					mesh_background_data_[i][j].n_ = closet_pnt_on_face - grid_position;
				}
		}, ap);
	}
	//===================================================================//
	void MeshBackground::ComputeCurvatureFromLevelSet(SPHBody &body)
	{
		Vecu number_of_operation = number_of_grid_points_;

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				mesh_background_data_[i][j].kappa_ = 0.0;
			}
		}

		parallel_for(blocked_range2d<size_t>
			(1, number_of_operation[0]-1, 1, number_of_operation[1]-1),
			[&](const blocked_range2d<size_t>& r) 
			{
			for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
				for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
				{
					Real grad_x  = 0.5 * (mesh_background_data_[i+1][j].phi_ - mesh_background_data_[i-1][j].phi_) / grid_spacing_;
					Real grad_y  = 0.5 * (mesh_background_data_[i][j+1].phi_ - mesh_background_data_[i][j-1].phi_) / grid_spacing_;

					Real grad_xy =0.25 * (mesh_background_data_[i+1][j+1].phi_ - mesh_background_data_[i-1][j+1].phi_
					 					- mesh_background_data_[i+1][j-1].phi_ + mesh_background_data_[i-1][j-1].phi_) 
											/ grid_spacing_ / grid_spacing_;
					Real grad_xx = 	(mesh_background_data_[i+1][j].phi_
									 - 2.0 * mesh_background_data_[i][j].phi_
								 	+ mesh_background_data_[i-1][j].phi_)/ grid_spacing_ / grid_spacing_;
					Real grad_yy = (mesh_background_data_[i][j+1].phi_ - 
									 2.0 * mesh_background_data_[i][j].phi_ + 
										mesh_background_data_[i][j-1].phi_ )/ grid_spacing_ / grid_spacing_;
					
					Real grad_phi = grad_x * grad_x + grad_y * grad_y;
					mesh_background_data_[i][j].kappa_ = (grad_xx * grad_y * grad_y - 2.0 * grad_x * grad_y * grad_xy + 
											grad_yy * grad_x * grad_x) / (grad_phi * sqrt(grad_phi) + 1.0e-15);
				}
		}, ap);
	}
	//===================================================================//
	Vecd MeshBackground::ProbeNormalDirection(Vecd Point)
	{
		Vec2u grid_idx = GridIndexesFromPosition(Point);
		Vec2d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec2d dis_grid = (Point - grid_pos) / grid_spacing_;

		Vec2d norm_to_face
			= mesh_background_data_[grid_idx[0]][grid_idx[1]].n_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]].n_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1].n_* (1.0 - dis_grid[0]) * dis_grid[1] 
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1] + 1].n_* dis_grid[0] * dis_grid[1];

		return  norm_to_face;
	}
	//===================================================================//
	Real MeshBackground::ProbeLevelSet(Vecd Point)
	{
		Vec2u grid_idx = GridIndexesFromPosition(Point);
		Vec2d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec2d dis_grid = (Point - grid_pos) / grid_spacing_;

		Real phi 
			= mesh_background_data_[grid_idx[0]][grid_idx[1]].phi_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]].phi_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1].phi_* (1.0 - dis_grid[0]) * dis_grid[1]
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1] + 1].phi_* dis_grid[0] * dis_grid[1];

		return  phi;
	}
	//===================================================================//
	Real MeshBackground::ProbeCurvature(Vecd Point)
	{
		Vec2u grid_idx = GridIndexesFromPosition(Point);
		Vec2d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec2d dis_grid = (Point - grid_pos) / grid_spacing_;

		Real kappa 
			= mesh_background_data_[grid_idx[0]][grid_idx[1]].kappa_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]].kappa_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1].kappa_* (1.0 - dis_grid[0]) * dis_grid[1]
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1] + 1].kappa_* dis_grid[0] * dis_grid[1];

		return  kappa;
	}
	//===================================================================//
	void MeshBackground::WriteMeshToVtuFile(ofstream &output_file)
	{
		cout << "\n This function WriteMeshToVtuFile is not done. Exit the program! \n";
		exit(0);

	}
	//===================================================================//
	void MeshBackground::WriteMeshToPltFile(ofstream &output_file)
	{
		Vecu number_of_operation = number_of_grid_points_;

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "phi, " << "n_x, " << "n_y, " << "kappa, " << "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd cell_position = GridPositionFromIndexes(Vecu(i, j));
				output_file << cell_position[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd cell_position = GridPositionFromIndexes(Vecu(i, j));
				output_file << cell_position[1] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << mesh_background_data_[i][j].phi_ << " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << mesh_background_data_[i][j].n_[0] << " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << mesh_background_data_[i][j].n_[1] << " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << mesh_background_data_[i][j].kappa_ << " ";

			}
			output_file << " \n";
		}

	}
	//===================================================================//
	void LevelSetDataPackage::AllocateMeshDataMatrix()
	{
		Allocate2dArray(phi_, number_of_grid_points_);
		Allocate2dArray(n_, number_of_grid_points_);
		Allocate2dArray(phi_addrs_, number_of_addrs_);
		Allocate2dArray(n_addrs_, number_of_addrs_);
	}
	//===================================================================//
	void LevelSetDataPackage::DeleteMeshDataMatrix()
	{
		Delete2dArray(phi_, number_of_grid_points_);
		Delete2dArray(n_, number_of_grid_points_);
		Delete2dArray(phi_addrs_, number_of_addrs_);
		Delete2dArray(n_addrs_, number_of_addrs_);
	}
	//===================================================================//
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (size_t i = 0; i != pkg_size_; ++i)
			for (size_t j = 0; j != pkg_size_; ++j) {
				phi_[i][j] = level_set;
				n_[i][j] = normal_direction;
			}
	}
	//=================================================================================================//
	void  LevelSetDataPackage
		::initializeDataPackage(SPHBody* sph_body)
	{
		for (size_t i = 0; i != pkg_size_; ++i)
			for (size_t j = 0; j != pkg_size_; ++j) {
				Vecd position = data_lower_bound_ + Vecd((Real)i * grid_spacing_, (Real)j * grid_spacing_);
				Vecd closet_pnt_on_face(0);
				sph_body->ClosestPointOnBodySurface(position, closet_pnt_on_face, phi_[i][j]);
				n_[i][j] = closet_pnt_on_face - position;
			}
	}
	//===================================================================//
	void LevelSet::initializeDataInACell(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		Vecd cell_position = CellPositionFromIndexes(cell_index);
		Vecd closet_pnt_on_face(0, 0);
		Real phi_from_surface = 0.0;
		sph_body_->ClosestPointOnBodySurface(cell_position, closet_pnt_on_face, phi_from_surface);
		Real measure = getMinAbslouteElement(closet_pnt_on_face - cell_position);
		if (measure < cell_spacing_) {
			LevelSetDataPackage* new_data_pkg = data_pakg_pool_.malloc();
			Vecd pkg_lower_bound = getGridPositionFromCellPosition(cell_position);
			new_data_pkg->initializePackageGoemetry(pkg_lower_bound, data_spacing_);
			new_data_pkg->initializeDataPackage(sph_body_);
			core_data_pkgs_.push_back(new_data_pkg);
			new_data_pkg->is_core_pkg_ = true;
			data_pkg_addrs_[i][j] = new_data_pkg;
		}
		else {
			data_pkg_addrs_[i][j] = phi_from_surface > 0 ?
			 singular_data_pkgs_addrs[0] : singular_data_pkgs_addrs[1];
		}
	}
	//===================================================================//
	void LevelSet::tagACellIsInnerPackage(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		bool is_inner_pkg = false;
		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				if (data_pkg_addrs_[l][m]->is_core_pkg_) is_inner_pkg = true;

		if (is_inner_pkg) {
			LevelSetDataPackage* current_data_pkg = data_pkg_addrs_[i][j];
			if (current_data_pkg->is_core_pkg_) {
				current_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(current_data_pkg);
			}
			else {
				LevelSetDataPackage* new_data_pkg = data_pakg_pool_.malloc();
				Vecd cell_position = CellPositionFromIndexes(cell_index);
				Vecd pkg_lower_bound = getGridPositionFromCellPosition(cell_position);
				new_data_pkg->initializePackageGoemetry(pkg_lower_bound, data_spacing_);
				new_data_pkg->initializeDataPackage(sph_body_);
				new_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(new_data_pkg);
				data_pkg_addrs_[i][j] = new_data_pkg;
			}
		}
	}
	//===================================================================//
	void LevelSet::WriteMeshToPltFile(ofstream& output_file)
	{
		Vecu number_of_operation = number_of_data_;

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "phi, " << "n_x, " << "n_y " << "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = getDataPositionFromIndex(Vecu(i, j));
				output_file << data_position[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = getDataPositionFromIndex(Vecu(i, j));
				output_file << data_position[1] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << getValueFromGlobalDataIndex<Real, &LevelSetDataPackage::phi_>(Vecu(i, j)) << " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << getValueFromGlobalDataIndex<Vecd, &LevelSetDataPackage::n_>(Vecu(i, j))[0]<< " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << getValueFromGlobalDataIndex<Vecd, &LevelSetDataPackage::n_>(Vecu(i, j))[1] << " ";

			}
			output_file << " \n";
		}
	}
	//===================================================================//
}
