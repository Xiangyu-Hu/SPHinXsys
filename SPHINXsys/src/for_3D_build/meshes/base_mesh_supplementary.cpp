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
	void MeshBackground
		::AllocateMeshDataMatrix()
	{
		Allocate3dArray(mesh_background_data_, number_of_grid_points_);
	}
	//=================================================================================================//
	void MeshBackground
		::DeleteMeshDataMatrix()
	{
		Delete3dArray(mesh_background_data_, number_of_grid_points_);
	}
	//=================================================================================================//
	void MeshBackground::InitializeLevelSetData(SPHBody &body)
	{
		Vecu number_of_operation = number_of_grid_points_;
		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					mesh_background_data_[i][j][k].phi_ = 0.0;
					mesh_background_data_[i][j][k].n_ = Vecd(0.0);
					mesh_background_data_[i][j][k].kappa_ = 0.0;
				}
			}
		}
		/** Compute the phi and norm form mesh. */
		parallel_for(blocked_range3d<size_t>
			(0, number_of_operation[0], 0, number_of_operation[1], 0, number_of_operation[2]),
			[&](const blocked_range3d<size_t>& r) {
			for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
				for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
					{
						Vecd grid_position = GridPositionFromIndexes(Vecu(i, j, k));
						Vecd closet_pnt_on_face(0, 0, 0);
						Real phi_from_surface = 0.0;

						body.ClosestPointOnBodySurface(grid_position, closet_pnt_on_face, phi_from_surface);

						mesh_background_data_[i][j][k].phi_ = phi_from_surface;
						mesh_background_data_[i][j][k].n_ = closet_pnt_on_face - grid_position;
					}
	 }, ap);
	}
	//=================================================================================================//
	Vecd MeshBackground::ProbeNormalDirection(Vecd Point)
	{
		Vec3u grid_idx = CellIndexesFromPosition(Point);
		Vec3d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec3d dis_grid = (Point - grid_pos) / grid_spacing_;

		Vec3d bilinear_1
			= mesh_background_data_[grid_idx[0]][grid_idx[1]][grid_idx[2]].n_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]].n_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]].n_* (1.0 - dis_grid[0]) * dis_grid[1] 
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]+1][grid_idx[2]].n_* dis_grid[0] * dis_grid[1];
		Vec3d bilinear_2
			= mesh_background_data_[grid_idx[0]][grid_idx[1]][grid_idx[2]+1].n_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]+1].n_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]+1].n_* (1.0 - dis_grid[0]) * dis_grid[1] 
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]+1][grid_idx[2]+1].n_* dis_grid[0] * dis_grid[1];
		return  bilinear_1 * (1.0 - dis_grid[2]) + bilinear_2 * dis_grid[2];
	}
	//=================================================================================================//
	Real MeshBackground::ProbeLevelSet(Vecd Point)
	{
		Vec3u grid_idx = GridIndexesFromPosition(Point);
		Vec3d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec3d dis_grid = (Point - grid_pos) / grid_spacing_;


		Real bilinear_1
			= mesh_background_data_[grid_idx[0]][grid_idx[1]][grid_idx[2]].phi_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]].phi_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]].phi_* (1.0 - dis_grid[0]) * dis_grid[1] 
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]+1][grid_idx[2]].phi_* dis_grid[0] * dis_grid[1];
		Real bilinear_2
			= mesh_background_data_[grid_idx[0]][grid_idx[1]][grid_idx[2]+1].phi_ * (1.0 - dis_grid[0]) * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]+1].phi_* dis_grid[0] * (1.0 - dis_grid[1])
			+ mesh_background_data_[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]+1].phi_* (1.0 - dis_grid[0]) * dis_grid[1] 
			+ mesh_background_data_[grid_idx[0] + 1][grid_idx[1]+1][grid_idx[2]+1].phi_* dis_grid[0] * dis_grid[1];
		return  bilinear_1 * (1.0 - dis_grid[2]) + bilinear_2 * dis_grid[2];
	}
	//=================================================================================================//
	void MeshBackground::ComputeCurvatureFromLevelSet(SPHBody &body)
	{

	}
	//=================================================================================================//
	Real MeshBackground::ProbeCurvature(Vecd Point)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);

		return  0.0;
	}
	//=================================================================================================//
	void MeshBackground::WriteMeshToVtuFile(ofstream &output_file)
	{
		cout << "\n This function WriteMeshToVtuFile is not done. Exit the program! \n";
		exit(0);

	}
	//=================================================================================================//
	void MeshBackground::WriteMeshToPltFile(ofstream &output_file)
	{
		Vecu number_of_operation = number_of_grid_points_;

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "z, " << "phi, " << "n_x, " << "n_y, " << "n_z, "<< "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[0] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[1] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[2] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].phi_ << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[0] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[1] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[2] << " ";
				}
				output_file << " \n";
			}
		}
	}
	//=================================================================================================//
	void LevelSetDataPackage::AllocateMeshDataMatrix()
	{
		Allocate3dArray(phi_, number_of_grid_points_);
		Allocate3dArray(n_, number_of_grid_points_);
		Allocate3dArray(phi_addrs_, number_of_addrs_);
		Allocate3dArray(n_addrs_, number_of_addrs_);
	}
	//=================================================================================================//
	void LevelSetDataPackage::DeleteMeshDataMatrix()
	{
		Delete3dArray(phi_, number_of_grid_points_);
		Delete3dArray(n_, number_of_grid_points_);
		Delete3dArray(phi_addrs_, number_of_addrs_);
		Delete3dArray(n_addrs_, number_of_addrs_);
	}
	//=================================================================================================//
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (size_t i = 0; i != 4; ++i)
			for (size_t j = 0; j != 4; ++j)
				for (size_t k = 0; k != 4; ++k)
				{
					phi_[i][j][k] = level_set;
					n_[i][j][k] = normal_direction;
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
					Vecd position = mesh_lower_bound_ + Vecd((Real)i * grid_spacing_, (Real)j * grid_spacing_);
					Vecd closet_pnt_on_face(0);
					sph_body->ClosestPointOnBodySurface(position, closet_pnt_on_face, phi_[i][j][k]);
					n_[i][j][k] = closet_pnt_on_face - position;
				}
	}
	//=================================================================================================//
	void LevelSet::initializeDataInACell(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

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
			data_pkg_addrs_[i][j][k] = new_data_pkg;
		}
		else {
			data_pkg_addrs_[i][j][k] = phi_from_surface > 0 ?
				singular_data_pkgs_addrs[0] : singular_data_pkgs_addrs[1];
		}
	}	
	//=================================================================================================//
	void LevelSet::tagACellIsInnerPackage(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++n)
					if (data_pkg_addrs_[l][m][n]->is_core_pkg_) goto endloop;
	endloop:
		{
			LevelSetDataPackage* current_data_pkg = data_pkg_addrs_[i][j][k];
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
				data_pkg_addrs_[i][j][k] = new_data_pkg;
			}
		}
	}
	//=================================================================================================//
	void LevelSet::WriteMeshToPltFile(ofstream& output_file)
	{
		cout << "\n This function LevelSet::WriteMeshToPltFile is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
}
