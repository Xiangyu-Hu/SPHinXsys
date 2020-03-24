#include "base_mesh.h"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "neighboring_particle.h"
#include "base_data_package.h"

namespace SPH {
	//===========================================================//
	Vecu Mesh::transfer1DtoMeshIndex(Vecu mesh_size, size_t i)
	{
		size_t row_times_column_size = mesh_size[1] * mesh_size[2];
		size_t page = i / row_times_column_size;
		size_t left_over = (i - page * row_times_column_size);
		size_t row_size = mesh_size[2];
		size_t column = left_over / row_size;
		return Vecu(page, column, left_over - column * row_size);
	}
	//===================================================================//
	size_t Mesh::transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index)
	{
		return mesh_index[0] * mesh_size[1] * mesh_size[2] 
			 + mesh_index[1] * mesh_size[2] 
			 + mesh_index[2];
	}
	//===================================================================//
	void MeshBackground
		::AllocateMeshDataMatrix()
	{
		Allocate3dArray(mesh_background_data_, number_of_grid_points_);
	}
	//===================================================================//
	void MeshBackground
		::DeleteMeshDataMatrix()
	{
		Delete3dArray(mesh_background_data_, number_of_grid_points_);
	}
	//===================================================================//
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
	//===================================================================//
	Vecd MeshBackground::ProbeNormalDirection(Vecd Point)
	{
		Vec3u grid_idx = GridIndexesFromPosition(Point);
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
	//===================================================================//
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
	//===================================================================//
	void MeshBackground::ComputeCurvatureFromLevelSet(SPHBody &body)
	{

	}
	//===================================================================//
	Real MeshBackground::ProbeCurvature(Vecd Point)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);

		return  0.0;
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
}