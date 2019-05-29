#include "base_mesh.h"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "neighboring_particle.h"
#include "base_data_package.h"

#include "math.h"

namespace SPH {

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
		//intialise the corresponding level set .
		Vecu number_of_operation = number_of_grid_points_;
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
}
