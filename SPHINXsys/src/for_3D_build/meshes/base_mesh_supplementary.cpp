#include "base_mesh.h"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "neighboring_particle.h"
#include "base_data_package.h"

namespace SPH {

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
		//intialise the corresponding level set .
		Vecu number_of_operation = number_of_grid_points_;
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
						if (abs(phi_from_surface) > 2.0 * grid_spacing_) {
							mesh_background_data_[i][j][k].n_ = (0.0, 0.0, 0.0);
						}
						else {
							mesh_background_data_[i][j][k].n_ = closet_pnt_on_face - grid_position;
						}
					}
		}, ap);
	}
	//===================================================================//
	Vecd MeshBackground::ProbeNormalDirection(Vecd Point)
	{
		Vecd normal_direction(0);

		return  normal_direction;
	}
	//===================================================================//
	Real MeshBackground::ProbeLevelSet(Vecd Point)
	{
		Real phi = 0.0;

		return  phi;
	}
	//===================================================================//
	void MeshBackground::ComputeCurvatureFromLevelSet(SPHBody &body)
	{

		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);
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
		output_file << "variables= " << "x, " << "y, " << "z " << "phi, " << "n_x, " << "n_y, " << "n_z" << "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[0]
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[0] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[1] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					Vecd cell_position = GridPositionFromIndexes(Vecu(i, j, k));
					output_file << cell_position[2] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].phi_ << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[0] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[1] << " ";
				}
				output_file << " \n";
			}
		}

		for (size_t k = 0; k != number_of_operation[0]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[2]; ++i)
				{
					output_file << mesh_background_data_[i][j][k].n_[2] << " ";
				}
				output_file << " \n";
			}
		}
	}
}