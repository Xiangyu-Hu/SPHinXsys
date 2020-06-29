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
						}
			}, ap);

		/** Compute the normal direction mesh. */
		parallel_for(blocked_range3d<size_t>
			(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1, 1, number_of_operation[2] - 1),
			[&](const blocked_range3d<size_t>& r) {
				for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
						for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
						{
							Real dphidx = mesh_background_data_[i + 1][j][k].phi_ - mesh_background_data_[i - 1][j][k].phi_;
							Real dphidy = mesh_background_data_[i][j + 1][k].phi_ - mesh_background_data_[i][j - 1][k].phi_;
							Real dphidz = mesh_background_data_[i][j][k + 1].phi_ - mesh_background_data_[i][j][k - 1].phi_;
							Vecd normal = Vecd(dphidx, dphidy, dphidz);
							mesh_background_data_[i][j][k].n_ = normal / (normal.norm() + 1.0e-16);
						}
			}, ap);
	}

	//===================================================================//
	void MeshBackground::ClearLevelSetData(SPHBody &body)
	{
		Real epsilon_ = 0.75*grid_spacing_;
		std::vector<Vecu> position_S_ = GetZeroLevelSetCutCell();
		std::vector<Vecu> position_P_ = GetPositiveCutCell();
		std::vector<Vecu> position_N_ = GetNegativeCutCell();
		Vecu number_of_operation = number_of_grid_points_;

		for (size_t x = 0; x < position_S_.size(); x++)
		{
			if (IfNeighborCellsInAuxiliayBandP(position_S_[x]) == 0) //function equal to zero means no neighborcells in band P 
			{
				std::vector<Real>d_P_;
				for (int i = 0; i < 9; i++)
				{
					for (int j = 0; j < 9; j++) 
					{
						for (int k = 0; k < 9; k++)
						{
							if (mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].gamma_P_)
							{
								Real phi_P_ = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].phi_;
								Vecd norm_to_face = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].n_
									/ mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].n_.norm();
								Vec3d grid_pos_S = GridPositionFromIndexes(position_S_[x]);
								Vec3d grid_pos_P = GridPositionFromIndexes(Vecu((position_S_[x][0] - 4 + i), (position_S_[x][1] - 4 + j), (position_S_[x][2] - 4 + k)));

								d_P_.push_back(sqrt(powern((grid_pos_S[0] - grid_pos_P[0] + (phi_P_ - epsilon_)*norm_to_face[0]), 2) +
									powern((grid_pos_S[1] - grid_pos_P[1] + (phi_P_ - epsilon_)*norm_to_face[1]), 2) + 
									powern((grid_pos_S[2] - grid_pos_P[2] + (phi_P_ - epsilon_)*norm_to_face[2]), 2)));
							}
						}
						
					}
				}
				Real d_min_P_;
				if (d_P_.size() != 0) {
					d_min_P_ = *min_element(d_P_.begin(), d_P_.end());
				}
				else {
					d_min_P_ = 2.0 * grid_spacing_;
				}
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]][position_S_[x][2]].phi_ = -(d_min_P_ - epsilon_);
				Real phi_temp_fortest = -(d_min_P_ - epsilon_);
				int a = position_S_[x][0];
				int b = position_S_[x][1];
				int c = position_S_[x][2];
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]][position_S_[x][2]].gamma_S_ = false;
			}

			if (IfNeighborCellsInAuxiliayBandN(position_S_[x]) == 0) //function equal to zero means no neighborcells in band N 
			{
				
				std::vector<Real>d_N_;
				for (int i = 0; i < 9; i++)
				{
					for (int j = 0; j < 9; j++) 
					{
						for (int k = 0; k < 9; k++)
						{
							if (mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].gamma_N_)
							{
								Real phi_N_ = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].phi_;
								Vecd norm_to_face = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].n_
									/ mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j][position_S_[x][2] - 4 + k].n_.norm();
								Vec3d grid_pos_S = GridPositionFromIndexes(position_S_[x]);
								Vec3d grid_pos_N = GridPositionFromIndexes(Vecu((position_S_[x][0] - 4 + i), (position_S_[x][1] - 4 + j), (position_S_[x][2] - 4 + k)));

								d_N_.push_back(sqrt(powern((grid_pos_S[0] - grid_pos_N[0] + (phi_N_ + epsilon_)*norm_to_face[0]), 2) +
									powern((grid_pos_S[1] - grid_pos_N[1] + (phi_N_ + epsilon_)*norm_to_face[1]), 2) +
									powern((grid_pos_S[2] - grid_pos_N[2] + (phi_N_ + epsilon_)*norm_to_face[2]), 2)));
							}
						}
					}
				}
				Real d_min_N_;
				if (d_N_.size() != 0) {
					d_min_N_ = *min_element(d_N_.begin(), d_N_.end());
				}
				else {
					d_min_N_ = 3.0 * grid_spacing_;
				}
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]][position_S_[x][2]].phi_ = d_min_N_ + epsilon_;
				Real phi_temp_fortest = d_min_N_ + epsilon_;
				int a = position_S_[x][0];
				int b = position_S_[x][1];
				int c = position_S_[x][2];
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]][position_S_[x][2]].gamma_S_ = false;
			}
		}
	}


	//===================================================================//
	void MeshBackground::MarkZeroLevelSetCutCell()
	{
		Vecu number_of_operation = number_of_grid_points_;


		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				for (int k = 1; k < number_of_operation[2] - 1; k++)
				{
					Real phi_1_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j][k + 1].phi_ + mesh_background_data_[i][j][k + 1].phi_ +
						mesh_background_data_[i][j + 1][k + 1].phi_ + mesh_background_data_[i + 1][j + 1][k + 1].phi_);
					Real phi_2_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_ +
						mesh_background_data_[i + 1][j + 1][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_);
					Real phi_3_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_);
					Real phi_4_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i + 1][j][k + 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_);
					Real phi_5_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k + 1].phi_ + mesh_background_data_[i][j + 1][k + 1].phi_);
					Real phi_6_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_);
					Real phi_7_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_);
					Real phi_8_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_);

					if (phi_1_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_2_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_3_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_4_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_5_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_6_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_7_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_8_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else {
						mesh_background_data_[i][j][k].gamma_S_ = false;
					}
				}
			}
		}
	}

	//===================================================================//
	vector<Vecu> MeshBackground::GetZeroLevelSetCutCell()
	{
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_S_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				for (int k = 1; k < number_of_operation[2] - 1; k++)
				{
					Real phi_1_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j][k + 1].phi_ + mesh_background_data_[i][j][k + 1].phi_ +
						mesh_background_data_[i][j + 1][k + 1].phi_ + mesh_background_data_[i + 1][j + 1][k + 1].phi_);
					Real phi_2_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_ +
						mesh_background_data_[i + 1][j + 1][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_);
					Real phi_3_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_);
					Real phi_4_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i + 1][j][k + 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_);
					Real phi_5_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k + 1].phi_ + mesh_background_data_[i][j + 1][k + 1].phi_);
					Real phi_6_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_);
					Real phi_7_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_);
					Real phi_8_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_);

					if (phi_1_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_2_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_3_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_4_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_5_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_6_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_7_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
					else if (phi_8_ * mesh_background_data_[i][j][k].phi_ <= 0) {
						position_S_.push_back(Vecu(i, j, k));
						mesh_background_data_[i][j][k].gamma_S_ = true;
					}
				}
			}
		}
		return position_S_;
	}

	//===================================================================//
	vector<Vecu> MeshBackground::GetPositiveCutCell()
	{
		Real epsilon_ = 0.75*grid_spacing_;
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_P_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				for (int k = 1; k < number_of_operation[2] - 1; k++)
				{
					Real phi_1_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j][k + 1].phi_ + mesh_background_data_[i][j][k + 1].phi_ +
						mesh_background_data_[i][j + 1][k + 1].phi_ + mesh_background_data_[i + 1][j + 1][k + 1].phi_);
					Real phi_2_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_ +
						mesh_background_data_[i + 1][j + 1][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_);
					Real phi_3_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_);
					Real phi_4_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i + 1][j][k + 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_);
					Real phi_5_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k + 1].phi_ + mesh_background_data_[i][j + 1][k + 1].phi_);
					Real phi_6_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_);
					Real phi_7_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_);
					Real phi_8_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_);

					if ((phi_1_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_2_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_3_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_4_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_5_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_6_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_7_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
					else if ((phi_8_ - epsilon_) * (mesh_background_data_[i][j][k].phi_ - epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_P_ = true;
						position_P_.push_back(Vecu(i, j, k));
					}
				}
			}
		}
		return position_P_;
	}

	//===================================================================//
	vector<Vecu> MeshBackground::GetNegativeCutCell()
	{
		Real epsilon_ = 0.75*grid_spacing_;
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_N_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				for (int k = 1; k < number_of_operation[2] - 1; k++)
				{
					Real phi_1_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j][k + 1].phi_ + mesh_background_data_[i][j][k + 1].phi_ +
						mesh_background_data_[i][j + 1][k + 1].phi_ + mesh_background_data_[i + 1][j + 1][k + 1].phi_);
					Real phi_2_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i + 1][j + 1][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_ +
						mesh_background_data_[i + 1][j + 1][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_);
					Real phi_3_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i + 1][j][k - 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_);
					Real phi_4_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i + 1][j][k].phi_ +
						mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i + 1][j][k + 1].phi_ +
						mesh_background_data_[i + 1][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_);
					Real phi_5_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k + 1].phi_ + mesh_background_data_[i][j + 1][k + 1].phi_);
					Real phi_6_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i - 1][j + 1][k].phi_ + mesh_background_data_[i][j + 1][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_ +
						mesh_background_data_[i - 1][j + 1][k - 1].phi_ + mesh_background_data_[i][j + 1][k - 1].phi_);
					Real phi_7_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k - 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_);
					Real phi_8_ = (1.0 / 8.0) * (mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j - 1][k].phi_ +
						mesh_background_data_[i - 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j][k].phi_ +
						mesh_background_data_[i][j][k + 1].phi_ + mesh_background_data_[i][j - 1][k + 1].phi_ +
						mesh_background_data_[i - 1][j - 1][k + 1].phi_ + mesh_background_data_[i - 1][j][k + 1].phi_);

					if ((phi_1_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_2_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_3_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_4_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_5_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_6_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_7_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
					else if ((phi_8_ + epsilon_) * (mesh_background_data_[i][j][k].phi_ + epsilon_) <= 0) {
						mesh_background_data_[i][j][k].gamma_N_ = true;
						position_N_.push_back(Vecu(i, j, k));
					}
				}
			}
		}
		return position_N_;
	}

	//===================================================================//
	int MeshBackground::IfNeighborCellsInAuxiliayBandP(Vecu position)
	{
		int p_ = 0;
		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int k = 0; k < 3; k++) 
				{
					if (mesh_background_data_[position[0] - 1 + i][position[1] - 1 + j][position[2] - 1 + k].gamma_P_)
					{
						p_ = p_ + 1;
					};
				}
				
			}
		}
		return p_;
	}
	//===================================================================//
	int MeshBackground::IfNeighborCellsInAuxiliayBandN(Vecu position)
	{
		int n_ = 0;

		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int k = 0; k < 3; k++)
				{
					if (mesh_background_data_[position[0] - 1 + i][position[1] - 1 + j][position[2] - 1 + k].gamma_N_)
					{
						n_ = n_ + 1;
					};
				}
			}
		}
		return n_;
	}

	//===================================================================//
	void MeshBackground::ReinitializeLevelSetData(SPHBody &body)
	{
		Real dtau = 0.5*grid_spacing_;
		Real s, ss;
		Real dv_xp, dv_x, dv_xn; //p for positive, n for negative
		Real dv_yp, dv_y, dv_yn;
		Real dv_zp, dv_z, dv_zn;
		size_t x = 0;
		Vecu number_of_operation = number_of_grid_points_;

		while (x <= 100)
		{
			for (size_t i = 1; i < number_of_operation[0] - 1; i++)
			{
				for (size_t j = 1; j < number_of_operation[1] - 1; j++)
				{
					for (size_t k = 1; k < number_of_operation[2] - 1; k++)
					{
						if (mesh_background_data_[i][j][k].gamma_S_ != true)
						{

							s = mesh_background_data_[i][j][k].phi_ /
								sqrt(mesh_background_data_[i][j][k].phi_*mesh_background_data_[i][j][k].phi_ + grid_spacing_ * grid_spacing_);
							dv_xp = (mesh_background_data_[i + 1][j][k].phi_ - mesh_background_data_[i][j][k].phi_);
							dv_xn = (mesh_background_data_[i][j][k].phi_ - mesh_background_data_[i - 1][j][k].phi_);
							dv_x = dv_xp;
							if (s * dv_xp >= 0.0 && s * dv_xn >= 0.0) dv_x = dv_xn;
							if (s * dv_xp <= 0.0 && s * dv_xn <= 0.0) dv_x = dv_xp;
							if (s * dv_xp > 0.0 && s * dv_xn < 0.0) dv_x = 0.0;
							if (s * dv_xp < 0.0 && s * dv_xn > 0.0)
							{
								ss = s * (fabs(dv_xp) - fabs(dv_xn)) / (dv_xp - dv_xn);
								if (ss > 0.0) dv_x = dv_xn;
							}

							dv_yp = (mesh_background_data_[i][j + 1][k].phi_ - mesh_background_data_[i][j][k].phi_);
							dv_yn = (mesh_background_data_[i][j][k].phi_ - mesh_background_data_[i][j - 1][k].phi_);
							dv_y = dv_yp;
							if (s * dv_yp >= 0.0 && s * dv_yn >= 0.0) dv_y = dv_yn;
							if (s * dv_yp <= 0.0 && s * dv_yn <= 0.0) dv_y = dv_yp;
							if (s * dv_yp > 0.0 && s * dv_yn < 0.0) dv_y = 0.0;
							if (s * dv_yp < 0.0 && s * dv_yn > 0.0)
							{
								ss = s * (fabs(dv_yp) - fabs(dv_yn)) / (dv_yp - dv_yn);
								if (ss > 0.0) dv_y = dv_yn;
							}

							dv_zp = (mesh_background_data_[i][j][k + 1].phi_ - mesh_background_data_[i][j][k].phi_);
							dv_zn = (mesh_background_data_[i][j][k].phi_ - mesh_background_data_[i][j][k - 1].phi_);
							dv_z = dv_zp;
							if (s * dv_zp >= 0.0 && s * dv_zn >= 0.0) dv_z = dv_zn;
							if (s * dv_zp <= 0.0 && s * dv_zn <= 0.0) dv_z = dv_zp;
							if (s * dv_zp > 0.0 && s * dv_zn < 0.0) dv_z = 0.0;
							if (s * dv_zp < 0.0 && s * dv_zn > 0.0)
							{
								ss = s * (fabs(dv_zp) - fabs(dv_zn)) / (dv_zp - dv_zn);
								if (ss > 0.0) dv_z = dv_zn;
							}


							//time step 
							mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][k].phi_ -
								dtau * s * (sqrt(dv_x * dv_x / grid_spacing_ / grid_spacing_ + dv_y * dv_y / grid_spacing_ / grid_spacing_ + +dv_z * dv_z / grid_spacing_ / grid_spacing_) - 1.0);
						}
					
					}
				}
			}

			//update bounday data
			for (size_t i = 0; i < 1; i++) 
			{
				for (size_t j = 0; j < number_of_operation[1]; j++)
				{
					for (size_t k = 0; k < number_of_operation[2]; k++) 
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[1][j][k].phi_;
					}					
				}
			}

			for (size_t i = number_of_operation[0] - 1; i < number_of_operation[0]; i++)
			{
				for (size_t j = 0; j < number_of_operation[1]; j++)
				{
					for (size_t k = 0; k < number_of_operation[2]; k++)
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[number_of_operation[0] - 2][j][k].phi_;
					}
				}
			}

			for (size_t j = 0; j < 1; j++) 
			{
				for (size_t k = 0; k < number_of_operation[2]; k++)
				{ 
					for (size_t i = 1; i < number_of_operation[0] - 1; i++)
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][1][k].phi_;
					}
				}
			}

			for (size_t j = number_of_operation[1] - 1; j < number_of_operation[1]; j++) 
			{
				for (size_t k = 0; k < number_of_operation[2]; k++)
				{

					for (size_t i = 1; i < number_of_operation[0] - 1; i++)
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][number_of_operation[1] - 2][k].phi_;
					}
				}
			}

			for (size_t k = 0; k < 1; k++)
			{
				for (size_t i = 1; i < number_of_operation[0] -1; i++)
				{
					for (size_t j = 1; j < number_of_operation[1] - 1; j++)
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][1].phi_;
					}
				}
			}

			for (size_t k = number_of_operation[2] - 1; k < number_of_operation[2]; k++)
			{
				for (size_t i = 0; i < number_of_operation[0] - 1; i++)
				{

					for (size_t j = 1; j < number_of_operation[1] - 1; j++)
					{
						mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][number_of_operation[2] - 2].phi_;
					}
				}
			}

			x = x + 1;
		}

		//update the normal direction
		for (size_t i = 2; i <= number_of_operation[0] - 2; i++)
		{
			for (size_t j = 2; j <= number_of_operation[1] - 2; j++)
			{
				for (size_t k = 2; k <= number_of_operation[2] - 2; k++)
				{
					mesh_background_data_[i][j][k].n_[0] = ((-1) * SGN(mesh_background_data_[i + 1][j][k].phi_) * mesh_background_data_[i + 1][j][k].phi_
						- (-1) * SGN(mesh_background_data_[i - 1][j][k].phi_) * mesh_background_data_[i - 1][j][k].phi_) / (2.0 * grid_spacing_);
					mesh_background_data_[i][j][k].n_[1] = ((-1) * SGN(mesh_background_data_[i][j + 1][k].phi_) * mesh_background_data_[i][j + 1][k].phi_
						- (-1) * SGN(mesh_background_data_[i][j - 1][k].phi_) * mesh_background_data_[i][j - 1][k].phi_) / (2.0 * grid_spacing_);
					mesh_background_data_[i][j][k].n_[2] = ((-1) * SGN(mesh_background_data_[i][j][k + 1].phi_) * mesh_background_data_[i][j][k + 1].phi_
						- (-1) * SGN(mesh_background_data_[i][j][k - 1].phi_) * mesh_background_data_[i][j][k - 1].phi_) / (2.0 * grid_spacing_);
				}
			}
		}
	}

	//===================================================================//
	void MeshBackground::SmoothLevelSetByCurvature(SPHBody &body, Real smoothing_coe)
	{
		Real lambda = 0.02;
		Real dtau = 0.02*grid_spacing_;
		size_t x = 0, square_width = 0;
		Vecu number_of_operation = number_of_grid_points_;
		/* The smoothing coefficient is a multiple of the grid length,
		which determines the minimum radius of curvature of the levelset that needs to be smoothed */
		Real coe_ = smoothing_coe;
		
		while (x <= 200)
		{
			parallel_for(blocked_range3d<size_t>
				(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1, 1, number_of_operation[2] - 1),
				[&](const blocked_range3d<size_t>& r) 
			{
				for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
				{
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					{
						for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
						{

							/*for (size_t i = 1; i < number_of_operation[0] - 1; i++)
							{
								for (size_t j = 1; j < number_of_operation[1] - 1; j++)
								{
									for (size_t k = 1; k < number_of_operation[2] - 1; k++)
									{*/
							Real grad_x = 0.5 * (mesh_background_data_[i + 1][j][k].phi_ - mesh_background_data_[i - 1][j][k].phi_) / grid_spacing_;
							Real grad_y = 0.5 * (mesh_background_data_[i][j + 1][k].phi_ - mesh_background_data_[i][j - 1][k].phi_) / grid_spacing_;
							Real grad_z = 0.5 * (mesh_background_data_[i][j][k + 1].phi_ - mesh_background_data_[i][j][k - 1].phi_) / grid_spacing_;

							Real grad_xy = 0.25 * (mesh_background_data_[i + 1][j + 1][k].phi_ - mesh_background_data_[i - 1][j + 1][k].phi_
								- mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j - 1][k].phi_)
								/ grid_spacing_ / grid_spacing_;
							Real grad_xz = 0.25 * (mesh_background_data_[i + 1][j][k + 1].phi_ - mesh_background_data_[i - 1][j][k + 1].phi_
								- mesh_background_data_[i + 1][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_)
								/ grid_spacing_ / grid_spacing_;
							Real grad_yz = 0.25 * (mesh_background_data_[i][j + 1][k + 1].phi_ - mesh_background_data_[i][j + 1][k - 1].phi_
								- mesh_background_data_[i][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_)
								/ grid_spacing_ / grid_spacing_;

							Real grad_xx = (mesh_background_data_[i + 1][j][k].phi_
								- 2.0 * mesh_background_data_[i][j][k].phi_
								+ mesh_background_data_[i - 1][j][k].phi_) / grid_spacing_ / grid_spacing_;
							Real grad_yy = (mesh_background_data_[i][j + 1][k].phi_ -
								2.0 * mesh_background_data_[i][j][k].phi_ +
								mesh_background_data_[i][j - 1][k].phi_) / grid_spacing_ / grid_spacing_;
							Real grad_zz = (mesh_background_data_[i][j][k + 1].phi_ -
								2.0 * mesh_background_data_[i][j][k].phi_ +
								mesh_background_data_[i][j][k - 1].phi_) / grid_spacing_ / grid_spacing_;


							Real grad_phi = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;


							Real kappa = (grad_x * grad_x *grad_yy - 2.0 * grad_x * grad_y * grad_xy + grad_y * grad_y *grad_xx +
								grad_x * grad_x *grad_zz - 2.0 * grad_x * grad_z * grad_xz + grad_z * grad_z *grad_xx + grad_y * grad_y *grad_zz -
								2.0 * grad_y * grad_z * grad_yz + grad_z * grad_z *grad_yy) / (grad_phi * sqrt(grad_phi) + 1.0e-15);

							mesh_background_data_[i][j][k].kappa_ = kappa;

							/* (min, grid_spacing), (max, grid_spacing)	*/
							mesh_background_data_[i][j][k].phi_temp_min = dtau * SMIN(kappa, 0.0) * lambda * sqrt(grad_phi);
							mesh_background_data_[i][j][k].phi_temp_max = dtau * SMAX(kappa, 0.0) * lambda * sqrt(grad_phi);

						}
					}
				}
			}, ap);

			parallel_for(blocked_range3d<size_t>
				(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1, 1, number_of_operation[2] - 1),
				[&](const blocked_range3d<size_t>& r)
			{
				for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
				{
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					{
						for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
						{
			
			/*for (size_t i = 1; i < number_of_operation[0] - 2; i++)
			{
				for (size_t j = 1; j < number_of_operation[1] - 2; j++)
				{
					for (size_t k = 1; k < number_of_operation[2] - 2; k++)
					{
						if (mesh_background_data_[i][j][k].gamma_S_)
						{*/
							Vecu position = Vecu(i, j, k);
							/*if levelset is negative inside the geometric*/
							if (AverageLevelSetValueForStencil(position, square_width) > 0)
							{
								if (mesh_background_data_[i][j][k].kappa_ > (1 / (coe_ * grid_spacing_))) 
								{
									//time step 
									mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j][k].phi_temp_min;
								}
							}
							else
							{
								if (mesh_background_data_[i][j][k].kappa_ < (-1 / (coe_ * grid_spacing_)))
								{
									mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j][k].phi_temp_max;
								}
							}

							///*if levelset is positive inside the geometric*/
							//if (AverageLevelSetValueForStencil(position, square_width) < 0)
							//{
							//	if (mesh_background_data_[i][j][k].kappa_ < (-1 / (coe_ * grid_spacing_)))
							//	{
							//		//time step 
							//		mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j][k].phi_temp_min;
							//	}
							//}
							//else
							//{
							//	if (mesh_background_data_[i][j][k].kappa_ > (1 / (coe_ * grid_spacing_)))
							//	{
							//		mesh_background_data_[i][j][k].phi_ = mesh_background_data_[i][j][k].phi_ + mesh_background_data_[i][j][k].phi_temp_max;
							//	}
							//}
						}
					}
				}
			}, ap);


			x = x + 1;
		}
	}

	//===================================================================//

	Real MeshBackground::AverageLevelSetValueForStencil(Vecu position, int square_width)
	{
		Real average_phi;
		Vecu pos = position;
		int stencil_width = square_width; // useless variable, too big if square_width > 0 

		
		Real phi_1_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] + 1][pos[2]].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1]][pos[2] + 1].phi_ + mesh_background_data_[pos[0]][pos[1]][pos[2] + 1].phi_ +
				mesh_background_data_[pos[0]][pos[1] + 1][pos[2] + 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1] + 1][pos[2] + 1].phi_);
		Real phi_2_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] + 1][pos[2]].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] - 1].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2] - 1].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] + 1][pos[2] - 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2] - 1].phi_);
		Real phi_3_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] - 1][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] - 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2] - 1].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] - 1][pos[2] - 1].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2] - 1].phi_);
		Real phi_4_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] - 1][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] + 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1]][pos[2] + 1].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] - 1][pos[2] + 1].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2] + 1].phi_);
		Real phi_5_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] + 1][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] + 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2] + 1].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] + 1][pos[2] + 1].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2] + 1].phi_);
		Real phi_6_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] + 1][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] - 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2] - 1].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] + 1][pos[2] - 1].phi_ + mesh_background_data_[pos[0]][pos[1] + 1][pos[2] - 1].phi_);
		Real phi_7_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2]].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] - 1][pos[2]].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] - 1].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2] - 1].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] - 1][pos[2] - 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2] - 1].phi_);
		Real phi_8_ = (1.0 / 8.0) * (mesh_background_data_[pos[0]][pos[1]][pos[2]].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2]].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] - 1][pos[2]].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2]].phi_ +
				mesh_background_data_[pos[0]][pos[1]][pos[2] + 1].phi_ + mesh_background_data_[pos[0]][pos[1] - 1][pos[2] + 1].phi_ +
				mesh_background_data_[pos[0] - 1][pos[1] - 1][pos[2] + 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]][pos[2] + 1].phi_);


		average_phi = 0.25 * (phi_1_ + phi_2_ + phi_3_ + phi_4_ + phi_5_ + phi_6_ + phi_7_ + phi_8_);
		
		return average_phi;
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
		Vecu number_of_operation = number_of_grid_points_;
		for (size_t k = 0; k != number_of_operation[2]; ++k)
		{
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					mesh_background_data_[i][j][k].kappa_ = 0.0;
				}
				
			}
		}

		parallel_for(blocked_range3d<size_t>
			(1, number_of_operation[0]-1, 1, number_of_operation[1]-1, 1, number_of_operation[2]-1),
			[&](const blocked_range3d<size_t>& r) {
			for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
				for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
					{
						Real grad_x = 0.5 * (mesh_background_data_[i + 1][j][k].phi_ - mesh_background_data_[i - 1][j][k].phi_) / grid_spacing_;
						Real grad_y = 0.5 * (mesh_background_data_[i][j + 1][k].phi_ - mesh_background_data_[i][j - 1][k].phi_) / grid_spacing_;
						Real grad_z = 0.5 * (mesh_background_data_[i][j][k + 1].phi_ - mesh_background_data_[i][j][k - 1].phi_) / grid_spacing_;

						Real grad_xy = 0.25 * (mesh_background_data_[i + 1][j + 1][k].phi_ - mesh_background_data_[i - 1][j + 1][k].phi_
							- mesh_background_data_[i + 1][j - 1][k].phi_ + mesh_background_data_[i - 1][j - 1][k].phi_)
							/ grid_spacing_ / grid_spacing_;
						Real grad_xz = 0.25 * (mesh_background_data_[i + 1][j][k + 1].phi_ - mesh_background_data_[i - 1][j][k + 1].phi_
							- mesh_background_data_[i + 1][j][k - 1].phi_ + mesh_background_data_[i - 1][j][k - 1].phi_)
							/ grid_spacing_ / grid_spacing_;
						Real grad_yz = 0.25 * (mesh_background_data_[i][j + 1][k + 1].phi_ - mesh_background_data_[i][j + 1][k - 1].phi_
							- mesh_background_data_[i][j - 1][k + 1].phi_ + mesh_background_data_[i][j - 1][k - 1].phi_)
							/ grid_spacing_ / grid_spacing_;

						Real grad_xx = (mesh_background_data_[i + 1][j][k].phi_
							- 2.0 * mesh_background_data_[i][j][k].phi_
							+ mesh_background_data_[i - 1][j][k].phi_) / grid_spacing_ / grid_spacing_;
						Real grad_yy = (mesh_background_data_[i][j + 1][k].phi_ -
							2.0 * mesh_background_data_[i][j][k].phi_ +
							mesh_background_data_[i][j - 1][k].phi_) / grid_spacing_ / grid_spacing_;
						Real grad_zz = (mesh_background_data_[i][j][k + 1].phi_ -
							2.0 * mesh_background_data_[i][j][k].phi_ +
							mesh_background_data_[i][j][k - 1].phi_) / grid_spacing_ / grid_spacing_;

						Real grad_phi = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;

						mesh_background_data_[i][j][k].kappa_ = (grad_x * grad_x *grad_yy - 2.0 * grad_x * grad_y * grad_xy + grad_y * grad_y *grad_xx +
							grad_x * grad_x *grad_zz - 2.0 * grad_x * grad_z * grad_xz + grad_z * grad_z *grad_xx + grad_y * grad_y *grad_zz - 
							2.0 * grad_y * grad_z * grad_yz + grad_z * grad_z *grad_yy) / (grad_phi * sqrt(grad_phi) + 1.0e-15);
					}
		}, ap);
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
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (size_t i = 0; i != 4; ++i)
			for (size_t j = 0; j != 4; ++j)
				for (size_t k = 0; k != 4; ++k)
				{
					pkg_data_[i][j][k].phi_ = level_set;
					pkg_data_[i][j][k].n_ = normal_direction;
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
					sph_body->ClosestPointOnBodySurface(position, closet_pnt_on_face, pkg_data_[i][j][k].phi_);
					pkg_data_[i][j][k].n_ = closet_pnt_on_face - position;
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
