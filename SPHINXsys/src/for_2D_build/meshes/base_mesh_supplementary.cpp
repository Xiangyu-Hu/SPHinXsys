/**
 * @file 	base_mesh_supplementary.cpp
 * @author	Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
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
//=================================================================================================//
namespace SPH {
	//=============================================================================================//
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt)
	{
		for (size_t i = index_begin[0]; i != index_end[0]; ++i)
			for (size_t j = index_begin[1]; j != index_end[1]; ++j) {
				mesh_functor(Vecu(i, j), dt);
			}
	}
	//=============================================================================================//
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
	//=============================================================================================//
	Vecu BaseMesh::transfer1DtoMeshIndex(Vecu mesh_size, size_t i)
	{
		int row_size = mesh_size[1];
		int column = i / row_size;
		return Vecu(column, i - column * row_size);
	}
	//=============================================================================================//
	size_t BaseMesh::transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index)
	{
		return mesh_index[0] * mesh_size[1] + mesh_index[1];
	}
	//=============================================================================================//
	void MeshBackground::AllocateMeshDataMatrix()
	{
		Allocate2dArray(mesh_background_data_, number_of_grid_points_);
	}
	//=============================================================================================//
	void MeshBackground::DeleteMeshDataMatrix()
	{
		Delete2dArray(mesh_background_data_, number_of_grid_points_);
	}
	//=============================================================================================//
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
				}
		}, ap);

		parallel_for(blocked_range2d<size_t>
			(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
					{
						Real dphidx = mesh_background_data_[i + 1][j].phi_ - mesh_background_data_[i - 1][j].phi_;
						Real dphidy = mesh_background_data_[i][j + 1].phi_ - mesh_background_data_[i][j - 1].phi_;
						Vecd normal = Vecd(dphidx, dphidy);
						mesh_background_data_[i][j].n_ = normal / (normal.norm() + 1.0e-16);
					}
			}, ap);
	}
	//=============================================================================================//
	void MeshBackground::ClearLevelSetData(SPHBody &body)
	{
		Real epsilon_ = 0.75*grid_spacing_;
		std::vector<Vecu> position_S_ = GetZeroLevelSetCutCell();
		std::vector<Vecu> position_P_ = GetPositiveCutCell();
		std::vector<Vecu> position_N_ = GetNegativeCutCell();
		Vecu number_of_operation = number_of_grid_points_;

		MarkZeroLevelSetCutCell();
		for (size_t x = 0; x < position_S_.size(); x++)
		{
			if (IfNeighborCellsInAuxiliayBandP(position_S_[x]) == 0) //function equal to zero means no neighborcells in band P 
			{
				std::vector<Real>d_P_;
				for (int i = 0; i < 9; i++)
				{
					for (int j = 0; j < 9; j++) {

						if (mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].gamma_P_)
						{
							Real phi_P_ = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].phi_;
							Vecd norm_to_face = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].n_
								/ mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].n_.norm();
							Vec2d grid_pos_S = GridPositionFromIndexes(position_S_[x]);
							Vec2d grid_pos_P = GridPositionFromIndexes(Vec2u((position_S_[x][0] - 4 + i), (position_S_[x][1] - 4 + j)));

							d_P_.push_back(sqrt(powern((grid_pos_S[0] - grid_pos_P[0] + (phi_P_ - epsilon_)*norm_to_face[0]), 2) +
								powern((grid_pos_S[1] - grid_pos_P[1] + (phi_P_ - epsilon_)*norm_to_face[1]), 2)));
						}
					}
				}
				Real d_min_P_;
				if (d_P_.size() != 0) {
					d_min_P_ = *min_element(d_P_.begin(), d_P_.end());
				}
				else {
					d_min_P_ = 3.0 * grid_spacing_;
				}
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]].phi_ = -(d_min_P_ - epsilon_);
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]].gamma_S_ = false;
			}

			if (IfNeighborCellsInAuxiliayBandN(position_S_[x]) == 0)  //function equal to zero means no neighborcells in band N 
			{
				std::vector<Real>d_N_;
				for (int i = 0; i < 9; i++)
				{
					for (int j = 0; j < 9; j++) {
						if (mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].gamma_N_)
						{
							Real phi_N_ = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].phi_;
							Vecd norm_to_face = mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].n_
								/ mesh_background_data_[position_S_[x][0] - 4 + i][position_S_[x][1] - 4 + j].n_.norm();
							Vec2d grid_pos_S = GridPositionFromIndexes(position_S_[x]);
							Vec2d grid_pos_N = GridPositionFromIndexes(Vec2u((position_S_[x][0] - 4 + i), (position_S_[x][1] - 4 + j)));

							d_N_.push_back(sqrt(powern((grid_pos_S[0] - grid_pos_N[0] + (phi_N_ + epsilon_)*norm_to_face[0]), 2) +
								powern((grid_pos_S[1] - grid_pos_N[1] + (phi_N_ + epsilon_)*norm_to_face[1]), 2)));
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
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]].phi_ = d_min_N_ + epsilon_;
				mesh_background_data_[position_S_[x][0]][position_S_[x][1]].gamma_S_ = false;
			}
		}
	}
	//=============================================================================================//
	void MeshBackground::MarkZeroLevelSetCutCell()
	{
		Vecu number_of_operation = number_of_grid_points_;
		parallel_for(blocked_range2d<size_t>
			(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1),
			[&](const blocked_range2d<size_t>& r)
		{
			for (size_t i = r.rows().begin(); i != r.rows().end(); i++)
			{
				for (size_t j = r.cols().begin(); j != r.cols().end(); j++)
				{
					Real phi_1_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
						mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i + 1][j + 1].phi_);
					Real phi_2_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
						mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i + 1][j - 1].phi_);
					Real phi_3_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
						mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i - 1][j + 1].phi_);
					Real phi_4_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
						mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i - 1][j - 1].phi_);

					if (phi_1_ * mesh_background_data_[i][j].phi_ <= 0) {
						mesh_background_data_[i][j].gamma_S_ = true;
					}
					else if (phi_2_ * mesh_background_data_[i][j].phi_ <= 0) {
						mesh_background_data_[i][j].gamma_S_ = true;
					}
					else if (phi_3_ * mesh_background_data_[i][j].phi_ <= 0) {
						mesh_background_data_[i][j].gamma_S_ = true;
					}
					else if (phi_4_ * mesh_background_data_[i][j].phi_ <= 0) {
						mesh_background_data_[i][j].gamma_S_ = true;
					}
					else {
						mesh_background_data_[i][j].gamma_S_ = false;
					}
				}
			}

		}, ap);
	}
	//=============================================================================================//
	vector<Vecu> MeshBackground::GetZeroLevelSetCutCell()
	{
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_S_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				Real phi_1_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i + 1][j + 1].phi_);
				Real phi_2_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i + 1][j - 1].phi_);
				Real phi_3_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i - 1][j + 1].phi_);
				Real phi_4_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i - 1][j - 1].phi_);

				if (phi_1_ * mesh_background_data_[i][j].phi_ <= 0) {
					position_S_.push_back(Vecu(i, j));
					mesh_background_data_[i][j].gamma_S_ = true;
				}
				else if (phi_2_ * mesh_background_data_[i][j].phi_ <= 0) {
					position_S_.push_back(Vecu(i, j));
					mesh_background_data_[i][j].gamma_S_ = true;
				}
				else if (phi_3_ * mesh_background_data_[i][j].phi_ <= 0) {
					position_S_.push_back(Vecu(i, j));
					mesh_background_data_[i][j].gamma_S_ = true;
				}
				else if (phi_4_ * mesh_background_data_[i][j].phi_ <= 0) {
					position_S_.push_back(Vecu(i, j));
					mesh_background_data_[i][j].gamma_S_ = true;
				}
			}
		}
		return position_S_;
	}
	//=============================================================================================//
	vector<Vecu> MeshBackground::GetPositiveCutCell()
	{
		Real epsilon_ = 0.75*grid_spacing_;
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_P_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				Real phi_1_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i + 1][j + 1].phi_);
				Real phi_2_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i + 1][j - 1].phi_);
				Real phi_3_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i - 1][j + 1].phi_);
				Real phi_4_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i - 1][j - 1].phi_);

				if ((phi_1_ - epsilon_) * (mesh_background_data_[i][j].phi_ - epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_P_ = true;
					position_P_.push_back(Vecu(i, j));
				}
				else if ((phi_2_ - epsilon_) * (mesh_background_data_[i][j].phi_ - epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_P_ = true;
					position_P_.push_back(Vecu(i, j));
				}
				else if ((phi_3_ - epsilon_) * (mesh_background_data_[i][j].phi_ - epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_P_ = true;
					position_P_.push_back(Vecu(i, j));
				}
				else if ((phi_4_ - epsilon_) * (mesh_background_data_[i][j].phi_ - epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_P_ = true;
					position_P_.push_back(Vecu(i, j));
				}
			}
		}
		return position_P_;
	}
	//=============================================================================================//
	vector<Vecu> MeshBackground::GetNegativeCutCell()
	{
		Real epsilon_ = 0.75*grid_spacing_;
		Vecu number_of_operation = number_of_grid_points_;
		std::vector<Vecu> position_N_;
		for (int i = 1; i < number_of_operation[0] - 1; i++)
		{
			for (int j = 1; j < number_of_operation[1] - 1; j++)
			{
				Real phi_1_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i + 1][j + 1].phi_);
				Real phi_2_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i + 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i + 1][j - 1].phi_);
				Real phi_3_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j + 1].phi_ + mesh_background_data_[i - 1][j + 1].phi_);
				Real phi_4_ = (1.0 / 4.0) * (mesh_background_data_[i][j].phi_ + mesh_background_data_[i - 1][j].phi_ +
					mesh_background_data_[i][j - 1].phi_ + mesh_background_data_[i - 1][j - 1].phi_);

				if ((phi_1_ + epsilon_) * (mesh_background_data_[i][j].phi_ + epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_N_ = true;
					position_N_.push_back(Vecu(i, j));
				}
				else if ((phi_2_ + epsilon_) * (mesh_background_data_[i][j].phi_ + epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_N_ = true;
					position_N_.push_back(Vecu(i, j));
				}
				else if ((phi_3_ + epsilon_) * (mesh_background_data_[i][j].phi_ + epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_N_ = true;
					position_N_.push_back(Vecu(i, j));
				}
				else if ((phi_4_ + epsilon_) * (mesh_background_data_[i][j].phi_ + epsilon_) <= 0) {
					mesh_background_data_[i][j].gamma_N_ = true;
					position_N_.push_back(Vecu(i, j));
				}
			}
		}
		return position_N_;
	}
	//=============================================================================================//
	int MeshBackground::IfNeighborCellsInAuxiliayBandP(Vecu position)
	{
		int p_ = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int a = (position[0] - 1 + i), b = (position[1] - 1 + j);
				bool P = mesh_background_data_[position[0] - 1 + i][position[1] - 1 + j].gamma_P_;
				if (mesh_background_data_[position[0] - 1 + i][position[1] - 1 + j].gamma_P_) 
				{
					p_ = p_ + 1;
				};
			}
		}
		return p_;
	}
	//=============================================================================================//
	int MeshBackground::IfNeighborCellsInAuxiliayBandN(Vecu position)
	{
		int n_ = 0;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (mesh_background_data_[position[0] - 1 + i][position[1] - 1 + j].gamma_N_)
				{
					n_ = n_ + 1;
				};
			}
		}
		return n_;
	}
	//=============================================================================================//
	void MeshBackground::ReinitializeLevelSetData(SPHBody &body)
	{
		Real dtau = 0.5*grid_spacing_;
		Real s, ss;
		Real dv_xp, dv_x, dv_xn;
		Real dv_yp, dv_y, dv_yn;
		size_t k = 0;
		Vecu number_of_operation = number_of_grid_points_;
		
		while (k <= 50)
		{
			for (size_t i = 1; i < number_of_operation[0] - 1; i++)
			{
				for (size_t j = 1; j < number_of_operation[1] - 1; j++)
				{
					if (mesh_background_data_[i][j].gamma_S_ != true)
					{

						s = mesh_background_data_[i][j].phi_ /
							sqrt(mesh_background_data_[i][j].phi_*mesh_background_data_[i][j].phi_ + grid_spacing_ * grid_spacing_);
						dv_xp = (mesh_background_data_[i + 1][j].phi_ - mesh_background_data_[i][j].phi_);
						dv_xn = (mesh_background_data_[i][j].phi_ - mesh_background_data_[i - 1][j].phi_);
						dv_x = dv_xp;
						if (s * dv_xp >= 0.0 && s * dv_xn >= 0.0) dv_x = dv_xn;
						if (s * dv_xp <= 0.0 && s * dv_xn <= 0.0) dv_x = dv_xp;
						if (s * dv_xp > 0.0 && s * dv_xn < 0.0) dv_x = 0.0;
						if (s * dv_xp < 0.0 && s * dv_xn > 0.0)
						{
							ss = s * (fabs(dv_xp) - fabs(dv_xn)) / (dv_xp - dv_xn);
							if (ss > 0.0) dv_x = dv_xn;
						}

						dv_yp = (mesh_background_data_[i][j + 1].phi_ - mesh_background_data_[i][j].phi_);
						dv_yn = (mesh_background_data_[i][j].phi_ - mesh_background_data_[i][j - 1].phi_);
						dv_y = dv_yp;
						if (s * dv_yp >= 0.0 && s * dv_yn >= 0.0) dv_y = dv_yn;
						if (s * dv_yp <= 0.0 && s * dv_yn <= 0.0) dv_y = dv_yp;
						if (s * dv_yp > 0.0 && s * dv_yn < 0.0) dv_y = 0.0;
						if (s * dv_yp < 0.0 && s * dv_yn > 0.0)
						{
							ss = s * (fabs(dv_yp) - fabs(dv_yn)) / (dv_yp - dv_yn);
							if (ss > 0.0) dv_y = dv_yn;
						}

						//time step 
						mesh_background_data_[i][j].phi_ = mesh_background_data_[i][j].phi_ -
							dtau * s * (sqrt(dv_x * dv_x / grid_spacing_ / grid_spacing_ + dv_y * dv_y / grid_spacing_ / grid_spacing_) - 1.0);
					}
				}
			}
			
			//update bounday data
			for (size_t i = 0; i < 1; i++) {
				for (size_t j = 0; j < number_of_operation[1]; j++)
				{
					mesh_background_data_[i][j].phi_ = mesh_background_data_[1][j].phi_;
				}
			}

			for (size_t i = number_of_operation[0] - 1; i < number_of_operation[0]; i++) {
				for (size_t j = 0; j < number_of_operation[1]; j++)
				{
					mesh_background_data_[i][j].phi_ = mesh_background_data_[number_of_operation[0] - 2][j].phi_;
				}
			}

			for (size_t j = 0; j < 1; j++) {
				for (size_t i = 1; i < number_of_operation[0] - 1; i++)
				{
					mesh_background_data_[i][j].phi_ = mesh_background_data_[i][1].phi_;
				}
			}

			for (size_t j = number_of_operation[1] - 1; j < number_of_operation[1]; j++) {
				for (size_t i = 1; i < number_of_operation[0] - 1; i++)
				{
					mesh_background_data_[i][j].phi_ = mesh_background_data_[i][number_of_operation[1] - 2].phi_;
				}
			}
			k = k + 1;
		}

		//update the normal direction
		for (size_t i = 2; i <= number_of_operation[0] - 2; i++)
		{
			for (size_t j = 2; j <= number_of_operation[1] - 2; j++)
			{
				mesh_background_data_[i][j].n_[0] = ((-1) * SGN(mesh_background_data_[i + 1][j].phi_) * mesh_background_data_[i + 1][j].phi_
					- (-1) * SGN(mesh_background_data_[i - 1][j].phi_) * mesh_background_data_[i - 1][j].phi_) / (2.0 * grid_spacing_);
				mesh_background_data_[i][j].n_[1] = ((-1) * SGN(mesh_background_data_[i][j + 1].phi_) * mesh_background_data_[i][j + 1].phi_
					- (-1) * SGN(mesh_background_data_[i][j - 1].phi_) * mesh_background_data_[i][j - 1].phi_) / (2.0 * grid_spacing_);
			}
		}
	}
	//=============================================================================================//
	void MeshBackground::SmoothLevelSetByCurvature(SPHBody &body, Real smoothing_coe)
	{
		Real lambda = 0.02;
		Real dtau = 0.02*grid_spacing_;
		size_t k = 0, square_width = 0;
		Vecu number_of_operation = number_of_grid_points_;
		/* The smoothing coefficient is a multiple of the grid length, 
		which determines the minimum radius of curvature of the levelset that needs to be smoothed */
		Real coe_ = smoothing_coe;

		while (k <= 200)
		{
			parallel_for(blocked_range2d<size_t>
				(1, number_of_operation[0] - 1, 1, number_of_operation[1] - 1),
				[&](const blocked_range2d<size_t>& r)
				{
					for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					{
						for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
						{
							Real grad_x = 0.5 * (mesh_background_data_[i + 1][j].phi_ - mesh_background_data_[i - 1][j].phi_) / grid_spacing_;
							Real grad_y = 0.5 * (mesh_background_data_[i][j + 1].phi_ - mesh_background_data_[i][j - 1].phi_) / grid_spacing_;
							Real grad_x1 = (mesh_background_data_[i + 1][j].phi_ - mesh_background_data_[i][j].phi_) / grid_spacing_;
							Real grad_y1 = (mesh_background_data_[i][j + 1].phi_ - mesh_background_data_[i][j].phi_) / grid_spacing_;
	
							Real grad_x2 = (mesh_background_data_[i][j].phi_ - mesh_background_data_[i - 1][j].phi_) / grid_spacing_;
							Real grad_y2 = (mesh_background_data_[i][j].phi_ - mesh_background_data_[i][j - 1].phi_) / grid_spacing_;

							Real delt_plus = sqrt(max(grad_x2, 0.0)*max(grad_x2, 0.0) + min(grad_x1, 0.0)*min(grad_x1, 0.0) + max(grad_y2, 0.0)*max(grad_y2, 0.0) + min(grad_y1, 0.0)*min(grad_y1, 0.0));
							Real delt_minus = sqrt(min(grad_x2, 0.0)*min(grad_x2, 0.0) + max(grad_x1, 0.0)*max(grad_x1, 0.0) + min(grad_y2, 0.0)*min(grad_y2, 0.0) + max(grad_y1, 0.0)*max(grad_y1, 0.0));

							Real grad_xy = 0.25 * (mesh_background_data_[i + 1][j + 1].phi_ - mesh_background_data_[i - 1][j + 1].phi_
								- mesh_background_data_[i + 1][j - 1].phi_ + mesh_background_data_[i - 1][j - 1].phi_)
								/ grid_spacing_ / grid_spacing_;
							Real grad_xx = (mesh_background_data_[i + 1][j].phi_
								- 2.0 * mesh_background_data_[i][j].phi_
								+ mesh_background_data_[i - 1][j].phi_) / grid_spacing_ / grid_spacing_;
							Real grad_yy = (mesh_background_data_[i][j + 1].phi_ -
								2.0 * mesh_background_data_[i][j].phi_ +
								mesh_background_data_[i][j - 1].phi_) / grid_spacing_ / grid_spacing_;

							Real grad_phi = grad_x * grad_x + grad_y * grad_y;

							Real kappa = (grad_xx * grad_y * grad_y - 2.0 * grad_x * grad_y * grad_xy +
								grad_yy * grad_x * grad_x) / (grad_phi * sqrt(grad_phi) + sqrt(grid_spacing_) * grid_spacing_);
							mesh_background_data_[i][j].kappa_ = kappa;

							/* (min, grid_spacing), (max, grid_spacing)	*/
							mesh_background_data_[i][j].phi_temp_min = dtau * SMIN(kappa, 0.0) * lambda * sqrt(grad_phi);
							mesh_background_data_[i][j].phi_temp_max = dtau * SMAX(kappa, 0.0) * lambda * sqrt(grad_phi);
						}
					}
			}, ap);

			parallel_for(blocked_range2d<size_t>
				(1, number_of_operation[0] - 2, 1, number_of_operation[1] - 2),
				[&](const blocked_range2d<size_t>& r)
				{
					for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					{
						for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
						{
							if (mesh_background_data_[i][j].gamma_S_)
							{
								Vecu position = Vecu(i, j);
								/*if levelset is negative inside the geometric*/
								if (AverageLevelSetValueForStencil(position, square_width) > 0)
								{
									if (mesh_background_data_[i][j].kappa_ > (1 / (coe_ * grid_spacing_)))
									{
										//time step 
										mesh_background_data_[i][j].phi_ = mesh_background_data_[i][j].phi_ + mesh_background_data_[i][j].phi_temp_min;
									}
								}
								else
								{
									if (mesh_background_data_[i][j].kappa_ < (-1 / (coe_ * grid_spacing_)))
									{
										mesh_background_data_[i][j].phi_ = mesh_background_data_[i][j].phi_ + mesh_background_data_[i][j].phi_temp_max;
									}
								}

							}
						}
					}
				}, ap);
			
			k = k + 1;
		}
	}
	//=============================================================================================//
	Real MeshBackground::AverageLevelSetValueForStencil(Vecu position, int square_width)
	{
		Real average_phi;
		Vecu pos = position;
		int stencil_width = square_width;

		if (stencil_width == 1)
		{
			average_phi = 1 / ((2.0 * Real(stencil_width) + 1) * (2.0 * Real(stencil_width) + 1)) *
				(mesh_background_data_[pos[0] - 1][pos[1] - 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]].phi_ +
					mesh_background_data_[pos[0] - 1][pos[1] + 1].phi_ + mesh_background_data_[pos[0]][pos[1] - 1].phi_ +
					mesh_background_data_[pos[0]][pos[1]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1].phi_ +
					mesh_background_data_[pos[0] + 1][pos[1] - 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1]].phi_ +
					mesh_background_data_[pos[0] + 1][pos[1] + 1].phi_);
		}

		if (stencil_width == 0)
		{
			Real phi_1 = 0.25 * (mesh_background_data_[pos[0] - 1][pos[1] - 1].phi_ + mesh_background_data_[pos[0] - 1][pos[1]].phi_ +
				mesh_background_data_[pos[0]][pos[1] - 1].phi_ + mesh_background_data_[pos[0]][pos[1]].phi_);
			Real phi_2 = 0.25 * (mesh_background_data_[pos[0] - 1][pos[1]].phi_ + mesh_background_data_[pos[0] - 1][pos[1] + 1].phi_ +
				mesh_background_data_[pos[0]][pos[1]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1].phi_);
			Real phi_3 = 0.25 * (mesh_background_data_[pos[0]][pos[1] - 1].phi_ + mesh_background_data_[pos[0]][pos[1]].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1] - 1].phi_ + mesh_background_data_[pos[0] + 1][pos[1]].phi_);
			Real phi_4 = 0.25 * (mesh_background_data_[pos[0]][pos[1]].phi_ + mesh_background_data_[pos[0]][pos[1] + 1].phi_ +
				mesh_background_data_[pos[0] + 1][pos[1]].phi_ + mesh_background_data_[pos[0] + 1][pos[1] + 1].phi_);

			average_phi = 0.25 * (phi_1 + phi_2 + phi_3 + phi_4);
		}
		return average_phi;
	}
	//=============================================================================================//
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
	//=============================================================================================//
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
	//=============================================================================================//
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
	//=============================================================================================//
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
	//=============================================================================================//
	void MeshBackground::WriteMeshToVtuFile(ofstream &output_file)
	{
		cout << "\n This function WriteMeshToVtuFile is not done. Exit the program! \n";
		exit(0);

	}
	//=============================================================================================//
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
	//=============================================================================================//
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (size_t i = 0; i != pkg_size_; ++i)
			for (size_t j = 0; j != pkg_size_; ++j) {
				pkg_data_[i][j].phi_ = level_set;
				pkg_data_[i][j].n_ = normal_direction;
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
				sph_body->ClosestPointOnBodySurface(position, closet_pnt_on_face, pkg_data_[i][j].phi_);
				pkg_data_[i][j].n_ = closet_pnt_on_face - position;
			}
	}
	//=============================================================================================//
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
	//=============================================================================================//
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
	//=============================================================================================//
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
				output_file << getValueFromGlobalDataIndex<LevelSetData>(Vecu(i, j)).phi_ << " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << getValueFromGlobalDataIndex<LevelSetData>(Vecu(i, j)).n_[0]<< " ";

			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << getValueFromGlobalDataIndex<LevelSetData>(Vecu(i, j)).n_[1] << " ";

			}
			output_file << " \n";
		}
	}
	//=============================================================================================//
}
//=============================================================================================//
