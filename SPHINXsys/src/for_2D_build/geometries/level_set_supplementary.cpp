/**
 * @file 	level_set_supplementary.cpp
 * @author	Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
 * @version	0.1
 */

#include "level_set.h"
#include "base_particles.h"
#include "base_body.h"

//=================================================================================================//
namespace SPH {
	//=============================================================================================//
	void LevelSetDataPackage::initializeWithUniformData(Real level_set, Vecd normal_direction)
	{
		for (int i = 0; i != PackageSize(); ++i)
			for (int j = 0; j != PackageSize(); ++j) {
				phi_[i][j] = level_set;
				n_[i][j] = normal_direction;
				near_interface_id_[i][j]= level_set < 0.0 ? -2 : 2;
			}
	}
	//=================================================================================================//
	void  LevelSetDataPackage::initializeDataPackage(ComplexShape& complex_shape)
	{
		for (int i = 0; i != PackageSize(); ++i)
			for (int j = 0; j != PackageSize(); ++j)
			{
				Vec2d position = data_lower_bound_ + Vec2d((Real)i * grid_spacing_, (Real)j * grid_spacing_);
				Real distance_from_surface
					= (complex_shape.findClosestPoint(position) - position).norm();
				phi_[i][j] = complex_shape.checkContain(position) ?
					-distance_from_surface : distance_from_surface;
			}
	}
	//=================================================================================================//
	void LevelSetDataPackage::stepReinitialization()
	{
		for (int i = AddressBufferWidth(); i != OperationUpperBound(); ++i)
			for (int j = AddressBufferWidth(); j != OperationUpperBound(); ++j)
			{
				//only reinitialize non cut cells
				if (*near_interface_id_addrs_[i][j] != 0)
				{
					Real phi_0 = *phi_addrs_[i][j];
					Real s = phi_0 / sqrt(phi_0 * phi_0 + grid_spacing_ * grid_spacing_);
					//x direction
					Real dv_xp = (*phi_addrs_[i + 1][j] - phi_0);
					Real dv_xn = (phi_0 - *phi_addrs_[i - 1][j]);
					Real dv_x = dv_xp;
					if (s * dv_xp >= 0.0 && s * dv_xn >= 0.0) dv_x = dv_xn;
					if (s * dv_xp <= 0.0 && s * dv_xn <= 0.0) dv_x = dv_xp;
					if (s * dv_xp > 0.0 && s * dv_xn < 0.0) dv_x = 0.0;
					if (s * dv_xp < 0.0 && s * dv_xn > 0.0)
					{
						Real ss = s * (fabs(dv_xp) - fabs(dv_xn)) / (dv_xp - dv_xn);
						if (ss > 0.0) dv_x = dv_xn;
					}
					//y direction
					Real dv_yp = (*phi_addrs_[i][j + 1] - phi_0);
					Real dv_yn = (phi_0 - *phi_addrs_[i][j - 1]);
					Real dv_y = dv_yp;
					if (s * dv_yp >= 0.0 && s * dv_yn >= 0.0) dv_y = dv_yn;
					if (s * dv_yp <= 0.0 && s * dv_yn <= 0.0) dv_y = dv_yp;
					if (s * dv_yp > 0.0 && s * dv_yn < 0.0) dv_y = 0.0;
					if (s * dv_yp < 0.0 && s * dv_yn > 0.0)
					{
						Real ss = s * (fabs(dv_yp) - fabs(dv_yn)) / (dv_yp - dv_yn);
						if (ss > 0.0) dv_y = dv_yn;
					}
					//time stepping
					*phi_addrs_[i][j] -= 0.5 * s * (sqrt(dv_x * dv_x + dv_y * dv_y) - grid_spacing_);
				}
			}
	}
	//=================================================================================================//
	void LevelSetDataPackage::markNearInterface()
	{
		Real small_shift = 0.75 * grid_spacing_;
		//corner averages, note that the first row and first column are not used 
		PackageTemporaryData<Real> corner_averages;
		for (int i = 1; i != AddressSize(); ++i)
			for (int j = 1; j != AddressSize(); ++j)
			{
				corner_averages[i][j] = CornerAverage(phi_addrs_, Veci(i, j), Veci(-1, -1));
			}

		for (int i = AddressBufferWidth(); i != OperationUpperBound(); ++i)
			for (int j = AddressBufferWidth(); j != OperationUpperBound(); ++j)
			{
				//first assume far cells
				Real phi_0 = *phi_addrs_[i][j];
				int near_interface_id = phi_0 > 0.0 ? 2 : -2;

				Real phi_average_0 = corner_averages[i][j];
				//find outer cut cells by comparing the sign of corner averages
				for (int l = 0; l != 2; ++l)
					for (int m = 0; m != 2; ++m)
					{
						int index_x = i + l;
						int index_y = j + m;
						Real phi_average = corner_averages[index_x][index_y];
						if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0) near_interface_id = 1;
						if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0) near_interface_id = -1;
					}

				//find zero cut cells by comparing the sign of corner averages
				for (int l = 0; l != 2; ++l)
					for (int m = 0; m != 2; ++m)
					{
						int index_x = i + l;
						int index_y = j + m;
						Real phi_average = corner_averages[index_x][index_y];
						if (phi_average_0 * phi_average < 0.0) near_interface_id = 0;
					}

				//find cells between cut cells
				if (fabs(phi_0) < small_shift && abs(near_interface_id) != 1)  near_interface_id = 0;

				//assign this to package
				*near_interface_id_addrs_[i][j] = near_interface_id;
			}
	}
	//=============================================================================================//
	void LevelSet::initializeDataInACell(Vecu cell_index, Real dt)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		Vecd cell_position = CellPositionFromIndexes(cell_index);
		Vecd closet_pnt_on_face = complex_shape_.findClosestPoint(cell_position);
		Real measure = getMinAbsoluteElement(closet_pnt_on_face - cell_position);
		if (measure < cell_spacing_) {
			mutex_my_pool.lock();
			LevelSetDataPackage* new_data_pkg = data_pkg_pool_.malloc();
			mutex_my_pool.unlock();
			Vecd pkg_lower_bound = GridPositionFromCellPosition(cell_position);
			new_data_pkg->initializePackageGeometry(pkg_lower_bound, data_spacing_);
			new_data_pkg->initializeDataPackage(complex_shape_);
			core_data_pkgs_.push_back(new_data_pkg);
			new_data_pkg->pkg_index_ = Vecu(i, j);
			new_data_pkg->is_core_pkg_ = true;
			data_pkg_addrs_[i][j] = new_data_pkg;
		}
		else {
			data_pkg_addrs_[i][j] = complex_shape_.checkContain(cell_position) ?
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

		if (is_inner_pkg) 
		{
			LevelSetDataPackage* current_data_pkg = data_pkg_addrs_[i][j];
			if (current_data_pkg->is_core_pkg_) {
				current_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(current_data_pkg);
			}
			else {
				mutex_my_pool.lock();
				LevelSetDataPackage* new_data_pkg = data_pkg_pool_.malloc();
				mutex_my_pool.unlock();
				Vecd cell_position = CellPositionFromIndexes(cell_index);
				Vecd pkg_lower_bound = GridPositionFromCellPosition(cell_position);
				new_data_pkg->initializePackageGeometry(pkg_lower_bound, data_spacing_);
				new_data_pkg->initializeDataPackage(complex_shape_);
				new_data_pkg->pkg_index_ = Vecu(i, j);
				new_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(new_data_pkg);
				data_pkg_addrs_[i][j] = new_data_pkg;
			}
		}
	}
	//=================================================================================================//
	void LevelSet::redistanceInterfaceForAPackage(LevelSetDataPackage* core_data_pkg, Real dt)
	{
		int l = (int)core_data_pkg->pkg_index_[0];
		int m = (int)core_data_pkg->pkg_index_[1];

		for (int i = pkg_addrs_buffer_; i != pkg_operations_; ++i)
			for (int j = pkg_addrs_buffer_; j != pkg_operations_; ++j)
			{
				int near_interface_id = *core_data_pkg->near_interface_id_addrs_[i][j];
				if (near_interface_id == 0)
				{
					bool positive_band = false;
					bool negative_band = false;
					for (int s = -1; s < 2; ++s)
						for (int t = -1; t < 2; ++t)
						{
							int neighbor_near_interface_id =
								*core_data_pkg->near_interface_id_addrs_[i + s][j + t];
							if (neighbor_near_interface_id >= 1) positive_band = true;
							if (neighbor_near_interface_id <= -1) negative_band = true;
						}
					if (positive_band == false)
					{
						Real min_distance_p = 5.0 * data_spacing_;
						for (int x = -4; x != 5; ++x)
							for (int y = -4; y != 5; ++y)
							{
								pair<int, int>  x_pair = CellShiftAndDataIndex(i + x);
								pair<int, int>  y_pair = CellShiftAndDataIndex(j + y);
								LevelSetDataPackage* neighbor_pkg
									= data_pkg_addrs_[l + x_pair.first][m + y_pair.first];
								int neighbor_near_interface_id
									= neighbor_pkg->near_interface_id_[x_pair.second][y_pair.second];
								if (neighbor_near_interface_id >= 1)
								{
									Real phi_p_ = neighbor_pkg->phi_[x_pair.second][y_pair.second];
									Vecd norm_to_face = neighbor_pkg->n_[x_pair.second][y_pair.second];
									min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y) * data_spacing_ + phi_p_ * norm_to_face).norm());
								}
							}
						*core_data_pkg->phi_addrs_[i][j] = -min_distance_p;
						// this immediate switch of near interface id 
						// does not intervenning with the identification of unresolved interface
						// based on the assumption that positive false_and negative bands are not close to each other
						*core_data_pkg->near_interface_id_addrs_[i][j] = -1;
					}
					if (negative_band == false)
					{
						Real min_distance_n = 5.0 * data_spacing_;
						for (int x = -4; x != 5; ++x)
							for (int y = -4; y != 5; ++y)
							{
								pair<int, int>  x_pair = CellShiftAndDataIndex(i + x);
								pair<int, int>  y_pair = CellShiftAndDataIndex(j + y);
								LevelSetDataPackage* neighbor_pkg
									= data_pkg_addrs_[l + x_pair.first][m + y_pair.first];
								int neighbor_near_interface_id
									= neighbor_pkg->near_interface_id_[x_pair.second][y_pair.second];
								if (neighbor_near_interface_id <= -1)
								{
									Real phi_n_ = neighbor_pkg->phi_[x_pair.second][y_pair.second];
									Vecd norm_to_face = neighbor_pkg->n_[x_pair.second][y_pair.second];
									min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y) * data_spacing_ - phi_n_ * norm_to_face).norm());
								}
							}
						*core_data_pkg->phi_addrs_[i][j] = min_distance_n;
						// this immediate switch of near interface id 
						// does not intervenning with the identification of unresolved interface
						// based on the assumption that positive false_and negative bands are not close to each other
						*core_data_pkg->near_interface_id_addrs_[i][j] = 1;
					}
				}
			}
	}
	//=============================================================================================//
	void LevelSet::writeMeshToPltFile(ofstream& output_file)
	{
		Vecu number_of_operation = total_number_of_data_points_;

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "phi, " << "n_x, " << "n_y " << "near_interface_id " << "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = DataPositionFromGlobalIndex(Vecu(i, j));
				output_file << data_position[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = DataPositionFromGlobalIndex(Vecu(i, j));
				output_file << data_position[1]<< " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << DataValueFromGlobalIndex<Real, LevelSetDataPackage::PackageData<Real>, 
					&LevelSetDataPackage::phi_>(Vecu(i, j)) << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << DataValueFromGlobalIndex<Vecd, LevelSetDataPackage::PackageData<Vecd>, 
					&LevelSetDataPackage::n_>(Vecu(i, j))[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << DataValueFromGlobalIndex<Vecd, LevelSetDataPackage::PackageData<Vecd>, 
					&LevelSetDataPackage::n_>(Vecu(i, j))[1] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << DataValueFromGlobalIndex<int, LevelSetDataPackage::PackageData<int>,
					&LevelSetDataPackage::near_interface_id_>(Vecu(i, j)) << " ";
			}
			output_file << " \n";
		}

	}
	//=============================================================================================//
	Vecd LevelSet::computeKernelIntegral(Vecd position, Kernel* kernel)
	{
		Vecd integral(0.0);
		Real delta = data_spacing_;
		Real phi = probeLevelSet(position);
		Real cutoff_radius = kernel->GetCutOffRadius();
		Real threshold = -cutoff_radius - delta;
		if (phi > threshold)
		{
			Vecu global_index_ = DataGlobalIndexFromPosition(position);
			for (int i = -2; i != 4; ++i)
				for (int j = -2; j != 4; ++j)
				{
					Vecu neighbor_index = Vecu(global_index_[0] + i, global_index_[1] + j);
					Real phi_neighbor = DataValueFromGlobalIndex<Real, LevelSetDataPackage::PackageData<Real>,
						&LevelSetDataPackage::phi_>(neighbor_index);
					if (phi_neighbor > -delta) {
						Vecd neighbor_position_ = DataPositionFromGlobalIndex(neighbor_index);
						Vecd displacement = position - neighbor_position_;
						if (displacement.norm() < cutoff_radius)
						integral += kernel->dW(displacement) 
								  * computeHeaviside(phi_neighbor, delta) * displacement.normalize();
					}
				}
		}
		return integral * data_spacing_ * data_spacing_;
	}
	//=============================================================================================//
}
//=============================================================================================//
