#include "level_set.h"
#include "mesh_with_data_packages.hpp"
#include "mesh_iterators.hpp"
#include "base_kernel.h"
#include "base_particles.h"
#include "base_particle_dynamics.h"
#include "base_body.h"

namespace SPH
{
	//=================================================================================================//
	void LevelSetDataPackage::initializeSingularData(Real far_field_level_set)
	{
		for (int i = 0; i != PackageSize(); ++i)
			for (int j = 0; j != PackageSize(); ++j)
				for (int k = 0; k != PackageSize(); ++k)
				{
					phi_[i][j][k] = far_field_level_set;
					phi_gradient_[i][j][k] = Vecd::Ones();
					kernel_weight_[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
					kernel_gradient_[i][j][k] = Vecd::Zero();
					near_interface_id_[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::initializeBasicData(Shape &shape)
	{
		for (int i = 0; i != PackageSize(); ++i)
			for (int j = 0; j != PackageSize(); ++j)
				for (int k = 0; k != PackageSize(); ++k)
				{
					Vecd position = DataLowerBound() + Vecd(i, j, k) * grid_spacing_;
					phi_[i][j][k] = shape.findSignedDistance(position);
					near_interface_id_[i][j][k] = phi_[i][j][k] < 0.0 ? -2 : 2;
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::computeKernelIntegrals(LevelSet &level_set)
	{
		for (int i = 0; i != PackageSize(); ++i)
			for (int j = 0; j != PackageSize(); ++j)
				for (int k = 0; k != PackageSize(); ++k)
				{
					Vecd position = DataLowerBound() + Vecd(i, j, k) * grid_spacing_;
					kernel_weight_[i][j][k] = level_set.computeKernelIntegral(position);
					kernel_gradient_[i][j][k] = level_set.computeKernelGradientIntegral(position);
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::stepReinitialization()
	{
		for (int i = AddressBufferWidth(); i != OperationUpperBound(); ++i)
			for (int j = AddressBufferWidth(); j != OperationUpperBound(); ++j)
				for (int k = AddressBufferWidth(); k != OperationUpperBound(); ++k)
				{
					// only reinitialize non cut cells
					if (*near_interface_id_addrs_[i][j][k] != 0)
					{
						Real phi_0 = *phi_addrs_[i][j][k];
						Real s = phi_0 / sqrt(phi_0 * phi_0 + grid_spacing_ * grid_spacing_);
						// x direction
						Real dv_xp = (*phi_addrs_[i + 1][j][k] - phi_0);
						Real dv_xn = (phi_0 - *phi_addrs_[i - 1][j][k]);
						Real dv_x = dv_xp;
						if (s * dv_xp >= 0.0 && s * dv_xn >= 0.0)
							dv_x = dv_xn;
						if (s * dv_xp <= 0.0 && s * dv_xn <= 0.0)
							dv_x = dv_xp;
						if (s * dv_xp > 0.0 && s * dv_xn < 0.0)
							dv_x = 0.0;
						if (s * dv_xp < 0.0 && s * dv_xn > 0.0)
						{
							Real ss = s * (fabs(dv_xp) - fabs(dv_xn)) / (dv_xp - dv_xn);
							if (ss > 0.0)
								dv_x = dv_xn;
						}
						// y direction
						Real dv_yp = (*phi_addrs_[i][j + 1][k] - phi_0);
						Real dv_yn = (phi_0 - *phi_addrs_[i][j - 1][k]);
						Real dv_y = dv_yp;
						if (s * dv_yp >= 0.0 && s * dv_yn >= 0.0)
							dv_y = dv_yn;
						if (s * dv_yp <= 0.0 && s * dv_yn <= 0.0)
							dv_y = dv_yp;
						if (s * dv_yp > 0.0 && s * dv_yn < 0.0)
							dv_y = 0.0;
						if (s * dv_yp < 0.0 && s * dv_yn > 0.0)
						{
							Real ss = s * (fabs(dv_yp) - fabs(dv_yn)) / (dv_yp - dv_yn);
							if (ss > 0.0)
								dv_y = dv_yn;
						}
						// z direction
						Real dv_zp = (*phi_addrs_[i][j][k + 1] - phi_0);
						Real dv_zn = (phi_0 - *phi_addrs_[i][j][k - 1]);
						Real dv_z = dv_zp;
						if (s * dv_zp >= 0.0 && s * dv_zn >= 0.0)
							dv_z = dv_zn;
						if (s * dv_zp <= 0.0 && s * dv_zn <= 0.0)
							dv_z = dv_zp;
						if (s * dv_zp > 0.0 && s * dv_zn < 0.0)
							dv_z = 0.0;
						if (s * dv_zp < 0.0 && s * dv_zn > 0.0)
						{
							Real ss = s * (fabs(dv_zp) - fabs(dv_zn)) / (dv_zp - dv_zn);
							if (ss > 0.0)
								dv_z = dv_zn;
						}
						// time stepping
						*phi_addrs_[i][j][k] -=
							0.3 * s * (sqrt(dv_x * dv_x + dv_y * dv_y + dv_z * dv_z) - grid_spacing_);
					}
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::stepDiffusionLevelSetSign()
	{
		for (int i = AddressBufferWidth(); i != OperationUpperBound(); ++i)
			for (int j = AddressBufferWidth(); j != OperationUpperBound(); ++j)
				for (int k = AddressBufferWidth(); k != OperationUpperBound(); ++k)
				{
					// near interface cells are not considered
					if (abs(*near_interface_id_addrs_[i][j][k]) > 1)
					{
						Real phi_0 = *phi_addrs_[i][j][k];
						for (int l = -1; l != 2; ++l)
							for (int m = -1; m != 2; ++m)
								for (int n = -1; n != 2; ++n)
								{
									int index_x = i + l;
									int index_y = j + m;
									int index_z = k + n;
									int near_interface_id = *near_interface_id_addrs_[index_x][index_y][index_z];
									if (abs(near_interface_id) == 1)
									{
										*near_interface_id_addrs_[i][j][k] = near_interface_id;
										*phi_addrs_[i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
										break;
									}
								}
					}
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::markNearInterface(Real small_shift_factor)
	{
		Real small_shift = small_shift_factor * grid_spacing_;
		// corner averages, note that the first row and first column are not used
		PackageTemporaryData<Real> corner_averages;
		for (int i = 1; i != AddressSize(); ++i)
			for (int j = 1; j != AddressSize(); ++j)
				for (int k = 1; k != AddressSize(); ++k)
				{
					corner_averages[i][j][k] = CornerAverage(phi_addrs_, Veci(i, j, k), Veci(-1, -1, -1));
				}

		for (int i = AddressBufferWidth(); i != OperationUpperBound(); ++i)
			for (int j = AddressBufferWidth(); j != OperationUpperBound(); ++j)
				for (int k = AddressBufferWidth(); k != OperationUpperBound(); ++k)
				{
					// first assume far cells
					Real phi_0 = *phi_addrs_[i][j][k];
					int near_interface_id = phi_0 > 0.0 ? 2 : -2;

					if (fabs(phi_0) < small_shift)
					{
						Real phi_average_0 = corner_averages[i][j][k];
						// find inner and outer cut cells
						for (int l = 0; l != 2; ++l)
							for (int m = 0; m != 2; ++m)
								for (int n = 0; n != 2; ++n)
								{
									int index_x = i + l;
									int index_y = j + m;
									int index_z = k + n;
									Real phi_average = corner_averages[index_x][index_y][index_z];
									if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
										near_interface_id = 1;
									if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
										near_interface_id = -1;
								}
						// find zero cut cells
						for (int l = 0; l != 2; ++l)
							for (int m = 0; m != 2; ++m)
								for (int n = 0; n != 2; ++n)
								{
									int index_x = i + l;
									int index_y = j + m;
									int index_z = k + n;
									Real phi_average = corner_averages[index_x][index_y][index_z];
									if (phi_average_0 * phi_average < 0.0)
										near_interface_id = 0;
								}
						// find cells between cut cells
						if (fabs(phi_0) < small_shift && abs(near_interface_id) != 1)
							near_interface_id = 0;
					}

					// assign this is to package
					*near_interface_id_addrs_[i][j][k] = near_interface_id;
				}
	}
	//=================================================================================================//
	LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing,
					   Shape &shape, SPHAdaptation &sph_adaptation)
		: LevelSet(tentative_bounds, data_spacing, 4, shape, sph_adaptation)
	{
		mesh_parallel_for(MeshRange(Vecu::Zero(), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  initializeDataInACell(Vecu(i, j, k));
						  });

		finishDataPackages();
	}
	//=================================================================================================//
	void LevelSet::finishDataPackages()
	{
		mesh_parallel_for(MeshRange(Vecu::Zero(), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  tagACellIsInnerPackage(Vecu(i, j, k));
						  });

		mesh_parallel_for(MeshRange(Vecu::Zero(), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  initializePackageAddressesInACell(Vecu(i, j, k));
						  });

		updateLevelSetGradient();
		updateKernelIntegrals();
	}
	//=================================================================================================//
	bool LevelSet::isWithinCorePackage(Vecd position)
	{
		Vecu cell_index = CellIndexFromPosition(position);
		return data_pkg_addrs_[cell_index[0]][cell_index[1]][cell_index[2]]->is_core_pkg_;
	}
	//=============================================================================================//
	bool LevelSet::isInnerPackage(const Vecu &cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

		bool is_inner_pkg = false;
		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++n)
					if (data_pkg_addrs_[l][m][n]->is_core_pkg_)
						is_inner_pkg = true;
		return is_inner_pkg;
	}
	//=================================================================================================//
	void LevelSet::redistanceInterfaceForAPackage(LevelSetDataPackage *core_data_pkg)
	{
		int l = (int)core_data_pkg->pkg_index_[0];
		int m = (int)core_data_pkg->pkg_index_[1];
		int n = (int)core_data_pkg->pkg_index_[2];

		for (int i = pkg_addrs_buffer_; i != pkg_operations_; ++i)
			for (int j = pkg_addrs_buffer_; j != pkg_operations_; ++j)
				for (int k = pkg_addrs_buffer_; k != pkg_operations_; ++k)
				{
					int near_interface_id = *core_data_pkg->near_interface_id_addrs_[i][j][k];
					if (near_interface_id == 0)
					{
						bool positive_band = false;
						bool negative_band = false;
						for (int r = -1; r < 2; ++r)
							for (int s = -1; s < 2; ++s)
								for (int t = -1; t < 2; ++t)
								{
									int neighbor_near_interface_id =
										*core_data_pkg->near_interface_id_addrs_[i + r][j + s][k + t];
									if (neighbor_near_interface_id >= 1)
										positive_band = true;
									if (neighbor_near_interface_id <= -1)
										negative_band = true;
								}
						if (positive_band == false)
						{
							Real min_distance_p = 5.0 * data_spacing_;
							for (int x = -4; x != 5; ++x)
								for (int y = -4; y != 5; ++y)
									for (int z = -4; z != 5; ++z)
									{
										std::pair<int, int> x_pair = CellShiftAndDataIndex(i + x);
										std::pair<int, int> y_pair = CellShiftAndDataIndex(j + y);
										std::pair<int, int> z_pair = CellShiftAndDataIndex(k + z);
										LevelSetDataPackage *neighbor_pkg = data_pkg_addrs_[l + x_pair.first][m + y_pair.first][n + z_pair.first];
										int neighbor_near_interface_id = neighbor_pkg->near_interface_id_[x_pair.second][y_pair.second][z_pair.second];
										if (neighbor_near_interface_id >= 1)
										{
											Real phi_p_ = neighbor_pkg->phi_[x_pair.second][y_pair.second][z_pair.second];
											Vecd norm_to_face = neighbor_pkg->phi_gradient_[x_pair.second][y_pair.second][z_pair.second];
											norm_to_face /= norm_to_face.norm() + TinyReal;
											min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
										}
									}
							*core_data_pkg->phi_addrs_[i][j][k] = -min_distance_p;
							/** This immediate switch of near interface id
							 * does not intervening with the identification of unresolved interface
							 * based on the assumption that positive false_and negative bands are not close to each other.
							  */
							*core_data_pkg->near_interface_id_addrs_[i][j][k] = -1;
						}
						if (negative_band == false)
						{
							Real min_distance_n = 5.0 * data_spacing_;
							for (int x = -4; x != 5; ++x)
								for (int y = -4; y != 5; ++y)
									for (int z = -4; z != 5; ++z)
									{
										std::pair<int, int> x_pair = CellShiftAndDataIndex(i + x);
										std::pair<int, int> y_pair = CellShiftAndDataIndex(j + y);
										std::pair<int, int> z_pair = CellShiftAndDataIndex(k + z);
										LevelSetDataPackage *neighbor_pkg = data_pkg_addrs_[l + x_pair.first][m + y_pair.first][n + z_pair.first];
										int neighbor_near_interface_id = neighbor_pkg->near_interface_id_[x_pair.second][y_pair.second][z_pair.second];
										if (neighbor_near_interface_id <= -1)
										{
											Real phi_n_ = neighbor_pkg->phi_[x_pair.second][y_pair.second][z_pair.second];
											Vecd norm_to_face = neighbor_pkg->phi_gradient_[x_pair.second][y_pair.second][z_pair.second];
											norm_to_face /= norm_to_face.norm() + TinyReal;
											min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
										}
									}
							*core_data_pkg->phi_addrs_[i][j][k] = min_distance_n;
							/** This immediate switch of near interface id
							 * does not intervening with the identification of unresolved interface
							 * based on the assumption that positive false_and negative bands are not close to each other. 
							 */
							*core_data_pkg->near_interface_id_addrs_[i][j][k] = 1;
						}
					}
				}
	}
	//=================================================================================================//
	void LevelSet::writeMeshFieldToPlt(std::ofstream &output_file)
	{
		Vecu number_of_operation = global_mesh_.NumberOfGridPoints();

		output_file << "\n";
		output_file << "title='View'"
					<< "\n";
		output_file << "variables= "
					<< "x, "
					<< "y, "
					<< "z, "
					<< "phi, "
					<< "n_x, "
					<< "n_y, "
					<< "n_z "
					<< "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
					<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = global_mesh_.GridPositionFromIndex(Vecu(i, j, k));
					output_file << data_position[0] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = global_mesh_.GridPositionFromIndex(Vecu(i, j, k));
					output_file << data_position[1] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = global_mesh_.GridPositionFromIndex(Vecu(i, j, k));
					output_file << data_position[2] << " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<Real, LevelSetDataPackage::PackageData<Real>,
															&LevelSetDataPackage::phi_>(Vecu(i, j, k))
								<< " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<Vecd, LevelSetDataPackage::PackageData<Vecd>,
															&LevelSetDataPackage::phi_gradient_>(Vecu(i, j, k))[0]
								<< " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<Vecd, LevelSetDataPackage::PackageData<Vecd>,
															&LevelSetDataPackage::phi_gradient_>(Vecu(i, j, k))[1]
								<< " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<Vecd, LevelSetDataPackage::PackageData<Vecd>,
															&LevelSetDataPackage::phi_gradient_>(Vecu(i, j, k))[2]
								<< " ";
				}
				output_file << " \n";
			}

		for (size_t k = 0; k != number_of_operation[2]; ++k)
			for (size_t j = 0; j != number_of_operation[1]; ++j)
			{
				for (size_t i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << DataValueFromGlobalIndex<int, LevelSetDataPackage::PackageData<int>,
															&LevelSetDataPackage::near_interface_id_>(Vecu(i, j, k))
								<< " ";
				}
				output_file << " \n";
			}
	}
	//=============================================================================================//
	Real LevelSet::computeKernelIntegral(const Vecd &position)
	{
		Real phi = probeSignedDistance(position);
		Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
		Real threshold = cutoff_radius + data_spacing_;

		Real integral(0);
		if (fabs(phi) < threshold)
		{
			Vecu global_index_ = global_mesh_.CellIndexFromPosition(position);
			for (int i = -3; i != 4; ++i)
				for (int j = -3; j != 4; ++j)
					for (int k = -3; k != 4; ++k)
					{
						Vecu neighbor_index = Vecu(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
						Real phi_neighbor = DataValueFromGlobalIndex<Real, LevelSetDataPackage::PackageData<Real>,
																	 &LevelSetDataPackage::phi_>(neighbor_index);
						if (phi_neighbor > -data_spacing_)
						{
							Vecd displacement = position - global_mesh_.GridPositionFromIndex(neighbor_index);
							Real distance = displacement.norm();
							if (distance < cutoff_radius)
								integral += kernel_.W(global_h_ratio_, distance, displacement) * computeHeaviside(phi_neighbor, data_spacing_);
						}
					}
		}
		return phi > threshold ? 1.0 : integral * data_spacing_ * data_spacing_ * data_spacing_;
	}
	//=============================================================================================//
	Vecd LevelSet::computeKernelGradientIntegral(const Vecd &position)
	{
		Real phi = probeSignedDistance(position);
		Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
		Real threshold = cutoff_radius + data_spacing_;

		Vecd integral = Vecd::Zero();
		if (fabs(phi) < threshold)
		{
			Vecu global_index_ = global_mesh_.CellIndexFromPosition(position);
			for (int i = -3; i != 4; ++i)
				for (int j = -3; j != 4; ++j)
					for (int k = -3; k != 4; ++k)
					{
						Vecu neighbor_index = Vecu(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
						Real phi_neighbor = DataValueFromGlobalIndex<Real, LevelSetDataPackage::PackageData<Real>,
																	 &LevelSetDataPackage::phi_>(neighbor_index);
						if (phi_neighbor > -data_spacing_)
						{
							Vecd displacement = position - global_mesh_.GridPositionFromIndex(neighbor_index);
							Real distance = displacement.norm();
							if (distance < cutoff_radius)
								integral += kernel_.dW(global_h_ratio_, distance, displacement) *
											computeHeaviside(phi_neighbor, data_spacing_) * displacement / (distance + TinyReal);
						}
					}
		}
		return integral * data_spacing_ * data_spacing_ * data_spacing_;
	}
	//=============================================================================================//
	RefinedLevelSet::RefinedLevelSet(BoundingBox tentative_bounds, LevelSet &coarse_level_set,
									 Shape &shape, SPHAdaptation &sph_adaptation)
		: RefinedMesh(tentative_bounds, coarse_level_set, 4, shape, sph_adaptation)
	{
		mesh_parallel_for(MeshRange(Vecu::Zero(), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  initializeDataInACellFromCoarse(Vecu(i, j, k));
						  });

		finishDataPackages();
	}
	//=============================================================================================//
}
