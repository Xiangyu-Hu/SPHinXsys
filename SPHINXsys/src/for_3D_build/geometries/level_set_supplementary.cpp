/**
 * @file 	level_set_supplementary.cpp
 * @author	Luhui Han, Chi Zhang Yongchuan YU and Xiangyu Hu
 */

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
		for (int i = 0; i != pkg_size; ++i)
			for (int j = 0; j != pkg_size; ++j)
				for (int k = 0; k != pkg_size; ++k)
				{
					phi_[i][j][k] = far_field_level_set;
					phi_gradient_[i][j][k] = Vecd(1.0);
					near_interface_id_[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::initializeBasicData(Shape &shape)
	{
		for (int i = 0; i != pkg_size; ++i)
			for (int j = 0; j != pkg_size; ++j)
				for (int k = 0; k != pkg_size; ++k)
				{
					Vec3d position = DataLowerBound() + Vec3d(i, j, k) * grid_spacing_;
					phi_[i][j][k] = shape.findSignedDistance(position);
					near_interface_id_[i][j][k] = phi_[i][j][k] < 0.0 ? -2 : 2;
				}
	}
	//=================================================================================================//
	void LevelSetDataPackage::stepReinitialization()
	{
		for_each3d<pkg_addrs_buffer, pkg_ops_end>(
			[&](int i, int j, int k)
			{
				// only reinitialize non cut cells
				if (*near_interface_id_addrs_[i][j][k] != 0)
				{
					Real phi_0 = *phi_addrs_[i][j][k];
					Real sign = phi_0 / sqrt(phi_0 * phi_0 + grid_spacing_ * grid_spacing_);
					Real dv_x = upwindDifference(sign, *phi_addrs_[i + 1][j][k] - phi_0, phi_0 - *phi_addrs_[i - 1][j][k]);
					Real dv_y = upwindDifference(sign, *phi_addrs_[i][j + 1][k] - phi_0, phi_0 - *phi_addrs_[i][j - 1][k]);
					Real dv_z = upwindDifference(sign, *phi_addrs_[i][j][k + 1] - phi_0, phi_0 - *phi_addrs_[i][j][k - 1]);
					*phi_addrs_[i][j][k] -= 0.3 * sign * (Vec3d(dv_x, dv_y, dv_z).norm() - grid_spacing_);
				}
			});
	}
	//=================================================================================================//
	void LevelSetDataPackage::stepDiffusionLevelSetSign()
	{
		for_each3d<pkg_addrs_buffer, pkg_ops_end>(
			[&](int i, int j, int k)
			{
				// near interface cells are not considered
				if (abs(*near_interface_id_addrs_[i][j][k]) > 1)
				{
					find_if3d<-1, 2>(
						[&](int l, int m, int n) -> bool
						{
							int near_interface_id = *near_interface_id_addrs_[i + l][j + m][k + n];
							bool is_found = abs(near_interface_id) == 1;
							if (is_found)
							{
								Real phi_0 = *phi_addrs_[i][j][k];
								*near_interface_id_addrs_[i][j][k] = near_interface_id;
								*phi_addrs_[i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
							}
							return is_found;
						});
				}
			});
	}
	//=================================================================================================//
	void LevelSetDataPackage::markNearInterface(Real small_shift_factor)
	{
		Real small_shift = small_shift_factor * grid_spacing_;
		// corner averages, note that the first row and first column are not used
		PackageTemporaryData<Real> corner_averages;
		for_each3d<1, pkg_addrs_size>(
			[&](int i, int j, int k)
			{
				corner_averages[i][j][k] = CornerAverage(phi_addrs_, Veci(i, j, k), Veci(-1, -1, -1));
			});

		for_each3d<pkg_addrs_buffer, pkg_ops_end>(
			[&](int i, int j, int k)
			{
				// first assume far cells
				Real phi_0 = *phi_addrs_[i][j][k];
				int near_interface_id = phi_0 > 0.0 ? 2 : -2;
				if (fabs(phi_0) < small_shift)
				{
					near_interface_id = 0;
					Real phi_average_0 = corner_averages[i][j][k];
					// find inner and outer cut cells
					for_each3d<0, 2>(
						[&](int l, int m, int n)
						{
							Real phi_average = corner_averages[i + l][j + m][k + n];
							if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
								near_interface_id = 1;
							if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
								near_interface_id = -1;
						});
					// find zero cut cells
					for_each3d<0, 2>(
						[&](int l, int m, int n)
						{
							Real phi_average = corner_averages[i + l][j + m][k + n];
							if (phi_average_0 * phi_average < 0.0)
								near_interface_id = 0;
						});
				}
				// assign this is to package
				*near_interface_id_addrs_[i][j][k] = near_interface_id;
			});
	}
	//=================================================================================================//
	LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing,
					   Shape &shape, SPHAdaptation &sph_adaptation)
		: LevelSet(tentative_bounds, data_spacing, 4, shape, sph_adaptation)
	{
		mesh_parallel_for(MeshRange(Vecu(0), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  initializeDataInACell(Vecu(i, j, k));
						  });

		finishDataPackages();
	}
	//=================================================================================================//
	void LevelSet::initializeSingularData(LevelSetDataPackage &data_pkg, Real far_field_level_set)
	{
		auto &kernel_weight = data_pkg.getPackageData(kernel_weight_);
		auto &kernel_gradient = data_pkg.getPackageData(kernel_gradient_);

		for (int i = 0; i != data_pkg.pkg_size; ++i)
			for (int j = 0; j != data_pkg.pkg_size; ++j)
				for (int k = 0; k != data_pkg.pkg_size; ++k)
				{
					data_pkg.phi_[i][j][k] = far_field_level_set;
					data_pkg.phi_gradient_[i][j][k] = Vecd(1.0);
					data_pkg.near_interface_id_[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
					kernel_weight[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
					kernel_gradient[i][j][k] = Vec3d(0);
				}
	}
	//=================================================================================================//
	void LevelSet::finishDataPackages()
	{
		mesh_parallel_for(MeshRange(Vecu(0), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  tagACellIsInnerPackage(Vecu(i, j, k));
						  });

		mesh_parallel_for(MeshRange(Vecu(0), number_of_cells_),
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
		return data_pkg_addrs_[cell_index[0]][cell_index[1]][cell_index[2]]->isCorePackage();
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
					if (data_pkg_addrs_[l][m][n]->isCorePackage())
						is_inner_pkg = true;
		return is_inner_pkg;
	}
	//=================================================================================================//
	void LevelSet::redistanceInterfaceForAPackage(LevelSetDataPackage *core_data_pkg)
	{
		int l = (int)core_data_pkg->CellIndexOnMesh()[0];
		int m = (int)core_data_pkg->CellIndexOnMesh()[1];
		int n = (int)core_data_pkg->CellIndexOnMesh()[2];

		for (int i = pkg_addrs_buffer; i != pkg_ops_end; ++i)
			for (int j = pkg_addrs_buffer; j != pkg_ops_end; ++j)
				for (int k = pkg_addrs_buffer; k != pkg_ops_end; ++k)
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
							// this immediate switch of near interface id
							// does not intervening with the identification of unresolved interface
							// based on the assumption that positive false_and negative bands are not close to each other
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
							// this immediate switch of near interface id
							// does not intervening with the identification of unresolved interface
							// based on the assumption that positive false_and negative bands are not close to each other
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

		Real integral(0.0);
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

		Vecd integral(0.0);
		if (fabs(phi) < threshold)
		{
			Vecu global_index_ = global_mesh_.CellIndexFromPosition(position);
			for (int i = -3; i != 4; ++i)
				for (int j = -3; j != 4; ++j)
					for (int k = -3; k != 4; ++k)
					{
						Vecu neighbor_index = Vecu(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
						Real phi_neighbor = DataValueFromGlobalIndex<
							Real, LevelSetDataPackage::PackageData<Real>, &LevelSetDataPackage::phi_>(neighbor_index);
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
		mesh_parallel_for(MeshRange(Vecu(0), number_of_cells_),
						  [&](size_t i, size_t j, size_t k)
						  {
							  initializeDataInACellFromCoarse(Vecu(i, j, k));
						  });

		finishDataPackages();
	}
	//=============================================================================================//
}
