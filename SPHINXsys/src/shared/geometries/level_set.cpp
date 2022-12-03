#include "level_set.h"
#include "adaptation.h"

#include "mesh_with_data_packages.hpp"

namespace SPH
{
	//=================================================================================================//
	LevelSetDataPackage::
		LevelSetDataPackage() : GridDataPackage<4, 6>(), is_core_pkg_(false) {}
	//=================================================================================================//
	void LevelSetDataPackage::registerAllVariables()
	{
		registerPackageData(phi_, phi_addrs_);
		registerPackageData(phi_gradient_, phi_gradient_addrs_);
		registerPackageData(kernel_weight_, kernel_weight_addrs_);
		registerPackageData(kernel_gradient_, kernel_gradient_addrs_);
		registerPackageData(near_interface_id_, near_interface_id_addrs_);
	}
	//=================================================================================================//
	void LevelSetDataPackage::computeLevelSetGradient()
	{
		computeGradient(phi_addrs_, phi_gradient_addrs_);
	}
	//=================================================================================================//
	BaseLevelSet ::BaseLevelSet(Shape &shape, SPHAdaptation &sph_adaptation)
		: BaseMeshField("LevelSet")
		, shape_(shape)
		, sph_adaptation_(sph_adaptation)
	{
		if (!shape_.isValid())
		{
			std::cout << "\n BaseLevelSet Error: shape_ is invalid." << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
	}
	//=================================================================================================//
	Real BaseLevelSet::computeHeaviside(Real phi, Real half_width)
	{
		Real heaviside = 0.0;
		Real normalized_phi = phi / half_width;
		if (phi < half_width && phi > -half_width)
			heaviside = (0.5 + 0.5 * normalized_phi) + 0.5 * sin(Pi * normalized_phi) / Pi;
		if (normalized_phi > 1.0)
			heaviside = 1.0;
		return heaviside;
	}
	//=================================================================================================//
	LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, Shape &shape, SPHAdaptation &sph_adaptation)
		: MeshWithGridDataPackages<BaseLevelSet, LevelSetDataPackage>(tentative_bounds, data_spacing, buffer_size, shape, sph_adaptation)
		, global_h_ratio_(sph_adaptation.ReferenceSpacing() / data_spacing)
		, kernel_(*sph_adaptation.getKernel())
	{
		Real far_field_distance = grid_spacing_ * (Real)buffer_width_;
		initializeASingularDataPackage(-far_field_distance);
		initializeASingularDataPackage(far_field_distance);
	}
	//=================================================================================================//
	void LevelSet::initializeAddressesInACell(const Vecu &cell_index)
	{
		initializePackageAddressesInACell(cell_index);
	}
	//=================================================================================================//
	void LevelSet::updateLevelSetGradient()
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { inner_data_pkgs_[i]->computeLevelSetGradient(); });
	}
	//=================================================================================================//
	void LevelSet::updateKernelIntegrals()
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { inner_data_pkgs_[i]->computeKernelIntegrals(*this); });
	}
	//=================================================================================================//
	Vecd LevelSet::probeNormalDirection(const Vecd &position)
	{
		Vecd probed_value = probeLevelSetGradient(position);

		Real threshold = 1.0e-2 * data_spacing_;
		while (probed_value.norm() < threshold)
		{
			Vecd jittered = position; // jittering
			for (int l = 0; l != position.size(); ++l)
				jittered[l] += (((Real)rand() / (RAND_MAX)) - 0.5) * 0.5 * data_spacing_;
			probed_value = probeLevelSetGradient(jittered);
		}
		return probed_value.normalized();
	}
	//=================================================================================================//
	Vecd LevelSet::probeLevelSetGradient(const Vecd &position)
	{
		return probeMesh<Vecd, LevelSetDataPackage::PackageDataAddress<Vecd>,
						 &LevelSetDataPackage::phi_gradient_addrs_>(position);
	}
	//=================================================================================================//
	Real LevelSet::probeSignedDistance(const Vecd &position)
	{
		return probeMesh<Real, LevelSetDataPackage::PackageDataAddress<Real>,
						 &LevelSetDataPackage::phi_addrs_>(position);
	}
	//=================================================================================================//
	Real LevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
	{
		return probeMesh<Real, LevelSetDataPackage::PackageDataAddress<Real>,
						 &LevelSetDataPackage::kernel_weight_addrs_>(position);
	}
	//=================================================================================================//
	Vecd LevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
	{
		return probeMesh<Vecd, LevelSetDataPackage::PackageDataAddress<Vecd>,
						 &LevelSetDataPackage::kernel_gradient_addrs_>(position);
	}
	//=============================================================================================//
	void LevelSet::reinitializeLevelSet()
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { inner_data_pkgs_[i]->stepReinitialization(); });
	}
	//=================================================================================================//
	void LevelSet::markNearInterface(Real small_shift_factor)
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { inner_data_pkgs_[i]->markNearInterface(small_shift_factor); });
	}
	//=================================================================================================//
	void LevelSet::redistanceInterface()
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { redistanceInterfaceForAPackage(inner_data_pkgs_[i]); });
	}
	//=================================================================================================//
	void LevelSet::diffuseLevelSetSign()
	{
		package_parallel_for(inner_data_pkgs_, [&](size_t i)
							 { inner_data_pkgs_[i]->stepDiffusionLevelSetSign(); });
	}
	//=================================================================================================//
	void LevelSet::cleanInterface(Real small_shift_factor)
	{
		markNearInterface(small_shift_factor);
		redistanceInterface();
		reinitializeLevelSet();
		updateLevelSetGradient();
		updateKernelIntegrals();
	}
	//=============================================================================================//
	void LevelSet::correctTopology(Real small_shift_factor)
	{
		markNearInterface(small_shift_factor);
		for (size_t i = 0; i != 10; ++i)
			diffuseLevelSetSign();
		updateLevelSetGradient();
		updateKernelIntegrals();
	}
	//=================================================================================================//
	bool LevelSet::probeIsWithinMeshBound(const Vecd &position)
	{
		bool is_bounded = true;
		Vecu cell_pos = CellIndexFromPosition(position);
		for (int i = 0; i != position.size(); ++i)
		{
			if (cell_pos[i] < 2)
				is_bounded = false;
			if (cell_pos[i] > (number_of_cells_[i] - 2))
				is_bounded = false;
		}
		return is_bounded;
	}
	//=================================================================================================//
	LevelSetDataPackage* LevelSet::createDataPackage(const Vecu &cell_index, const Vecd &cell_position)
	{
		mutex_my_pool.lock();
		LevelSetDataPackage* new_data_pkg = data_pkg_pool_.malloc();
		mutex_my_pool.unlock();

		new_data_pkg->registerAllVariables();
		Vecd pkg_lower_bound = GridPositionFromCellPosition(cell_position);
		new_data_pkg->initializePackageGeometry(pkg_lower_bound, data_spacing_);
		new_data_pkg->initializeBasicData(shape_);
		new_data_pkg->pkg_index_ = cell_index;
		assignDataPackageAddress(cell_index, new_data_pkg);
		return new_data_pkg;
	}
	//=================================================================================================//
	void LevelSet::initializeDataInACell(const Vecu &cell_index)
	{
		Vecd cell_position = CellPositionFromIndex(cell_index);
		Real signed_distance = shape_.findSignedDistance(cell_position);
		Vecd normal_direction = shape_.findNormalDirection(cell_position);
		Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
		
		if (measure < grid_spacing_)
		{
			LevelSetDataPackage *new_data_pkg = createDataPackage(cell_index, cell_position);
			new_data_pkg->is_core_pkg_ = true;
			core_data_pkgs_.push_back(new_data_pkg);
		}
		else
		{
			LevelSetDataPackage *singular_data_pkg =
				shape_.checkContain(cell_position) ? singular_data_pkgs_addrs_[0] : singular_data_pkgs_addrs_[1];
			assignDataPackageAddress(cell_index, singular_data_pkg);
		}
	}
	//=============================================================================================//
	void LevelSet::tagACellIsInnerPackage(const Vecu &cell_index)
	{
		bool is_inner_pkg = isInnerPackage(cell_index);

		if (is_inner_pkg)
		{
			LevelSetDataPackage *current_data_pkg = DataPackageFromCellIndex(cell_index);
			if (current_data_pkg->is_core_pkg_)
			{
				current_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(current_data_pkg);
			}
			else
			{
				Vecd cell_position = CellPositionFromIndex(cell_index);
				LevelSetDataPackage *new_data_pkg = createDataPackage(cell_index, cell_position);
				new_data_pkg->is_inner_pkg_ = true;
				inner_data_pkgs_.push_back(new_data_pkg);
			}
		}
	}
	//=============================================================================================//
	void RefinedLevelSet::initializeDataInACellFromCoarse(const Vecu &cell_index)
	{
		Vecd cell_position = CellPositionFromIndex(cell_index);
		LevelSetDataPackage *singular_data_pkg = coarse_mesh_.probeSignedDistance(cell_position) < 0.0
													 ? singular_data_pkgs_addrs_[0]
													 : singular_data_pkgs_addrs_[1];
		assignDataPackageAddress(cell_index, singular_data_pkg);
		if (coarse_mesh_.isWithinCorePackage(cell_position))
		{
			Real signed_distance = shape_.findSignedDistance(cell_position);
			Vecd normal_direction = shape_.findNormalDirection(cell_position);
			Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
			if (measure < grid_spacing_)
			{
				LevelSetDataPackage *new_data_pkg = createDataPackage(cell_index, cell_position);
				new_data_pkg->is_core_pkg_ = true;
				core_data_pkgs_.push_back(new_data_pkg);
			}
		}
	}
	//=============================================================================================//
	MultilevelLevelSet::MultilevelLevelSet(BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation)
		: MultilevelMesh<BaseLevelSet, LevelSet, RefinedLevelSet>(tentative_bounds, reference_data_spacing, total_levels, shape, sph_adaptation) 
	{}
	//=================================================================================================//
	size_t MultilevelLevelSet::getCoarseLevel(Real h_ratio)
	{
		for (size_t level = total_levels_; level != 0; --level)
			if (h_ratio > mesh_levels_[level - 1]->global_h_ratio_)
				return level - 1; // jump out the loop!

		std::cout << "\n Error: LevelSet level searching out of bound!" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
		return 999; // means an error in level searching
	};
	//=================================================================================================//
	void MultilevelLevelSet::cleanInterface(Real small_shift_factor)
	{
		mesh_levels_.back()->cleanInterface(small_shift_factor);
	}
	//=============================================================================================//
	void MultilevelLevelSet::correctTopology(Real small_shift_factor)
	{
		mesh_levels_.back()->correctTopology(small_shift_factor);
	}
	//=============================================================================================//
	Real MultilevelLevelSet::probeSignedDistance(const Vecd &position)
	{
		return mesh_levels_[getProbeLevel(position)]->probeSignedDistance(position);
	}
	//=============================================================================================//
	Vecd MultilevelLevelSet::probeNormalDirection(const Vecd &position)
	{
		return mesh_levels_[getProbeLevel(position)]->probeNormalDirection(position);
	}
	//=============================================================================================//
	Vecd MultilevelLevelSet::probeLevelSetGradient(const Vecd &position)
	{
		return mesh_levels_[getProbeLevel(position)]->probeLevelSetGradient(position);
	}
	//=============================================================================================//
	size_t MultilevelLevelSet::getProbeLevel(const Vecd &position)
	{
		for (size_t level = total_levels_; level != 0; --level)
			if (mesh_levels_[level - 1]->isWithinCorePackage(position))
				return level - 1; // jump out of the loop!
		return 0;
	}
	//=================================================================================================//
	Real MultilevelLevelSet::probeKernelIntegral(const Vecd &position, Real h_ratio)
	{
		size_t coarse_level = getCoarseLevel(h_ratio);
		Real alpha = (mesh_levels_[coarse_level + 1]->global_h_ratio_ - h_ratio) /
					 (mesh_levels_[coarse_level + 1]->global_h_ratio_ - mesh_levels_[coarse_level]->global_h_ratio_);
		Real coarse_level_value = mesh_levels_[coarse_level]->probeKernelIntegral(position);
		Real fine_level_value = mesh_levels_[coarse_level + 1]->probeKernelIntegral(position);

		return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
	}
	//=================================================================================================//
	Vecd MultilevelLevelSet::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
	{
		size_t coarse_level = getCoarseLevel(h_ratio);
		Real alpha = (mesh_levels_[coarse_level + 1]->global_h_ratio_ - h_ratio) /
					 (mesh_levels_[coarse_level + 1]->global_h_ratio_ - mesh_levels_[coarse_level]->global_h_ratio_);
		Vecd coarse_level_value = mesh_levels_[coarse_level]->probeKernelGradientIntegral(position);
		Vecd fine_level_value = mesh_levels_[coarse_level + 1]->probeKernelGradientIntegral(position);

		return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
	}
	//=================================================================================================//
	bool MultilevelLevelSet::probeIsWithinMeshBound(const Vecd &position)
	{
		bool is_bounded = true;
		for (size_t l = 0; l != total_levels_; ++l)
		{
			if (!mesh_levels_[l]->probeIsWithinMeshBound(position))
			{
				is_bounded = false;
				break;
			};
		}
		return is_bounded;
	}
	//=============================================================================================//
}
