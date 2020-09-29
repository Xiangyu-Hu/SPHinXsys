/**
 * @file 	level_set.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "level_set.h"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	LevelSetDataPackage::
		LevelSetDataPackage() : BaseDataPackage<4, 6>(), is_core_pkg_(false)
	{
		initializePackageDataAddress(phi_, phi_addrs_);
		initializePackageDataAddress(n_, n_addrs_);
		initializePackageDataAddress(kappa_, kappa_addrs_);
		initializePackageDataAddress(near_interface_id_, near_interface_id_addrs_);
	}
	//=================================================================================================//
	void LevelSetDataPackage::
		assignAllPackageDataAddress(Vecu addrs_index, LevelSetDataPackage* src_pkg, Vecu data_index)
	{
		assignPackageDataAddress(phi_addrs_, addrs_index, src_pkg->phi_, data_index);
		assignPackageDataAddress(n_addrs_, addrs_index, src_pkg->n_, data_index);
		assignPackageDataAddress(kappa_addrs_, addrs_index, src_pkg->kappa_, data_index);
		assignPackageDataAddress(near_interface_id_addrs_, addrs_index, src_pkg->near_interface_id_, data_index);
	}
	//=================================================================================================//
	void LevelSetDataPackage::computeNormalDirection()
	{
		computeNormalizedGradient(phi_addrs_, n_addrs_);
	}
	//=================================================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd lower_bound,
			Vecd upper_bound, Real grid_spacing, size_t buffer_width)
		: BaseMeshWithDataPackages(lower_bound, upper_bound,
			grid_spacing, buffer_width) {}
	//=================================================================================================//
	BaseLevelSet
		::BaseLevelSet(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: BaseMeshWithDataPackages(mesh_lower_bound, number_of_cells,
			cell_spacing) {}
	//=================================================================================================//
	Real BaseLevelSet::computeHeaviside(Real phi, Real half_width)
	{
		Real heaviside = 0.0;
		Real normalized_phi = phi / half_width;
		if (phi < half_width && phi > - half_width)
			heaviside = (0.5 + 0.5 * normalized_phi) + 0.5 * sin(Pi * normalized_phi) / Pi;
		if (normalized_phi > 1.0) heaviside = 1.0;
		return heaviside;
	} 
	//=================================================================================================//
	LevelSet
		::LevelSet(ComplexShape& complex_shape, 
			Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_width)
		: MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>(lower_bound,
			upper_bound, grid_spacing, buffer_width), complex_shape_(complex_shape)
	{
		Real far_field_distance = cell_spacing_ * (Real)buffer_width * 2.0;
		LevelSetDataPackage* negative_far_field = new LevelSetDataPackage();
		negative_far_field->initializeWithUniformData(-far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(negative_far_field);
		LevelSetDataPackage* positive_far_field = new LevelSetDataPackage();
		positive_far_field->initializeWithUniformData(far_field_distance, Vecd(0));
		singular_data_pkgs_addrs.push_back(positive_far_field);
		initializeDataPackages();
	}
	//=================================================================================================//
	void LevelSet::initializeDataPackages()
	{
		MeshFunctor initialize_data_in_a_cell = std::bind(&LevelSet::initializeDataInACell, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, initialize_data_in_a_cell);
		MeshFunctor tag_a_cell_inner_pkg = std::bind(&LevelSet::tagACellIsInnerPackage, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, tag_a_cell_inner_pkg);
		MeshFunctor initial_address_in_a_cell = std::bind(&LevelSet::initializeAddressesInACell, this, _1, _2);
		MeshIterator_parallel(Vecu(0), number_of_cells_, initial_address_in_a_cell);
		updateNormalDirection();
	}
	//=================================================================================================//
	void LevelSet::initializeAddressesInACell(Vecu cell_index, Real dt)
	{
		initializePackageAddressesInACell(cell_index);
	}
	//=================================================================================================//
	void LevelSet::updateNormalDirection()
	{
		PackageFunctor<void, LevelSetDataPackage> update_normal_diraction
			= std::bind(&LevelSet::updateNormalDirectionForAPackage, this, _1, _2);
		PackageIterator_parallel<LevelSetDataPackage>(inner_data_pkgs_, update_normal_diraction);
	}
	//=================================================================================================//
	Vecd LevelSet::probeNormalDirection(Vecd position)
	{
		return probeMesh<Vecd, LevelSetDataPackage::PackageDataAddress<Vecd>, &LevelSetDataPackage::n_addrs_>(position);
	}
	//=================================================================================================//
	Real LevelSet::probeLevelSet(Vecd position)
	{
		return probeMesh<Real, LevelSetDataPackage::PackageDataAddress<Real>, &LevelSetDataPackage::phi_addrs_>(position);
		
	}
	//=================================================================================================//
	void LevelSet::
		updateNormalDirectionForAPackage(LevelSetDataPackage* inner_data_pkg, Real dt)
	{
		inner_data_pkg->computeNormalDirection();
	}
	//=================================================================================================//
	void LevelSet::
		stepReinitializationForAPackage(LevelSetDataPackage* inner_data_pkg, Real dt)
	{
		inner_data_pkg->stepReinitialization();
	}
	//=============================================================================================//
	void LevelSet::reinitializeLevelSet()
	{
		PackageFunctor<void, LevelSetDataPackage> reinitialize_levelset
			= std::bind(&LevelSet::stepReinitializationForAPackage, this, _1, _2);
		for (size_t i = 0; i < 50; ++i)
			PackageIterator_parallel<LevelSetDataPackage>(inner_data_pkgs_, reinitialize_levelset);
	}
	//=================================================================================================//
	void LevelSet::markNearInterface()
	{
		PackageFunctor<void, LevelSetDataPackage> mark_cutcell_by_levelset
			= std::bind(&LevelSet::markNearInterfaceForAPackage, this, _1, _2);
		PackageIterator_parallel<LevelSetDataPackage>(core_data_pkgs_, mark_cutcell_by_levelset);
	}
	//=================================================================================================//
	void LevelSet::markNearInterfaceForAPackage(LevelSetDataPackage* core_data_pkg, Real dt)
	{
		core_data_pkg->markNearInterface();
	}
	//=================================================================================================//
	void LevelSet::redistanceInterface()
	{
		PackageFunctor<void, LevelSetDataPackage> clean_levelset
			= std::bind(&LevelSet::redistanceInterfaceForAPackage, this, _1, _2);
		PackageIterator_parallel<LevelSetDataPackage>(core_data_pkgs_, clean_levelset);
	}
	//=================================================================================================//
	void LevelSet::cleanInterface(bool isSmoothed)
	{
		markNearInterface();
		redistanceInterface();
		reinitializeLevelSet();
		updateNormalDirection();
	}
	//=============================================================================================//
	
}
