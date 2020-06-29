/**
* @file 	base_mesh.hpp
* @brief 	This is the implementation of the template function and class for base mesh
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_mesh.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template<class PackageDataType>
	void BaseDataPackage<PackageDataType>::AllocateMeshDataMatrix()
	{
		Allocate3dArray(pkg_data_, number_of_grid_points_);
		Allocate3dArray(pkg_data_addrs_, number_of_addrs_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	void BaseDataPackage<PackageDataType>::DeleteMeshDataMatrix()
	{
		Delete3dArray(pkg_data_, number_of_grid_points_);
		Delete3dArray(pkg_data_addrs_, number_of_grid_points_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	template<class DataType, DataType PackageDataType:: * MemPtr>
	DataType BaseDataPackage<PackageDataType>
		::ProbeDataPackage(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Vecd& position)
	{
		Vec3u grid_idx = GridIndexesFromPosition(position);
		Vec3d grid_pos = GridPositionFromIndexes(grid_idx);
		Vec3d alpha = (position - grid_pos) / grid_spacing_;
		Vec3d beta = Vec3d(1.0) - alpha;

		DataType bilinear_1
			= pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2]]->*MemPtr			* beta[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]]->*MemPtr		* alpha[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]]->*MemPtr		* beta[0] * alpha[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2]]->*MemPtr	* alpha[0] * alpha[1];
		DataType bilinear_2
			= pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2] + 1]->*MemPtr * beta[0]		* beta[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2] + 1]->*MemPtr			* alpha[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2] + 1]->*MemPtr			* beta[0] * alpha[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2] + 1]->*MemPtr		* alpha[0] * alpha[1];
		return  bilinear_1 * beta[2] + bilinear_2 * alpha[2];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
	PackageDataType MeshWithDataPackages<BaseMeshType, DataPackageType>::getValueFromGlobalDataIndex(Vecu global_data_index)
	{
		cout << "\n This function MeshWithDataPackages<DataPackageType>::getValueFromGlobalDataIndex is not done in 3D. Exit the program! \n";
		exit(0);
		return PackageDataType();
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
		void MeshWithDataPackages<BaseMeshType, DataPackageType>::initializePackageAdressesInACell(Vecu cell_index)
	{
		cout << "\n This function MeshWithDataPackages<DataPackageType>::initializeOneVariableAdressesInACell is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::AllocateMeshDataMatrix()
	{
		Allocate3dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::DeleteMeshDataMatrix()
	{
        Delete3dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, class PackageDataType, DataType PackageDataType::* MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::probeMesh(Vecd& position)
	{
        Vecu grid_index = BaseMeshType::GridIndexesFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];
		size_t k = grid_index[2];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j][k];
		if (data_pkg->is_inner_pkg_) {
			return data_pkg->DataPackageType::template ProbeDataPackage<DataType, MemPtr>(data_pkg->pkg_data_addrs_, position);
		}
		else return data_pkg->pkg_data_addrs_[0][0][0]->*MemPtr;
	}
	//=================================================================================================//
}
//=================================================================================================//
