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
		Allocate2dArray(pkg_data_, number_of_grid_points_);
		Allocate2dArray(pkg_data_addrs_, number_of_addrs_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	void BaseDataPackage<PackageDataType>::DeleteMeshDataMatrix()
	{
		Delete2dArray(pkg_data_, number_of_grid_points_);
		Delete2dArray(pkg_data_addrs_, number_of_grid_points_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	template<class DataType, DataType PackageDataType:: * MemPtr>
	DataType BaseDataPackage<PackageDataType>
		::ProbeDataPackage(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Vecd& position)
	{
		Vecu grid_idx = GridIndexesFromPosition(position);
		Vecd grid_pos = GridPositionFromIndexes(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear
			= pkg_data_addrs[grid_idx[0]][grid_idx[1]]->*MemPtr				* beta[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]]->*MemPtr			* alpha[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1]->*MemPtr			* beta[0] * alpha[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1]->*MemPtr		* alpha[0] * alpha[1];

		return  bilinear;
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
	PackageDataType MeshWithDataPackages<BaseMeshType, DataPackageType>::getValueFromGlobalDataIndex(Vecu global_data_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
			for (int n = 0; n != 2; n++)
			{
				int cell_index_in_this_direction = (int)global_data_index[n] / (int)pkg_size_;
				pkg_index_[n] = cell_index_in_this_direction;
				local_data_index[n]	= (int)global_data_index[n] - cell_index_in_this_direction * (int)pkg_size_;
			}
		return data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]]->pkg_data_[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::initializePackageAdressesInACell(Vecu cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->is_inner_pkg_) {
			size_t bf_sz = data_pkg->getBufferSize();
			MeshDataMatrix<PackageDataType*>& pkg_data_addrs = data_pkg->pkg_data_addrs_;
			/** inner data*/
			MeshDataMatrix<PackageDataType>& pkg_data = data_pkg->pkg_data_;
			for (size_t l = 0; l != pkg_size_; ++l)
				for (size_t m = 0; m != pkg_size_; ++m) {
					pkg_data_addrs[l + bf_sz][m + bf_sz] = &(pkg_data[l][m]);
				}
			size_t bd_index = pkg_size_ + bf_sz;
			size_t bd_data_index = pkg_size_ - 1;
			/** left and right sides*/
			MeshDataMatrix<PackageDataType>& pkg_data_l = data_pkg_addrs_[i - 1][j]->pkg_data_;
			MeshDataMatrix<PackageDataType>& pkg_data_r = data_pkg_addrs_[i + 1][j]->pkg_data_;
			for (size_t bf = 0; bf != bf_sz; ++bf)
				for (size_t m = 0; m != pkg_size_; ++m) {
					//left
					pkg_data_addrs[bf][m + bf_sz] = &(pkg_data_l[bd_data_index - bf][m]);
					//right
					pkg_data_addrs[bd_index + bf][m + bf_sz] = &(pkg_data_r[bf][m]);
				}
			/** bottom and upper sides*/
			MeshDataMatrix<PackageDataType>& pkg_data_b = data_pkg_addrs_[i][j - 1]->pkg_data_;
			MeshDataMatrix<PackageDataType>& pkg_data_u = data_pkg_addrs_[i][j + 1]->pkg_data_;
			for (size_t bf = 0; bf != bf_sz; ++bf)
				for (size_t l = 0; l != pkg_size_; ++l) {
					//bottom
					pkg_data_addrs[l + bf_sz][bf] = &(pkg_data_b[l][bd_data_index - bf]);
					//upper
					pkg_data_addrs[l + bf_sz][bd_index + bf] = &(pkg_data_u[l][bf]);
				}
			/** corners*/
			MeshDataMatrix<PackageDataType>& pkg_data_lb = data_pkg_addrs_[i - 1][j - 1]->pkg_data_;
			MeshDataMatrix<PackageDataType>& pkg_data_rb = data_pkg_addrs_[i + 1][j - 1]->pkg_data_;
			MeshDataMatrix<PackageDataType>& pkg_data_lu = data_pkg_addrs_[i - 1][j + 1]->pkg_data_;
			MeshDataMatrix<PackageDataType>& pkg_data_ru = data_pkg_addrs_[i + 1][j + 1]->pkg_data_;
			for (size_t bf_lf = 0; bf_lf != bf_sz; ++bf_lf)
				for (size_t bf_bu = 0; bf_bu != bf_sz; ++bf_bu) {
					//left-lower
					pkg_data_addrs[bf_lf][bf_bu] = &(pkg_data_lb[bd_data_index - bf_lf][bd_data_index - bf_bu]);
					//right-lower
					pkg_data_addrs[bd_index + bf_lf][bf_bu] = &(pkg_data_rb[bf_lf][bd_data_index - bf_bu]);
					//left-upper
					pkg_data_addrs[bf_lf][bd_index + bf_bu] = &(pkg_data_lu[bd_data_index - bf_lf][bf_bu]);
					//right-upper
					pkg_data_addrs[bd_index + bf_lf][bd_index + bf_bu] = &(pkg_data_ru[bf_lf][bf_bu]);
				}
		}
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::AllocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::DeleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, class PackageDataType, DataType PackageDataType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::probeMesh(Vecd& position)
	{
		Vecu grid_index = BaseMeshType::GridIndexesFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];
	
		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		return data_pkg->is_inner_pkg_ ?
			 data_pkg->DataPackageType::template ProbeDataPackage<DataType, MemPtr>(data_pkg->pkg_data_addrs_, position)
			: data_pkg->pkg_data_addrs_[0][0]->*MemPtr;
	}
	//=================================================================================================//
}
//=================================================================================================//
