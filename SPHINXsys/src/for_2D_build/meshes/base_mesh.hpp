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
	template<class DataType>
	DataType BaseDataPackage
		::ProbeDataPackage(MeshDataMatrix<DataType*> data_addrs, Vecd& position)
	{
		Vecu grid_idx = GridIndexesFromPosition(position);
		Vecd grid_pos = GridPositionFromIndexes(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear
			= *(data_addrs[grid_idx[0]][grid_idx[1]]) * beta[0] * beta[1]
			+ *(data_addrs[grid_idx[0] + 1][grid_idx[1]]) * alpha[0] * beta[1]
			+ *(data_addrs[grid_idx[0]][grid_idx[1] + 1]) * beta[0] * alpha[1]
			+ *(data_addrs[grid_idx[0] + 1][grid_idx[1] + 1]) * alpha[0] * alpha[1];

		return  bilinear;
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::getValueFromGlobalDataIndex(Vecu global_data_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
			for (int n = 0; n != 2; n++)
			{
				int cell_index_in_this_direction = (int)global_data_index[n] / (int)pkg_size_;
				pkg_index_[n] = cell_index_in_this_direction;
				local_data_index[n]	= (int)global_data_index[n] - cell_index_in_this_direction * (int)pkg_size_;
			}
		return (data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]]->*MemPtr)[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, MeshDataMatrix<DataType*> DataPackageType:: * MemPtrAddrss,
		MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
		void MeshWithDataPackages<BaseMeshType, DataPackageType>::initializeOneVariableAdressesInACell(Vecu cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->is_inner_pkg_) {
			size_t bf_sz = data_pkg->getBufferSize();
			/** inner data*/
			for (size_t l = 0; l != pkg_size_; ++l)
				for (size_t m = 0; m != pkg_size_; ++m) {
					(data_pkg->*MemPtrAddrss)[l + bf_sz][m + bf_sz] = &(data_pkg->*MemPtr)[l][m];
				}
			size_t bd_index = pkg_size_ + bf_sz;
			size_t bd_data_index = pkg_size_ - 1;
			/** left and right sides*/
			for (size_t bf = 0; bf != bf_sz; ++bf)
				for (size_t m = 0; m != pkg_size_; ++m) {
					//left
					(data_pkg->*MemPtrAddrss)[bf][m + bf_sz] = &(data_pkg_addrs_[i - 1][j]->*MemPtr)[bd_data_index - bf][m];
					//right
					(data_pkg->*MemPtrAddrss)[bd_index + bf][m + bf_sz] = &(data_pkg_addrs_[i + 1][j]->*MemPtr)[bf][m];
				}
			/** bottom and upper sides*/
			for (size_t bf = 0; bf != bf_sz; ++bf)
				for (size_t l = 0; l != pkg_size_; ++l) {
					//bottom
					(data_pkg->*MemPtrAddrss)[l + bf_sz][bf] = &(data_pkg_addrs_[i][j - 1]->*MemPtr)[l][bd_data_index - bf];
					//upper
					(data_pkg->*MemPtrAddrss)[l + bf_sz][bd_index + bf] = &(data_pkg_addrs_[i][j + 1]->*MemPtr)[l][bf];
				}
			/** corners*/
			for (size_t bf_lf = 0; bf_lf != bf_sz; ++bf_lf)
				for (size_t bf_bu = 0; bf_bu != bf_sz; ++bf_bu) {
					//left-lower
					(data_pkg->*MemPtrAddrss)[bf_lf][bf_bu] = &(data_pkg_addrs_[i - 1][j - 1]->*MemPtr)[bd_data_index - bf_lf][bd_data_index - bf_bu];
					//right-lower
					(data_pkg->*MemPtrAddrss)[bd_index + bf_lf][bf_bu] = &(data_pkg_addrs_[i + 1][j - 1]->*MemPtr)[bf_lf][bd_data_index - bf_bu];
					//left-upper
					(data_pkg->*MemPtrAddrss)[bf_lf][bd_index + bf_bu] = &(data_pkg_addrs_[i - 1][j + 1]->*MemPtr)[bd_data_index - bf_lf][bf_bu];
					//right-upper
					(data_pkg->*MemPtrAddrss)[bd_index + bf_lf][bd_index + bf_bu] = &(data_pkg_addrs_[i + 1][j + 1]->*MemPtr)[bf_lf][bf_bu];
				}
		}
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::AllocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::DeleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, MeshDataMatrix<DataType*> DataPackageType:: * MemPtrAddrss,
		MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::probeMesh(Vecd& position)
	{
		Vecu grid_index = GridIndexesFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];
	
		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		return data_pkg->is_inner_pkg_ ?
			 data_pkg->ProbeDataPackage(data_pkg->*MemPtrAddrss, position)
			: (data_pkg->*MemPtr)[0][0];
	}
	//=================================================================================================//
}
//=================================================================================================//
