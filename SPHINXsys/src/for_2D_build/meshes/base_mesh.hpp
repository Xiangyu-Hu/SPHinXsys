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
	void BaseDataPackage<PackageDataType>::allocateMeshDataMatrix()
	{
		Allocate2dArray(pkg_data_, Vecu(pkg_size_));
		Allocate2dArray(pkg_data_addrs_, number_of_grid_points_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	void BaseDataPackage<PackageDataType>::deleteMeshDataMatrix()
	{
		Delete2dArray(pkg_data_, Vecu(pkg_size_));
		Delete2dArray(pkg_data_addrs_, number_of_grid_points_);
	}
	//=================================================================================================//
	template<class PackageDataType>
	template<class DataType, DataType PackageDataType:: * MemPtr>
	DataType BaseDataPackage<PackageDataType>
		::probeDataPackage(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Vecd& position)
	{
		Vecu grid_idx = GridIndexFromPosition(position);
		Vecd grid_pos = GridPositionFromIndex(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear
			= pkg_data_addrs[grid_idx[0]][grid_idx[1]]->*MemPtr * beta[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]]->*MemPtr * alpha[0] * beta[1]
			+ pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1]->*MemPtr * beta[0] * alpha[1]
			+ pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1]->*MemPtr * alpha[0] * alpha[1];

		return  bilinear;
	}
	//=================================================================================================//
	template<class PackageDataType>
	template<Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
	void BaseDataPackage<PackageDataType>
		::computeGradient(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Real dt)
	{
		for (size_t i = 1; i != pkg_size_ + 1; ++i)
			for (size_t j = 1; j != pkg_size_ + 1; ++j)
			{
				Real dphidx = (pkg_data_addrs[i + 1][j]->*MemPtrSrc - pkg_data_addrs[i - 1][j]->*MemPtrSrc);
				Real dphidy = (pkg_data_addrs[i][j + 1]->*MemPtrSrc - pkg_data_addrs[i][j - 1]->*MemPtrSrc);
				pkg_data_addrs[i][j]->*MemPtrTrg = Vecd(dphidx, dphidy);
			}
	}
	//=================================================================================================//
	template<class PackageDataType>
	template<Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
	void BaseDataPackage<PackageDataType>
		::computeNormalizedGradient(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Real dt)
	{
		for (size_t i = 1; i != pkg_size_ + 1; ++i)
			for (size_t j = 1; j != pkg_size_ + 1; ++j)
			{
				Real dphidx = (pkg_data_addrs[i + 1][j]->*MemPtrSrc - pkg_data_addrs[i - 1][j]->*MemPtrSrc);
				Real dphidy = (pkg_data_addrs[i][j + 1]->*MemPtrSrc - pkg_data_addrs[i][j - 1]->*MemPtrSrc);
				Real norm = sqrt(dphidx * dphidx + dphidy * dphidy) + TinyReal;
				pkg_data_addrs[i][j]->*MemPtrTrg = Vecd(dphidx, dphidy) / norm;
			}
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
	PackageDataType MeshWithDataPackages<BaseMeshType, DataPackageType>::DataValueFromGlobalIndex(Vecu global_data_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
		for (int n = 0; n != 2; n++)
		{
			size_t cell_index_in_this_direction = global_data_index[n] / pkg_size_;
			pkg_index_[n] = cell_index_in_this_direction;
			local_data_index[n] = global_data_index[n] - cell_index_in_this_direction * pkg_size_;
		}
		return data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]]->pkg_data_[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> template<class PackageDataType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::initializePackageAddressesInACell(Vecu cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->is_inner_pkg_) {
			MeshDataMatrix<PackageDataType*>& pkg_data_addrs = data_pkg->pkg_data_addrs_;
			for (int l = 0; l != pkg_addrs_size_; ++l)
				for (int m = 0; m != pkg_addrs_size_; ++m) {
					pair<int, int>  x_pair = CellShiftAndDataIndex(l);
					pair<int, int>  y_pair = CellShiftAndDataIndex(m);
					pkg_data_addrs[l][m]
						= &(data_pkg_addrs_[i + x_pair.first][j + y_pair.first]->pkg_data_[x_pair.second][y_pair.second]);
				}
		}
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::allocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::deleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, class PackageDataType, DataType PackageDataType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::probeMesh(Vecd& position)
	{
		Vecu grid_index = BaseMeshType::GridIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		return data_pkg->is_inner_pkg_ ?
			data_pkg->DataPackageType::template probeDataPackage<DataType, MemPtr>(data_pkg->pkg_data_addrs_, position)
			: data_pkg->pkg_data_addrs_[0][0]->*MemPtr;
	}
	//=================================================================================================//
}
//=================================================================================================//
