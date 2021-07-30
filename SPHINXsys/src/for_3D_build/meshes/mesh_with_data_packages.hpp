/**
* @file 	base_mesh.hpp
* @brief 	This is the implementation of the template function and class for base mesh
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef MESH_WITH_DATA_PACKAGES_3D_HPP
#define MESH_WITH_DATA_PACKAGES_3D_HPP

#include "mesh_with_data_packages.h"

//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<class DataType>
	DataType BaseDataPackage<PKG_SIZE, ADDRS_SIZE>
		::probeDataPackage(PackageDataAddress<DataType>& pkg_data_addrs, const Vecd& position)
	{
		Vec3u grid_idx = CellIndexFromPosition(position);
		Vec3d grid_pos = GridPositionFromIndex(grid_idx);
		Vec3d alpha = (position - grid_pos) / grid_spacing_;
		Vec3d beta = Vec3d(1.0) - alpha;

		DataType bilinear_1
			= *pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2]]			* beta[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]]		* alpha[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]]		* beta[0] * alpha[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2]]	* alpha[0] * alpha[1];
		DataType bilinear_2
			= *pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2] + 1]				* beta[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2] + 1]			* alpha[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2] + 1]			* beta[0] * alpha[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2] + 1]		* alpha[0] * alpha[1];
		return  bilinear_1 * beta[2] + bilinear_2 * alpha[2];
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename InDataType, typename OutDataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeGradient(PackageDataAddress<InDataType>& in_pkg_data_addrs,
			PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
				for (int k = 1; k != PKG_SIZE + 1; ++k)
				{
					Real dphidx = (*in_pkg_data_addrs[i + 1][j][k] - *in_pkg_data_addrs[i - 1][j][k]);
					Real dphidy = (*in_pkg_data_addrs[i][j + 1][k] - *in_pkg_data_addrs[i][j - 1][k]);
					Real dphidz = (*in_pkg_data_addrs[i][j][k + 1] - *in_pkg_data_addrs[i][j][k - 1]);
					*out_pkg_data_addrs[i][j][k] = Vecd(dphidx, dphidy, dphidz);
				}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename InDataType, typename OutDataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeNormalizedGradient(PackageDataAddress<InDataType>& in_pkg_data_addrs,
			PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
				for (int k = 1; k != PKG_SIZE + 1; ++k)
				{
					Real dphidx = (*in_pkg_data_addrs[i + 1][j][k] - *in_pkg_data_addrs[i - 1][j][k]);
					Real dphidy = (*in_pkg_data_addrs[i][j + 1][k] - *in_pkg_data_addrs[i][j - 1][k]);
					Real dphidz = (*in_pkg_data_addrs[i][j][k + 1] - *in_pkg_data_addrs[i][j][k - 1]);
					Vecd normal = Vecd(dphidx, dphidy, dphidz);
					*out_pkg_data_addrs[i][j][k] = normal / (normal.norm() + TinyReal);
				}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		initializePackageDataAddress(PackageData<DataType>& pkg_data,
			PackageDataAddress<DataType>& pkg_data_addrs)
	{
		for (int i = 0; i != ADDRS_SIZE; ++i)
			for (int j = 0; j != ADDRS_SIZE; ++j)
				for (int k = 0; k != ADDRS_SIZE; ++k)
			{
				pkg_data_addrs[i][j][k] = &pkg_data[0][0][0];
			}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		assignPackageDataAddress(PackageDataAddress<DataType>& pkg_data_addrs, Vecu& addrs_index,
			PackageData<DataType>& pkg_data, Vecu& data_index)
	{
		pkg_data_addrs[addrs_index[0]][addrs_index[1]][addrs_index[2]] 
			= &pkg_data[data_index[0]][data_index[1]][data_index[2]];
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	DataType  BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		CornerAverage(PackageDataAddress<DataType>& pkg_data_addrs, Veci addrs_index, Veci corner_direction)
	{
		DataType average(0);
		for (int i = 0; i != 2; ++i)
			for (int j = 0; j != 2; ++j)
				for (int k = 0; k != 2; ++k)
				{
					int x_index = addrs_index[0] + i * corner_direction[0];
					int y_index = addrs_index[1] + j * corner_direction[1];
					int z_index = addrs_index[2] + k * corner_direction[2];
					average += *pkg_data_addrs[x_index][y_index][z_index];
				}
		return average * 0.125;
	}
	//=================================================================================================//
	template<class MeshFieldType, class DataPackageType>
	template<typename DataType, typename PackageDataType, PackageDataType DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<MeshFieldType, DataPackageType>::
		DataValueFromGlobalIndex(Vecu global_grid_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
		for (int n = 0; n != 3; n++)
		{
			size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size_;
			pkg_index_[n] = cell_index_in_this_direction;
			local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size_;
		}
		PackageDataType& data = data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]][pkg_index_[2]]->*MemPtr;
		return data[local_data_index[0]][local_data_index[1]][local_data_index[2]];
	}
	//=================================================================================================//
	template<class MeshFieldType, class DataPackageType>
	void MeshWithDataPackages<MeshFieldType, DataPackageType>::initializePackageAddressesInACell(Vecu cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];
		int k = (int)cell_index[2];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j][k];
		if (data_pkg->is_inner_pkg_) {
			for (int l = 0; l != pkg_addrs_size_; ++l)
				for (int m = 0; m != pkg_addrs_size_; ++m)
					for (int n = 0; n != pkg_addrs_size_; ++n) {
						std::pair<int, int>  x_pair = CellShiftAndDataIndex(l);
						std::pair<int, int>  y_pair = CellShiftAndDataIndex(m);
						std::pair<int, int>  z_pair = CellShiftAndDataIndex(n);

						data_pkg->assignAllPackageDataAddress(Vecu(l, m, n),
							data_pkg_addrs_[i + x_pair.first][j + y_pair.first][k + z_pair.first],
							Vecu(x_pair.second, y_pair.second, z_pair.second));
					}
		}
	}
	//=================================================================================================//
	template<class MeshFieldType, class DataPackageType>
	void MeshWithDataPackages<MeshFieldType, DataPackageType>::allocateMeshDataMatrix()
	{
		Allocate3dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template<class MeshFieldType, class DataPackageType>
	void MeshWithDataPackages<MeshFieldType, DataPackageType>::deleteMeshDataMatrix()
	{
        Delete3dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template<class MeshFieldType, class DataPackageType>
	template<class DataType, typename PackageDataAddressType, PackageDataAddressType DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<MeshFieldType, DataPackageType>::probeMesh(const Vecd& position)
	{
        Vecu grid_index = CellIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];
		size_t k = grid_index[2];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j][k];
		PackageDataAddressType& pkg_data_addrs = data_pkg->*MemPtr;
		return data_pkg->is_inner_pkg_ ?
			data_pkg->DataPackageType::template probeDataPackage<DataType>(pkg_data_addrs, position)
			: *pkg_data_addrs[0][0][0];
	}
	//=================================================================================================//
}
//=================================================================================================//
#endif //MESH_WITH_DATA_PACKAGES_3D_HPP