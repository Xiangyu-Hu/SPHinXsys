/**
 * @file 	base_mesh.hpp
 * @brief 	This is the implementation of the template function and class for base mesh
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_2D_HPP
#define MESH_WITH_DATA_PACKAGES_2D_HPP

#include "mesh_with_data_packages.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename FunctionOnData>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		for_each_data(const FunctionOnData &function)
	{
		for (int i = 0; i != pkg_size; ++i)
			for (int j = 0; j != pkg_size; ++j)
			{
				function(i, j);
			}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename FunctionOnAddress>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		for_each_addrs(const FunctionOnAddress &function)
	{
		for (int i = pkg_addrs_buffer; i != pkg_ops_end; ++i)
			for (int j = pkg_addrs_buffer; j != pkg_ops_end; ++j)
			{
				function(i, j);
			}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <class DataType>
	DataType GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position)
	{
		Vecu grid_idx = CellIndexFromPosition(position);
		Vecd grid_pos = GridPositionFromIndex(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear = *pkg_data_addrs[grid_idx[0]][grid_idx[1]] * beta[0] * beta[1] +
							*pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]] * alpha[0] * beta[1] +
							*pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1] * beta[0] * alpha[1] +
							*pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1] * alpha[0] * alpha[1];

		return bilinear;
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename InDataType, typename OutDataType>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		computeGradient(const DiscreteVariable<InDataType> &in_variable,
						const DiscreteVariable<OutDataType> &out_variable)
	{
		auto &in_variable_addrs = getPackageDataAddress(in_variable);
		auto &out_variable_addrs = getPackageDataAddress(out_variable);

		for_each_addrs(
			[&](int i, int j)
			{
				Real dphidx = (*in_variable_addrs[i + 1][j] - *in_variable_addrs[i - 1][j]);
				Real dphidy = (*in_variable_addrs[i][j + 1] - *in_variable_addrs[i][j - 1]);
				*out_variable_addrs[i][j] = 0.5 * Vecd(dphidx, dphidy) / grid_spacing_;
			});
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename DataType, typename FunctionByPosition>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		assignByPosition(const DiscreteVariable<DataType> &discrete_variable,
						 const FunctionByPosition &function_by_position)
	{
		auto &pkg_data = getPackageData(discrete_variable);
		for (int i = 0; i != pkg_size; ++i)
			for (int j = 0; j != pkg_size; ++j)
			{
				Vec2d position = DataPosition(Vec2d(i, j));
				pkg_data[i][j] = function_by_position(position);
			}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename DataType>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::initializeAllPackageDataAddress<DataType>::
	operator()(DataContainerAssemble<PackageData> &all_pkg_data,
			   DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs)
	{
		constexpr int type_index = DataTypeIndex<DataType>::value;
		for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
		{
			PackageData<DataType> &pkg_data = std::get<type_index>(all_pkg_data)[l];
			PackageDataAddress<DataType> &pkg_data_addrs = std::get<type_index>(all_pkg_data_addrs)[l];
			for (int i = 0; i != pkg_addrs_size; ++i)
				for (int j = 0; j != pkg_addrs_size; ++j)
				{
					pkg_data_addrs[i][j] = &pkg_data[0][0];
				}
		}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename DataType>
	DataType GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
		CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Veci addrs_index, Veci corner_direction)
	{
		DataType average(0);
		for (int i = 0; i != 2; ++i)
			for (int j = 0; j != 2; ++j)
			{
				int x_index = addrs_index[0] + i * corner_direction[0];
				int y_index = addrs_index[1] + j * corner_direction[1];
				average += *pkg_data_addrs[x_index][y_index];
			}
		return average * 0.25;
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_BUFFER>
	template <typename DataType>
	void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::AssignAllPackageDataAddress<DataType>::
	operator()(DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs,
			   const Vecu &addrs_index,
			   DataContainerAssemble<PackageData> &all_pkg_data,
			   const Vecu &data_index)
	{
		constexpr int type_index = DataTypeIndex<DataType>::value;
		for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
		{
			PackageData<DataType> &pkg_data = std::get<type_index>(all_pkg_data)[l];
			PackageDataAddress<DataType> &pkg_data_addrs = std::get<type_index>(all_pkg_data_addrs)[l];
			pkg_data_addrs[addrs_index[0]][addrs_index[1]] = &pkg_data[data_index[0]][data_index[1]];
		}
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	template <typename DataType>
	DataType MeshWithGridDataPackages<GridDataPackageType>::
		DataValueFromGlobalIndex(const DiscreteVariable<DataType> &discrete_variable,
								 const Vecu &global_grid_index)
	{
		Vecu cell_index_on_mesh_(0);
		Vecu local_data_index(0);
		for (int n = 0; n != 2; n++)
		{
			size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
			cell_index_on_mesh_[n] = cell_index_in_this_direction;
			local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
		}
		auto &data = data_pkg_addrs_[cell_index_on_mesh_[0]][cell_index_on_mesh_[1]]->getPackageData(discrete_variable);
		return data[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	void MeshWithGridDataPackages<GridDataPackageType>::
		initializePackageAddressesInACell(const Vecu &cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		GridDataPackageType *data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->isInnerPackage())
		{
			for (int l = 0; l != pkg_addrs_size; ++l)
				for (int m = 0; m != pkg_addrs_size; ++m)
				{
					std::pair<int, int> x_pair = CellShiftAndDataIndex(l);
					std::pair<int, int> y_pair = CellShiftAndDataIndex(m);
					data_pkg->assignAllPackageDataAddress(
						Vecu(l, m),
						data_pkg_addrs_[i + x_pair.first][j + y_pair.first],
						Vecu(x_pair.second, y_pair.second));
				}
		}
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	void MeshWithGridDataPackages<GridDataPackageType>::allocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	void MeshWithGridDataPackages<GridDataPackageType>::deleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	void MeshWithGridDataPackages<GridDataPackageType>::
		assignDataPackageAddress(const Vecu &cell_index, GridDataPackageType *data_pkg)
	{
		data_pkg_addrs_[cell_index[0]][cell_index[1]] = data_pkg;
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	GridDataPackageType *MeshWithGridDataPackages<GridDataPackageType>::
		DataPackageFromCellIndex(const Vecu &cell_index)
	{
		return data_pkg_addrs_[cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	template <class DataType, typename PackageDataAddressType, PackageDataAddressType GridDataPackageType::*MemPtr>
	DataType MeshWithGridDataPackages<GridDataPackageType>::probeMesh(const Vecd &position)
	{
		Vecu grid_index = CellIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];

		GridDataPackageType *data_pkg = data_pkg_addrs_[i][j];
		PackageDataAddressType &pkg_data_addrs = data_pkg->*MemPtr;
		return data_pkg->isInnerPackage() ? data_pkg->GridDataPackageType::template probeDataPackage<DataType>(pkg_data_addrs, position)
										  : *pkg_data_addrs[0][0];
	}
	//=================================================================================================//
	template <class GridDataPackageType>
	template <class DataType>
	DataType MeshWithGridDataPackages<GridDataPackageType>::
		probeMesh(const DiscreteVariable<DataType> &discrete_variable, const Vecd &position)
	{
		Vecu grid_index = CellIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];

		GridDataPackageType *data_pkg = data_pkg_addrs_[i][j];
		auto &pkg_data_addrs = data_pkg->getPackageDataAddress(discrete_variable);
		return data_pkg->isInnerPackage() ? data_pkg->GridDataPackageType::template probeDataPackage<DataType>(pkg_data_addrs, position)
										  : *pkg_data_addrs[0][0];
	}
	//=================================================================================================//
}
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_2D_HPP