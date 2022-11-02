/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	mesh_with_data_packages.h
 * @brief 	This class is designed to save memory and increase computational efficiency
 *	on mesh. //TODO: the connection between successive meshes in refined mesh should enhanced.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_H
#define MESH_WITH_DATA_PACKAGES_H

#include "base_mesh.h"
#include "my_memory_pool.h"

#include <fstream>
#include <algorithm>
#include <mutex>
#include <functional>
using namespace std::placeholders;

namespace SPH
{
	/** Iterator on a collection of mesh data packages. sequential computing. */
	template <class DataPackageType, typename LocalFunction, typename... Args>
	void package_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
					 const LocalFunction &local_function, Args &&...args)
	{
		for (size_t i = 0; i != data_pkgs.size(); ++i)
			local_function(i);
	};
	/** Iterator on a collection of mesh data packages. parallel computing. */
	template <class DataPackageType, typename LocalFunction, typename... Args>
	void package_parallel_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
							  const LocalFunction &local_function, Args &&...args)
	{
		parallel_for(
			blocked_range<size_t>(0, data_pkgs.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_function(i);
				}
			},
			ap);
	};

	/**
	 * @class BaseDataPackage
	 * @brief Abstract base class for a data package,
	 * by which the data in a derived class can be on- or off-grid.
	 * The data package can be defined in a cell of a background mesh so the pkg_index is
	 * the cell location on the mesh.
	 * TODO: The class will be enriched with general methods for all data packages.
	 */
	class BaseDataPackage
	{
	public:
		Vecu pkg_index_;	/**< index of this data package on the background mesh, Vecu(0) if it is not on the mesh. */
		bool is_inner_pkg_; /**< If true, its data package is on the background mesh. */

		BaseDataPackage() : pkg_index_(0), is_inner_pkg_(false){};
		virtual ~BaseDataPackage(){};
	};

	/**
	 * @class GridDataPackage
	 * @brief Abstract base class for a data package
	 * whose data are defined on the grids of a small mesh patch.
	 * note tha ADDRS_SIZE = PKG_SIZE + 2 * pkg_addrs_buffer_;
	 * Also note that, while the mesh lower bound locates the first data address,
	 * the data lower bound locates the first data.
	 */
	template <int PKG_SIZE, int ADDRS_SIZE>
	class GridDataPackage : public BaseDataPackage, public BaseMesh
	{
	public:
		template <typename DataType>
		using PackageData = PackageDataMatrix<DataType, PKG_SIZE>;
		template <typename DataType>
		using PackageDataAddress = PackageDataMatrix<DataType *, ADDRS_SIZE>;
		GeneralDataAssemble<PackageData> all_pkg_data_;
		GeneralDataAssemble<PackageDataAddress> all_pkg_data_addrs_;

		/** Matrix data for temporary usage. Note that it is array with ADDRS_SIZE.  */
		template <typename DataType>
		using PackageTemporaryData = PackageDataMatrix<DataType, ADDRS_SIZE>;

		GridDataPackage() : BaseDataPackage(), BaseMesh(Vecu(ADDRS_SIZE)){};
		virtual ~GridDataPackage(){};

		constexpr int PackageSize() { return PKG_SIZE; };
		constexpr int AddressSize() { return ADDRS_SIZE; };
		constexpr int AddressBufferWidth() { return (ADDRS_SIZE - PKG_SIZE) / 2; };
		constexpr int OperationUpperBound() { return PKG_SIZE + AddressBufferWidth(); };
		/** lower bound coordinate for the data as reference */
		Vecd DataLowerBound() { return mesh_lower_bound_ + Vecd(grid_spacing_) * (Real)AddressBufferWidth(); };
		/** initialize package mesh geometric information. */
		void initializePackageGeometry(const Vecd &pkg_lower_bound, Real data_spacing)
		{
			mesh_lower_bound_ = pkg_lower_bound - Vecd(data_spacing) * ((Real)AddressBufferWidth() - 0.5);
			grid_spacing_ = data_spacing;
		};
		/** This function probes by applying bi and tri-linear interpolation within the package. */
		template <typename DataType>
		DataType probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position);
		/** This function compute gradient transform within data package */
		template <typename InDataType, typename OutDataType>
		void computeGradient(PackageDataAddress<InDataType> &in_pkg_data_addrs,
							 PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt = 0.0);
		/** This function compute normalized gradient transform within data package  */
		template <typename InDataType, typename OutDataType>
		void computeNormalizedGradient(PackageDataAddress<InDataType> &in_pkg_data_addrs,
									   PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt = 0.0);

	protected:
		/** register a variable defined in a class (can be non-particle class) */
		template <typename DataType>
		void registerPackageData(PackageData<DataType> &pkg_data,
								 PackageDataAddress<DataType> &pkg_data_addrs)
		{
			constexpr int type_index = DataTypeIndex<DataType>::value;
			std::get<type_index>(all_pkg_data_).push_back(&pkg_data);
			std::get<type_index>(all_pkg_data_addrs_).push_back(&pkg_data_addrs);
		};
		/** set the initial package data address within a derived class constructor */
		template <typename DataType>
		struct initializePackageDataAddress
		{
			void operator()(GeneralDataAssemble<PackageData> &all_pkg_data,
							GeneralDataAssemble<PackageDataAddress> &all_pkg_data_addrs);
		};
		DataAssembleOperation<initializePackageDataAddress> initialize_pkg_data_addrs_;
		/** assign address for a package data when the package is an inner one */
		template <typename DataType>
		struct assignPackageDataAddress
		{
			void operator()(GeneralDataAssemble<PackageDataAddress> &all_pkg_data_addrs,
							const Vecu &addrs_index,
							GeneralDataAssemble<PackageData> &all_pkg_data,
							const Vecu &data_index);
		};
		DataAssembleOperation<assignPackageDataAddress> assign_pkg_data_addrs_;

		/** obtain averaged value at a corner of a data cell */
		template <typename DataType>
		DataType CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Veci addrs_index, Veci corner_direction);

	public:
		void initializeSingularDataAddress()
		{
			initialize_pkg_data_addrs_(all_pkg_data_, all_pkg_data_addrs_);
		};

		void assignAllPackageDataAddress(const Vecu &addrs_index, GridDataPackage *src_pkg, const Vecu &data_index)
		{
			assign_pkg_data_addrs_(all_pkg_data_addrs_, addrs_index, src_pkg->all_pkg_data_, data_index);
		};
	};

	/**
	 * @class MeshWithGridDataPackages
	 * @brief Abstract class for mesh with data packages
	 * @details The idea is to save sparse data on a cell-based mesh.
	 * We say sparse data, it means that only in some inner mesh cells there are no trivial data.
	 * A typical example is a level set field which only has meaningful values near the interface,
	 * while the latter is in the inner region of a mesh.
	 * In this class, only some inner mesh cells are filled with data packages.
	 * Each data package is again a mesh, but grid based, where two sets of data are saved on its grid points.
	 * One is the field data of matrices with PKG_SIZE, the other is corresponding address data of matrices with ADDRS_SIZE.
	 * For two neighboring data packages, they share the data in the buffer which is in the overlap region.
	 * The filling of field data is achieved first by the data matrices by the function initializeDataInACell
	 * and then the address matrix by the function initializeAddressesInACell.
	 * All these data packages are indexed by a concurrent vector inner_data_pkgs_.
	 * Note that a data package should be not near the mesh bound, otherwise one will encounter the error "out of range".
	 */
	template <class MeshFieldType, class GridDataPackageType>
	class MeshWithGridDataPackages : public MeshFieldType, public Mesh
	{
	public:
		MyMemoryPool<GridDataPackageType> data_pkg_pool_;	   /**< memory pool for all packages in the mesh. */
		MeshDataMatrix<GridDataPackageType *> data_pkg_addrs_; /**< Address of data packages. */
		ConcurrentVec<GridDataPackageType *> inner_data_pkgs_; /**< Inner data packages which is able to carry out spatial operations. */

		template <typename... Args>
		explicit MeshWithGridDataPackages(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, Args &&...args)
			: MeshFieldType(std::forward<Args>(args)...),
			  Mesh(tentative_bounds, GridDataPackageType().PackageSize() * data_spacing, buffer_size),
			  data_spacing_(data_spacing),
			  global_mesh_(this->mesh_lower_bound_ + Vecd(data_spacing) * 0.5, data_spacing, this->number_of_cells_ * pkg_size_)
		{
			allocateMeshDataMatrix();
		};
		virtual ~MeshWithGridDataPackages() { deleteMeshDataMatrix(); };

		void allocateMeshDataMatrix(); /**< allocate memories for addresses of data packages. */
		void deleteMeshDataMatrix();   /**< delete memories for addresses of data packages. */

		/** This function probe a mesh value */
		template <class DataType, typename PackageDataAddressType, PackageDataAddressType GridDataPackageType::*MemPtr>
		DataType probeMesh(const Vecd &position);
		virtual Real DataSpacing() override { return data_spacing_; };

	protected:
		Real data_spacing_;														  /**< spacing of data in the data packages*/
		const int pkg_size_ = GridDataPackageType().PackageSize();				  /**< the size of the data package matrix*/
		const int pkg_addrs_buffer_ = GridDataPackageType().AddressBufferWidth(); /**< the size of address buffer, a value less than the package size. */
		const int pkg_operations_ = pkg_size_ + pkg_addrs_buffer_;				  /**< the size of operation loops. */
		const int pkg_addrs_size_ = pkg_size_ + 2 * pkg_addrs_buffer_;			  /**< the size of address matrix in the data packages. */
		std::mutex mutex_my_pool;												  /**< mutex exclusion for memory pool */
		BaseMesh global_mesh_;													  /**< the mesh for the locations of all possible data points. */
		/** Singular data packages. provided for far field condition with usually only two values.
		 * For example, when level set is considered. The first value for inner far-field and second for outer far-field */
		StdVec<GridDataPackageType *> singular_data_pkgs_addrs_;

		template <typename... ConstructorArgs>
		void initializeASingularDataPackage(ConstructorArgs &&...args)
		{
			GridDataPackageType *new_data_pkg = data_pkg_pool_.malloc();
			new_data_pkg->registerAllVariables();
			new_data_pkg->initializeSingularData(std::forward<ConstructorArgs>(args)...);
			new_data_pkg->initializeSingularDataAddress();
			singular_data_pkgs_addrs_.push_back(new_data_pkg);
		};

		void assignDataPackageAddress(const Vecu &cell_index, GridDataPackageType *data_pkg);
		GridDataPackageType *DataPackageFromCellIndex(const Vecu &cell_index);
		/** This function find the value of data from its index from global mesh. */
		template <typename DataType, typename PackageDataType, PackageDataType GridDataPackageType::*MemPtr>
		DataType DataValueFromGlobalIndex(const Vecu &global_grid_index);
		void initializePackageAddressesInACell(const Vecu &cell_index);
		/** find related cell index and data index for a data package address matrix */
		std::pair<int, int> CellShiftAndDataIndex(int data_addrs_index_component)
		{
			std::pair<int, int> shift_and_index;
			int signed_date_index = data_addrs_index_component - pkg_addrs_buffer_;
			shift_and_index.first = (signed_date_index + pkg_size_) / pkg_size_ - 1;
			shift_and_index.second = signed_date_index - shift_and_index.first * pkg_size_;
			return shift_and_index;
		}
	};
}
#endif // MESH_WITH_DATA_PACKAGES_H