/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	base_mesh.h
* @brief 	This is the base classes of mesh, which describe ordered and indexed
*			data sets.  Depending on application, there are different data 
* 			saved on the mesh. The intersection points of mesh lines are called 
*			grid points, the element enclosed by mesh lines (2D) or faces (3D) called 
*			cells. The mesh line or face are also called cell faces. Grid points are
*			also called cell corners.
* @author	Chi ZHang and Xiangyu Hu
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
	class Kernel;

	/** Functor for operation on the mesh data package. */
	template <class ReturnType, class DataPackageType>
	using PackageFunctor = std::function<ReturnType(DataPackageType *, Real)>;
	/** Iterator on a collection of mesh data packages. sequential computing. */
	template <class DataPackageType>
	void PackageIterator(ConcurrentVector<DataPackageType *> &data_pkgs,
						 PackageFunctor<void, DataPackageType> &pkg_functor, Real dt = 0.0)
	{
		for (size_t i = 0; i != data_pkgs.size(); ++i)
			pkg_functor(data_pkgs[i], dt);
	};
	/** Iterator on a collection of mesh data packages. parallel computing. */
	template <class DataPackageType>
	void PackageIterator_parallel(ConcurrentVector<DataPackageType *> &data_pkgs,
								  PackageFunctor<void, DataPackageType> &pkg_functor, Real dt = 0.0)
	{
		parallel_for(
			blocked_range<size_t>(0, data_pkgs.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					pkg_functor(data_pkgs[i], dt);
				}
			},
			ap);
	};
	/** Package iterator for reducing. sequential computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReducePackageIterator(ConcurrentVector<DataPackageType *> &data_pkgs, ReturnType temp,
									 PackageFunctor<ReturnType, DataPackageType> &reduce_pkg_functor,
									 ReduceOperation &reduce_operation, Real dt = 0.0)
	{
		for (size_t i = 0; i < data_pkgs.size(); ++i)
		{
			temp = reduce_operation(temp, reduce_pkg_functor(data_pkgs[i], dt));
		}
		return temp;
	};
	/** Package iterator for reducing. parallel computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReducePackageIterator_parallel(ConcurrentVector<DataPackageType *> &data_pkgs, ReturnType temp,
											  PackageFunctor<ReturnType, DataPackageType> &reduce_pkg_functor,
											  ReduceOperation &reduce_operation, Real dt = 0.0)
	{
		return parallel_reduce(
			blocked_range<size_t>(0, data_pkgs.size()),
			temp, [&](const blocked_range<size_t> &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = reduce_operation(temp, reduce_pkg_functor(data_pkgs[i], dt));
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y) -> ReturnType
			{
				return reduce_operation(x, y);
			});
	};

	/**
	 * @class BaseDataPackage
	 * @brief Abstract base class for a data package
	 * whose data are defined on the grids of a small mesh patch.
	 * note tha ADDRS_SIZE = PKG_SIZE + 2 * pkg_addrs_buffer_;
	 * Also note that, while the mesh lower bound locates the first data address,
	 * the data lower bound locates the first data.
	 */
	template <int PKG_SIZE, int ADDRS_SIZE>
	class BaseDataPackage : public BaseMesh
	{
	public:
		Vecd data_lower_bound_; /**< lower bound coordinate for the data as reference */
		Vecu pkg_index_;		/**< index of the inner packages in the mesh, 0 for far-field packages. */
		bool is_inner_pkg_;		/**< If true, its data saved in memory pool. */
		/** define package data type */
		template <typename DataType>
		using PackageData = PackageDataMatrix<DataType, PKG_SIZE>;
		/** define package data address type */
		template <typename DataType>
		using PackageDataAddress = PackageDataMatrix<DataType *, ADDRS_SIZE>;
		/** define matrix data for temporary usage*/
		template <typename DataType>
		using PackageTemporaryData = PackageDataMatrix<DataType, ADDRS_SIZE>;

		BaseDataPackage() : BaseMesh(Vecu(ADDRS_SIZE)),
							data_lower_bound_(0), pkg_index_(0), is_inner_pkg_(false){};
		virtual ~BaseDataPackage(){};

		constexpr int PackageSize() { return PKG_SIZE; };
		constexpr int AddressSize() { return ADDRS_SIZE; };
		constexpr int AddressBufferWidth() { return (ADDRS_SIZE - PKG_SIZE) / 2; };
		constexpr int OperationUpperBound() { return PKG_SIZE + AddressBufferWidth(); };
		/** initialize package mesh geometric information. */
		void initializePackageGeometry(Vecd &pkg_lower_bound, Real data_spacing)
		{
			mesh_lower_bound_ = pkg_lower_bound - Vecd(data_spacing) * ((Real)AddressBufferWidth() - 0.5);
			grid_spacing_ = data_spacing;
			data_lower_bound_ = pkg_lower_bound + Vecd(data_spacing) * 0.5;
		};
		/** This function probes by applying Bi and tri-linear interpolation within the package. */
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
		/** initialize package data address within a derived class constructor */
		template <typename DataType>
		void initializePackageDataAddress(PackageData<DataType> &pkg_data,
										  PackageDataAddress<DataType> &pkg_data_addrs);
		/** assign address for a package data when the package is an inner one */
		template <typename DataType>
		void assignPackageDataAddress(PackageDataAddress<DataType> &pkg_data_addrs, Vecu &addrs_index,
									  PackageData<DataType> &pkg_data, Vecu &data_index);
		/** obtain averaged value at a corner of a data cell */
		template <typename DataType>
		DataType CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Veci addrs_index, Veci corner_direction);
	};

    /**
     * @class MeshWithDataPackages
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
	 * Note that a data package should be not near the mesh bound, otherwise one will encouter the error "out of range". 
     */ 
	template <class MeshFieldType, class DataPackageType>
	class MeshWithDataPackages : public MeshFieldType, public Mesh
	{
	public:
		MyMemoryPool<DataPackageType> data_pkg_pool_;		  /**< memory pool for all packages in the mesh. */
		MeshDataMatrix<DataPackageType *> data_pkg_addrs_;	  /**< Address of data packages. */
		ConcurrentVector<DataPackageType *> inner_data_pkgs_; /**< Inner data packages which is able to carry out spatial operations. */

		virtual void allocateMeshDataMatrix() override; /**< allocate memories for addresses of data packages. */
		virtual void deleteMeshDataMatrix() override;	/**< delete memories for addresses of data packages. */

		template <typename... Args>
		explicit MeshWithDataPackages(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, Args &&...args)
			: MeshFieldType(std::forward<Args>(args)...),
			  Mesh(tentative_bounds, DataPackageType().PackageSize() * data_spacing, buffer_size),
			  data_spacing_(data_spacing),
			  pkg_size_((int)DataPackageType().PackageSize()),
			  pkg_addrs_buffer_((int)DataPackageType().AddressBufferWidth()),
			  pkg_operations_(pkg_size_ + pkg_addrs_buffer_),
			  pkg_addrs_size_(pkg_size_ + 2 * pkg_addrs_buffer_),
			  global_mesh_(this->mesh_lower_bound_ + Vecd(data_spacing) * 0.5, data_spacing, this->number_of_cells_ * pkg_size_)
		{
			allocateMeshDataMatrix();
		};
		virtual ~MeshWithDataPackages() { deleteMeshDataMatrix(); };

		/** This function probe a mesh value */
		template <class DataType, typename PackageDataAddressType, PackageDataAddressType DataPackageType::*MemPtr>
		DataType probeMesh(const Vecd &position);

	protected:
		Real data_spacing_;									/**< spacing of data in the data packages*/
		int pkg_size_;										/**< the size of the data package matrix*/
		int pkg_addrs_buffer_;								/**< the size of address buffer, a value less than the package size. */
		int pkg_operations_;								/**< the size of operation loops. */
		int pkg_addrs_size_;								/**< the size of address matrix in the data packages. */
		StdVec<DataPackageType *> singular_data_pkgs_addrs; /**< singular data packages. prodvied for far field condition. */
		std::mutex mutex_my_pool;							/**< mutex exclusion for memory pool */
		BaseMesh global_mesh_;								/**< the mesh for the locations of all possible data points. */

		virtual void initializeDataInACell(const Vecu &cell_index, Real dt) = 0;
		virtual void initializeAddressesInACell(const Vecu &cell_index, Real dt) = 0;
		/** This function tag if a data package is inner package. */
		virtual void tagACellIsInnerPackage(const Vecu &cell_index, Real dt) = 0;
		/** This function initialize the data packages with external information */
		virtual void initializeDataPackages() = 0;

		/** This function find the value of data from its index from global mesh. */
		template <typename DataType, typename PackageDataType, PackageDataType DataPackageType::*MemPtr>
		DataType DataValueFromGlobalIndex(Vecu global_grid_index);
		void initializePackageAddressesInACell(Vecu cell_index);
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
#endif //MESH_WITH_DATA_PACKAGES_H