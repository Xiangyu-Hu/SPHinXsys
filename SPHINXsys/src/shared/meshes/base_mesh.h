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
* @version	0.1
* @version  0.2.0
* 			Now narrow bounded levelset mesh is added to replace the whole domain background levelset mesh. 
*/

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "my_memory_pool.h"
#include <fstream>
#include <algorithm>

#include <functional>
using namespace std::placeholders;

using namespace std;

namespace SPH
{
	/** Functor for operation on the mesh. */
	typedef std::function<void(Vecu, Real)> MeshFunctor;
	/** Functor for operation on the mesh data package. */
	template <class ReturnType, class DataPackageType>
	using PackageFunctor = std::function<ReturnType(DataPackageType*, Real)>;
	/** Iterator on the mesh by looping index. sequential computing. */
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);
	/** Iterator on the mesh by looping index. parallel computing. */
	void MeshIterator_parallel(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);
	/** Iterator on a collection of mesh data packages. sequential computing. */
	template <class DataPackageType>
	void PackageIterator(ConcurrentVector<DataPackageType*> inner_data_pkgs,
		PackageFunctor<void, DataPackageType>& pkg_functor, Real dt = 0.0)
	{
		for (size_t i = 0; i != inner_data_pkgs.size(); ++i)
			pkg_functor(inner_data_pkgs[i], dt);

	};
	/** Iterator on a collection of mesh data packages. parallel computing. */
	template <class DataPackageType>
	void PackageIterator_parallel(ConcurrentVector<DataPackageType*> inner_data_pkgs,
		PackageFunctor<void, DataPackageType>& pkg_functor, Real dt = 0.0)
	{
		parallel_for(blocked_range<size_t>(0, inner_data_pkgs.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					pkg_functor(inner_data_pkgs[i], dt);
				}
			}, ap);
	};
	/** Package iterator for reducing. sequential computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReducePackageIterator(ConcurrentVector<DataPackageType*> inner_data_pkgs, ReturnType temp,
		PackageFunctor<ReturnType, DataPackageType>& reduce_pkg_functor, ReduceOperation& reduce_operation, Real dt = 0.0)
	{
		for (size_t i = 0; i < inner_data_pkgs.size(); ++i)
		{
			temp = reduce_operation(temp, reduce_functor(inner_data_pkgs[i], dt));
		}
		return temp;
	};
	/** Package iterator for reducing. parallel computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReducePackageIterator_parallel(ConcurrentVector<DataPackageType*> inner_data_pkgs, ReturnType temp,
		PackageFunctor<ReturnType, DataPackageType>& reduce_pkg_functor, ReduceOperation& reduce_operation, Real dt = 0.0) {
		return parallel_reduce(blocked_range<size_t>(0, inner_data_pkgs.size()),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i) {
					temp0 = reduce_operation(temp, reduce_functor(inner_data_pkgs[i], dt));
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y)->ReturnType {
				return reduce_operation(x, y);
			}
			);
	};

	/**
	 * @class BaseMesh
	 * @brief Abstract base class for all meshes
	 */
	class BaseMesh
	{
	protected:
		Vecd mesh_lower_bound_; 		/**< mesh lower bound as reference coordinate */
		Real grid_spacing_; 			/**< grid_spacing */
		Vecu number_of_grid_points_; 	/**< number of grid points by dimension */

	public:
		/** default constructors */
		BaseMesh() : mesh_lower_bound_(0),
			grid_spacing_(1.0), number_of_grid_points_(1) {};
		/** Constructors */
		BaseMesh(Vecu number_of_grid_points /** Number of grid points in each direction. */
		)
			: mesh_lower_bound_(0), grid_spacing_(1.0),
			number_of_grid_points_(number_of_grid_points) {};
		virtual ~BaseMesh() {};

		/** Return the mesh lower bound. */
		Vecd MeshLowerBound() { return mesh_lower_bound_; };
		/** Return the grid spacing. */
		Real GridSpacing() { return grid_spacing_; };
		/** Return the number of grid points in each direction. */
		Vecu NumberOfGridPoints() { return number_of_grid_points_; };

		/** find grid indexes from point position */
		/**
		 *@brief This function find grid indexes from point position/
		 *@param[in] position(Vecd) inquiry position.
		 *@param[out] (Vecu) grid index.
		 */
		Vecu GridIndexFromPosition(Vecd& position);
		/**
		 *@brief This function find grid position from indexes
		 *@param[in] position(Vecd) inquiry grid index
		 *@param[out] (Vecd) grid position.
		 */
		Vecd GridPositionFromIndex(Vecu grid_index);
		/**
		 *@brief This function convert 1d vector index to mesh index.
		 *@param[in] mesh_size(Vecu) mesh_size in each direction.
		 *@param[in] i(size_t) 1D index
		 *@param[out] (Vecu) 2 or 3 D index
		 */
		Vecu transfer1DtoMeshIndex(Vecu mesh_size, size_t i);
		/** convert mesh index to 1d vector index. */
		/**
		 *@brief This function convert mesh index to 1d vector index.
		 *@param[in] mesh_size(Vecu) Mesh size in each direction.
		 *@param[in] mesh_index Mesh index in each direction.
		 *@param[out] (size_t) 1D index.
		 */
		size_t transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index);

		/** allcate memories for the mesh data matrix*/
		virtual void allocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() = 0;
	};

	/**
	 * @class Mesh
	 * @brief Abstract base class for defining basic mesh properties.
	 * The mesh is proposed for several functions.
	 * First, it is used in cell linked list for neighbor search.
	 * Second, it is used for background maps such as level sets.
	 * This class is the counterpart of the class particles.
	 */
	class Mesh : public BaseMesh
	{
	protected:
		/** buffer size to avoid bound check.*/
		size_t buffer_size_;	/**< buffer size to avoid bound check.*/
		Real cell_spacing_; 	/**< cell_spacing */
		Vecu number_of_cells_; 	/**< number of cells by dimension */
		/**
		 *@brief This function set the mesh lower bound including the buffer region.
		 *@param[in] lower_bound(Vecd) Input mesh lower bound.
		 *@param[in] grid_spacing(Real) Input grid spacing.
		 *@param[in] buffer_size(Real)  Buffersize to extened the mesh from the physical bound of a body.
		 */
		void setMeshLowerBound(Vecd lower_bound, Real grid_spacing, size_t buffer_size);
		/**
		 *@brief This function compute number of total cells
		 *@param[in] lower_bound(Vecd) Input mesh lower bound.
		 *@param[in] grid_spacing(Real) Input grid spacing.
		 *@param[in] buffer_size(Real)  Buffersize to exted the mesh from physical domain of a body or something.
		 */
		Vecu calcNumberOfCells(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size);
		/**
		 *@brief This function compute number of total grid points form total cells
		 *@param[in] number_of_cells(Vecu) Number of cell in each direction.
		 *@param[out] (Vecu)  number of total grid points
		 */
		Vecu NumberOfGridPointsFromNumberOfCells(Vecu number_of_cells) { return number_of_cells + Vecu(1); };
		/**
		 *@brief This function compute number of total cells form total grid points
		 *@param[in] number_of_grid_points(Vecu) Number of grid points in each direction.
		 *@param[in] number of total cells.
		 */
		Vecu NumberOfCellsFromNumberOfGridPoints(Vecu number_of_grid_points) { return number_of_grid_points - Vecu(1); };

		/** copy mesh properties from another mesh. */
		void copyMeshProperties(Mesh* another_mesh);
		/**
		 *@brief This function shift position between cell and grid positions.
		 *@param[in] cell_position(Vecd) Cell position.
		 *@param[out] (Vecd) shifted position.
		 */
		Vecd GridPositionFromCellPosition(Vecd& cell_position)
		{
			return cell_position - Vecd(0.5 * cell_spacing_);
		};
	public:
		/** Constructor using domain information. */
		Mesh(Vecd lower_bound, 		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing,  	/**< Grid spacing. */
			size_t buffer_size = 0 /**< Buffer size. */
		);
		/** Constructor using mesh information directly. */
		Mesh(Vecd mesh_lower_bound, /**< Mesh lower bound. */
			Vecu number_of_cells,  /**< Mesh upper bound. */
			Real cell_spacing 		/**< Mesh cell spacing. */
		);
		virtual ~Mesh() {};

		/** Return the cell spacing. */
		Real CellSpacing() { return cell_spacing_; };
		/** Return the number of cells. */
		Vecu NumberOfCells() { return number_of_cells_; };
		/** Retrun the buffer size. */
		size_t MeshBufferSize() { return buffer_size_; };
		/**
		 * @brief This function check whether a position well within in the mesh bounds
		 * @param[in] position(Vecd) Input position.
		 */
		bool isWithinMeshBound(Vecd position);
		/** find cell indexes from point position */
		Vecu CellIndexesFromPosition(Vecd& position);
		/** Find cell position from indexes.
		  * It is the position shift in the upper-right direction half grid size */
		Vecd CellPositionFromIndexes(Vecu cell_indexes);

		/** output mesh data for Paraview visualization */
		virtual void writeMeshToVtuFile(ofstream& output_file) = 0;
		/** output mesh data for Tecplot visualization */
		virtual void writeMeshToPltFile(ofstream& output_file) = 0;
	};

	/**
	 * @class LevelSetData
	 * @brief Level set is for describing complex geometrics,
	 * It is the distance to the surface of the geometry
	 * and the direction leads to the nearest point on the surface
	 */
	class LevelSetData
	{
	public:
		/** level set is the signed distance to
		  * an interface, here, the surface of a body */
		Real phi_;
		/** level set unit gradient, to approximate interface normal direction */
		Vecd n_;
		/** level set curvature, to approximate interface curvature */
		Real kappa_;
		/*mark the cells cut by near interface level sets, P for positive; N for negative; S for 0 levelset */
		bool gamma_P_, gamma_N_, gamma_S_;
		/** Default constructor */
		LevelSetData() : phi_(0), n_(0), 
			kappa_(0), gamma_P_(false), gamma_N_(false), gamma_S_(false) {};
		virtual ~LevelSetData() {};
	};
	/**
	 * @class BaseDataPackage
	 * @brief Abstract base class for a data package
	 * which is given by a small mesh patch.
	 */
	template<class PackageDataType>
	class BaseDataPackage : public BaseMesh
	{
	public:
		Vecd data_lower_bound_;	/**< lower bound coordinate for the data as reference */
		size_t pkg_size_; 		/**< the size of data package. */
		size_t pkg_addrs_buffer_; 	/**< the size of data package address buffer, should be less than package_size */
		bool is_inner_pkg_; 	/**< If true, its data saved in memory pool. */

		/** A block of data using continuous memory. */
		MeshDataMatrix<PackageDataType> pkg_data_; 			
		/** address of the mesh data for spatial operatings, only valid for inner packages.
		 * The extra layer of data is from the neighbor packages. */
		MeshDataMatrix<PackageDataType*> pkg_data_addrs_; 	

		/** Constructor with package size information.  */
		BaseDataPackage(size_t pkg_size, size_t addrs_buffer)
			: BaseMesh(Vecu(pkg_size + 2 * addrs_buffer)), data_lower_bound_(0),
			pkg_size_(pkg_size), pkg_addrs_buffer_(addrs_buffer), is_inner_pkg_(false) {
			allocateMeshDataMatrix();
		};
		virtual ~BaseDataPackage() {
			deleteMeshDataMatrix();
		};

		/** This function allocate memory for package data */
		virtual void allocateMeshDataMatrix() override;
		/** This function delete memory for package data */
		virtual void deleteMeshDataMatrix() override;
		/** This function return the size of package*/
		size_t PackageSize() { return pkg_size_; };
		/** This function return the size of package fringe*/
		size_t PackageAddressBuffer() { return pkg_addrs_buffer_; };
		/** initialize package mesh information. */
		void initializePackageGeometry(Vecd& pkg_lower_bound, Real data_spacing) {
			mesh_lower_bound_ = pkg_lower_bound - Vecd(data_spacing * 0.5);;
			grid_spacing_ = data_spacing;
			data_lower_bound_ = pkg_lower_bound + Vecd(data_spacing * 0.5);
		};
		/**
		 *@brief This function initialize the package constructed by default.
		 *@param[in] sph_body(SPHBody) A pointer of SPH Body.
		 */
		virtual void initializeDataPackage(SPHBody* sph_body) = 0;
		/**
		 *@brief This function probe Bi and tri-linear interpolation within the package.
		 *@param[in] pkg_data_addrs(MeshDataMatrix<PackageDataType*>) The data matrix.
		 *@param[in] position(Vecd) The inquiry postion.
		 */
		template<class DataType, DataType PackageDataType:: * MemPtr>
		DataType probeDataPackage(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Vecd& position);
		/**
		 *@brief This function compute gradient transform within data package
		 *@param[in] pkg_data_addrs(MeshDataMatrix<PackageDataType*>) The data matrix.
		 *@param[in] dt(Real) Time step (Not used)
		 */
		template<Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
		void computeGradient(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Real dt = 0.0);
		/**
		 *@brief This function compute normalized gradient transform within data package
		 *@param[in] pkg_data_addrs(MeshDataMatrix<PackageDataType*>) The data matrix.
		 *@param[in] dt(Real) Time step (Not used)
		 */
		template<Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
		void computeNormalizedGradient(MeshDataMatrix<PackageDataType*> pkg_data_addrs, Real dt = 0.0);
	};

	/**
	 * @class LevelSetDataPackage
	 * @brief Fixed memory level set data packed in a package.
	 */
	class LevelSetDataPackage : public BaseDataPackage<LevelSetData>
	{
	protected:
		/**
		 *@brief This Heavside function compute volume fraction is a cell.
		 *@param[in] phi(Real) Level set value.
		 *@param[out] Volume of the cell.
		 */
		Real getHeaviside(Real phi);
		/**
		 *@brief This function get aperture at a grid point for computing patch area in a cell.
		 *@param[in] grid_index(Vecu) Index of grid
		 *@param[out] Aperture.
		 */
		Vecd getAperture(Vecu grid_index);
	public:
		bool is_core_pkg_;	/**< If true, the package is near to zero level set. */
		/** default constructor */
		LevelSetDataPackage() : BaseDataPackage<LevelSetData>(4, 1), is_core_pkg_(false) {}
		virtual ~LevelSetDataPackage() {};
		/** This function initialize the level set package constructed by default. */
		virtual void initializeDataPackage(SPHBody* sph_body) override;
		/**
		 *@brief This function initialize with uniform level set field
		 *@param[in] level_set(Real) Level set value
		 *@param[in] normal_direction(Vecd) Normal direction of levelset field.
		 */
		void initializeWithUniformData(Real level_set, Vecd normal_direction);
		/**
		 *@brief This function get curvation at a grid point.
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[out] return the curvature.
		 */
		Real getCurvature(Vecu grid_index);
		/**
		 *@brief This function get normal direction at a grid point.
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[out] return the nurm.
		 */
		Vecd getNormalDirection(Vecu grid_index);
		/**
		 *@brief This function get area of a surface patch in a cell.
		 *@param[in] cell_index(Vecu) Cell index.
		 *@param[out] return the area.
		 */
		Real getSurfaceInCell(Vecu cell_index);
		/**
		 *@brief This function get volume fraction in a cell
		 *@param[in] cell_index(Vecu) Cell index.
		 *@param[out] return  volume fraction
		 */
		Real getVolumeInCell(Vecu cell_index);
		/**
		 *@brief This function indicate if a grid point near the surafce
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[out] ture or false
		 */
		bool isNearSurfaceGrid(Vecu grid_index);
		/**
		 *@brief This function indicate if a grid point within the narrow band.
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[out] ture or false
		 */
		bool isNarrowBandGrid(Vecu grid_index);
		/**
		 *@brief This function check whether apply reinitialization step
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[out] ture or false
		 */
		bool stepReinitialization(Vecu grid_index);
		/**
		 *@brief This function check whether apply reset levelset by computing distance to another levelset
		 *@param[in] grid_index(Vecu) Grid index.
		 *@param[in] phi(Real) Levelset value
		 *@param[out] ture or false
		 */
		bool redistanceToLevelset(Vecu grid_index, Real phi);
	};

	/**
	 * @class BaseMeshWithDataPackages
	 * @brief Abstract class for a mesh with data packages
	 */
	class BaseMeshWithDataPackages : public Mesh
	{
	protected:
		/**
		 *@brief This function initialize a cell with package data.
		 *@param[in] cell_index(Vecu) Cell index
		 *@param[in] dt(Real) Time step size
		 */
		virtual void initializeDataInACell(Vecu cell_index, Real dt) = 0;
		/**
		 *@brief This function initialize the addresses in a data package.
		 *@param[in] cell_index(Vecu) Cell index
		 *@param[in] dt(Real) Time step size
		 */
		virtual void initializeAddressesInACell(Vecu cell_index, Real dt) = 0;
		/**
		 *@brief This function tag if a data package is inner package.
		 *@param[in] cell_index(Vecu) Cell index
		 *@param[in] dt(Real) Time step size
		 */
		virtual void tagACellIsInnerPackage(Vecu cell_index, Real dt) = 0;
	public:
		/** Constructor using domain information. */
		BaseMeshWithDataPackages(Vecd lower_bound, 		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spacing. */
			size_t buffer_size = 0	/**< Buffer size. */
		)
			: Mesh(lower_bound, upper_bound, grid_spacing, buffer_size) {};
		/** Constructor using mesh information directly. */
		BaseMeshWithDataPackages(Vecd mesh_lower_bound, /**< Lower bound. */
			Vecu number_of_cells, 	/**< Upper bound. */
			Real cell_spacing		/**< Cell spacing. */
		)
			: Mesh(mesh_lower_bound, number_of_cells, cell_spacing) {};
		virtual ~BaseMeshWithDataPackages() {};
		/**
		 *@brief This function initialize the data packages with external information
		 */
		virtual void initializeDataPackages() = 0;
	};

	/**
	  * @class BaseLevelSet
	  * @brief A abstract describes a mesh with level set data packages.
	  */
	class BaseLevelSet : public BaseMeshWithDataPackages
	{
	protected:
		SPHBody* sph_body_;
	public:
		/** Constructor using domain information. */
		BaseLevelSet(Vecd lower_bound, 		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spcaing. */
			size_t buffer_size = 0	/**< BUffer size. */
		);
		/** Constructor using mesh information directly. */
		BaseLevelSet(Vecd mesh_lower_bound, /**< Lower bound. */
			Vecu number_of_cells,  /**< Upper bound. */
			Real cell_spacing		/**< Cell spcaing. */
		);
		virtual ~BaseLevelSet() {};
		/**
		 *@brief This function set the SPH body externally
		 *@param[in] sph_body(SPHBody) Pointer to SPH Body
		 */
		void setSPHBody(SPHBody* sph_body) { sph_body_ = sph_body; };
		/**
		 *@brief This function probe the level set at a off-grid position.
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Real probeLevelSet(Vecd position) = 0;
		/**
		 *@brief This function probe the normal direction at a off-grid position
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Vecd probeNormalDirection(Vecd position) = 0;
		/**update the normal direction */
		virtual void updateNormalDirection() = 0;

	};

	/**
	 * @class MeshWithDataPackages
	 * @brief Abstract class fpr mesh with data packages
	 */
	template<class BaseMeshType = BaseLevelSet, class DataPackageType = LevelSetDataPackage>
	class MeshWithDataPackages : public BaseMeshType
	{
	public:
		Mypool<DataPackageType> data_pkg_pool_; 			 /**< memory pool for all packages in the mesh. */
		MeshDataMatrix<DataPackageType*> data_pkg_addrs_; 	 /**< Address of mesh date packages. */
		ConcurrentVector<DataPackageType*> inner_data_pkgs_; /**< Inner data packages which is able to carry out spatial operations. */

		virtual void allocateMeshDataMatrix() override;	/**< allocate memories for addresses of data packages. */
		virtual void deleteMeshDataMatrix() override; 	/**< delete memories for addresses of data packages. */

		/** Constructor using domain information. */
		MeshWithDataPackages(Vecd lower_bound,		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spacing. */
			size_t buffer_size = 0 /**< Buffer size. */
		)
			: BaseMeshType(lower_bound, upper_bound, grid_spacing, buffer_size)
		{
			pkg_size_ = (int)DataPackageType().PackageSize();
			pkg_addrs_buffer_ = (int)DataPackageType().PackageAddressBuffer();
			pkg_addrs_size_ = pkg_size_ + 2 * pkg_addrs_buffer_;
			data_spacing_ = grid_spacing / (Real)pkg_size_;
			total_number_of_data_points_ = this->number_of_cells_ * pkg_size_;
			allocateMeshDataMatrix();
		}
		/** Constructor using mesh information directly. */
		MeshWithDataPackages(Vecd mesh_lower_bound, /**< Lower bound. */
			Vecu number_of_cells, 	/**< Upper bound. */
			Real cell_spacing		/**< Cell spacing. */
		)
			: BaseMeshType(mesh_lower_bound, number_of_cells, cell_spacing)
		{
			pkg_size_ = (int)DataPackageType().PackageSize();
			pkg_addrs_buffer_ = (int)DataPackageType().PackageAddressBuffer();
			pkg_addrs_size_ = pkg_size_ + 2 * pkg_addrs_buffer_;
			data_spacing_ = this->grid_spacing_ / (Real)pkg_size_;
			total_number_of_data_points_ = this->number_of_cells_ * pkg_size_;
			allocateMeshDataMatrix();
		};
		virtual ~MeshWithDataPackages() { deleteMeshDataMatrix(); };
		/**
		 *@brief This function probe a mesh value
		 *@param[in]  position(Vecd) input position.
		 *@param[in]  (DataType) return the probe data.
		 */
		template<class DataType, class PackageDataType, DataType PackageDataType:: * MemPtr>
		DataType probeMesh(Vecd& position);
		/**
		 *@brief This function compute gradient using central difference scheme.
		 *@param[in]  inner_data_pkg(DataPackageType) input data package.
		 *@param[in]  dt(Real) Not used for spatial gradient.
		 */
		template<class PackageDataType, Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
		void computeGradient(DataPackageType* inner_data_pkg, Real dt = 0.0)
		{
			inner_data_pkg->DataPackageType::template computeGradient<MemPtrSrc, MemPtrTrg>(inner_data_pkg->pkg_data_addrs_, dt);
		};
		/**
		 *@brief This function compute normalized gradient using central difference scheme.
		 *@param[in]  inner_data_pkg(DataPackageType) input data package.
		 *@param[in]  dt(Real) Not used for spatial gradient.
		 */
		template<class PackageDataType, Real PackageDataType:: * MemPtrSrc, Vecd PackageDataType:: * MemPtrTrg>
		void computeNormalizedGradient(DataPackageType* inner_data_pkg, Real dt = 0.0)
		{
			inner_data_pkg->DataPackageType::template computeNormalizedGradient<MemPtrSrc, MemPtrTrg>(inner_data_pkg->pkg_data_addrs_, dt);
		};

	protected:
		/** spacing of data in the data packages*/
		Real data_spacing_;
		/** the size of the data package matrix*/
		int pkg_size_;
		/** the size of address buffer, a value less than the package size. */
		int pkg_addrs_buffer_;
		/** the size of address matrix in the data packages. */
		int pkg_addrs_size_;
		/** total numer of data points in all the packages. */
		Vecu total_number_of_data_points_;
		/** singular data packages. prodvied for far field condition. */
		StdVec<DataPackageType*> singular_data_pkgs_addrs;

		/** find the position data from its global index */
		Vecd DataPositionFromGlobalIndex(Vecu global_data_index)
		{
			Vecd data_position;
			for (int n = 0; n < data_position.size(); n++)
			{
				data_position[n] = this->mesh_lower_bound_[n] + 0.5 * data_spacing_
					+ Real(global_data_index[n]) * data_spacing_;
			}
			return data_position;
		};
		/**
		 *@brief This function find the value of data from its global index
		 *@param[in]  global_data_index(Vecu) input global data index.
		 *@param[out] PackageDataType return the corresponding data package.
		 */
		template<class PackageDataType>
		PackageDataType DataValueFromGlobalIndex(Vecu global_data_index);
		/** initialize the addresses in a data package for one variable. */
		template<class PackageDataType>
		void initializePackageAddressesInACell(Vecu cell_index);
		/** find related cell index and data index for a data package address matrix */
		pair<int, int> CellShiftAndDataIndex(size_t data_addrs_index_component)
		{
			pair<int, int> shift_and_index;
			int signed_date_index = (int)data_addrs_index_component - pkg_addrs_buffer_;
			shift_and_index.first = (signed_date_index + pkg_size_) / pkg_size_ - 1;
			shift_and_index.second = signed_date_index - shift_and_index.first * pkg_size_;
			return shift_and_index;
		}
	};

	/**
	 * @class LevelSet
	 * @brief Mesh with level set data as packages.
	 */
	class LevelSet
		: public MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>
	{
	protected:
		/**
		 *@brief This function initialize level set in a cell.
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void initializeDataInACell(Vecu cell_index, Real dt) override;
		/**
		 *@brief This function initialize the addresses in a data package
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void initializeAddressesInACell(Vecu cell_index, Real dt) override;
		/**
		 *@brief This function tag if a data package is inner package.
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void tagACellIsInnerPackage(Vecu cell_index, Real dt) override;
	public:
		/** Core packages which are near to zero level set. */
		ConcurrentVector<LevelSetDataPackage*> core_data_pkgs_;

		/** Constructor using domain information. */
		LevelSet(Vecd mesh_lower_bound, /**< Lower bound. */
			Vecu number_of_cells,  /**< Number of cells. */
			Real cell_spacing);
		/** Constructor using mesh information directly. */
		LevelSet(Vecd lower_bound, 		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spcaing. */
			size_t buffer_size = 0  /**< Buffer size. */
		);
		/** Constructor using domain and sph body information. */
		LevelSet(SPHBody* sph_body,  	/**< Link to SPH body. */
			Vecd lower_bound,      /**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spcaing. */
			size_t buffer_size = 0 /**< Buffer size. */
		);
		virtual ~LevelSet() {};

		/**
		 *@brief This function initialize the Levelset data package.
		 */
		virtual void initializeDataPackages() override;
		/**
		 *@brief This function probe the level set at a off-grid position.
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Real probeLevelSet(Vecd position) override;
		/**
		 *@brief This function probe the normal direction at a off-grid position
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Vecd probeNormalDirection(Vecd position) override;
		/**
		 *@brief This function update the norm of levelset field using central difference scheme.
		 */
		virtual void updateNormalDirection() override;
		/**
		 *@brief This function output mesh data for Paraview
		 *@param[out] output_file(ofstream) output ofstream.
		 */
		virtual void writeMeshToVtuFile(ofstream& output_file) override {};
		/**
		 *@brief This function output mesh data for Tecplot visualization
		 *@param[out] output_file(ofstream) output ofstream.
		 */
		virtual void writeMeshToPltFile(ofstream& output_file) override;
	};
	/**
	 * @class MultilevelMesh
	 * @brief Multi level Meshes with multi resolution mesh data
	 */
	template<class BaseMeshType = BaseLevelSet, class MeshLevelType = LevelSet>
	class MultilevelMesh : public BaseMeshType
	{
	protected:
		/** total levels of the storage pyramid.*/
		size_t total_levels_;
		/**cell spacing for the mesh levels*/
		StdVec<Real> cell_spacing_levels_;
		/** number of cells by dimension */
		StdVec<Vecu> number_of_cells_levels_;

	public:
		/** point to every level. */
		StdVec<MeshLevelType*> mesh_levels_;

		MultilevelMesh(Vecd lower_bound, Vecd upper_bound,
			Real reference_cell_spacing, size_t total_levels, size_t buffer_size = 0)
			: BaseMeshType(lower_bound, upper_bound, reference_cell_spacing, buffer_size),
			total_levels_(total_levels) {

			/** build the zero level mesh first.*/
			int middle_level = ((int)total_levels - 1) / 2;
			Real zero_level_cell_spacing = reference_cell_spacing * powern(2.0, middle_level);
			cell_spacing_levels_.push_back(zero_level_cell_spacing);
			MeshLevelType* zero_level_mesh
				= new MeshLevelType(lower_bound, upper_bound, zero_level_cell_spacing, buffer_size);
			mesh_levels_.push_back(zero_level_mesh);
			Vecu zero_level_number_of_cells = zero_level_mesh->NumberOfCells();
			/** copy zero evel mesh perperties to this. */
			this->copyMeshProperties(zero_level_mesh);

			/** other levels. */
			for (size_t level = 1; level != total_levels; ++level) {
				Real cell_spacing = zero_level_cell_spacing * powern(0.5, (int)level);
				cell_spacing_levels_.push_back(cell_spacing);
				Vecu number_of_cells = zero_level_number_of_cells * powern(2, (int)level);
				MeshLevelType* mesh_level
					= new MeshLevelType(lower_bound, number_of_cells, cell_spacing);
				mesh_levels_.push_back(mesh_level);
			}
		};
		virtual ~MultilevelMesh() {};
	};

	/**
	 * @class MultiresolutionLevelSet
	 * @brief Multi level Meshes for level set data packages
	 */
	class MultiresolutionLevelSet : public MultilevelMesh<BaseLevelSet, LevelSet>
	{
	protected:
		/**the body whose geometry is described by the level set. */
		SPHBody* sph_body_;
	public:
		MultiresolutionLevelSet(SPHBody* sph_body, Vecd lower_bound, Vecd upper_bound,
			Real reference_cell_spacing, size_t total_levels, size_t buffer_size = 0);
		virtual ~MultiresolutionLevelSet() {};
	};

}

