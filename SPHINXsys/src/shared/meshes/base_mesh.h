/**
* @file 	base_mesh.h
* @brief 	This is the base classes of mesh, which describe ordered and indexed
*			data sets.  Depending on application, there are different data 
* 			saved on the mesh. The intersection points of mesh lines are called 
*			grid points, the element enclosed by mesh lines (2D) or faces (3D) called 
*			cells. The mesh line or face are also called cell faces. Grid points are
*			also called cell corners.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
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
	/** Functor for operation oon the mesh. */
	typedef std::function<void(Vecu, Real)> MeshFunctor;
	/** Functor for operation on the mesh data package. */
	template <class ReturnType>
	using PacakgeFunctor = std::function<ReturnType(Real)>;
	/** Iterator on the mesh by looping index. sequential computing. */
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);
	/** Iterator on the mesh by looping index. parallel computing. */
	void MeshIterator_parallel(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);
	/** Iterator on a collection of mesh data packages. sequential computing. */
	template <class DataPackageType>
	void PackageIterator(ConcurrentVector<DataPackageType*> inner_data_pkgs,
		PacakgeFunctor<void>& pkg_functor, Real dt = 0.0) 
	{
		for (size_t i = 0; i != inner_data_pkgs.size(); ++i)
			inner_data_pkgs[i]->pkg_functor(dt);

	};
	/** Iterator on a collection of mesh data packages. parallel computing. */
	template <class DataPackageType>
	void PackageIterator_parallel(ConcurrentVector<DataPackageType*> inner_data_pkgs,
		PacakgeFunctor<void>& pkg_functor, Real dt = 0.0) 
	{
		parallel_for(blocked_range<size_t>(0, inner_data_pkgs.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					inner_data_pkgs[i]->pkg_functor(dt);
				}
			}, ap);
	};
	/** Package iterator for reducing. sequential computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReducePacakageIterator(ConcurrentVector<DataPackageType*> inner_data_pkgs, ReturnType temp,
		PacakgeFunctor<ReturnType>& reduce_pkg_functor, ReduceOperation& ruduce_operation, Real dt = 0.0) 
	{
		for (size_t i = 0; i < inner_data_pkgs.size(); ++i)
		{
			temp = reduce_operation(temp, inner_data_pkgs[i]->reduce_functor(dt));
		}
		return temp;
	};
	/** Package iterator for reducing. parallel computing. */
	template <class ReturnType, typename ReduceOperation, class DataPackageType>
	ReturnType ReduceMeshIterator_parallel(ConcurrentVector<DataPackageType*> inner_data_pkgs, ReturnType temp,
		PacakgeFunctor<ReturnType>& reduce_fupkg_functor, ReduceOperation& ruduce_operation, Real dt = 0.0) {
		return parallel_reduce(blocked_range<size_t>(0, inner_data_pkgs.size()),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType 
			{
				for (size_t i = r.begin(); i != r.end(); ++i) {
					temp0 = reduce_operation(temp0, inner_data_pkgs[i]->reduce_functor(dt));
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
		/** mesh lower bound as reference coordiante */
		Vecd mesh_lower_bound_;
		/** buffer size to avoid bound check.*/
		size_t buffer_size_;
		/** grid_spacing */
		Real grid_spacing_;
		/** number of grid points by dimension */
		Vecu number_of_grid_points_;

	public:
		/** default constructors */
		BaseMesh() : mesh_lower_bound_(0), 
			buffer_size_(0), grid_spacing_(1.0), number_of_grid_points_(1) {};
		/** Constructors */
		BaseMesh(Vecu number_of_grid_points, size_t buffer_size = 0)
			: mesh_lower_bound_(0), buffer_size_(buffer_size), grid_spacing_(1.0),
			number_of_grid_points_(number_of_grid_points) {};
		virtual ~BaseMesh() {};

		/** accesss protected variables */
		Vecd getMeshLowerBound() { return mesh_lower_bound_; };
		Real getGridSpacing() { return grid_spacing_; };
		Vecu getNumberOfGridPoints() { return number_of_grid_points_; };
		size_t getBufferSize() { return buffer_size_; };

		/** find grid indexes from point poistion */
		Vecu GridIndexesFromPosition(Vecd& position);
		/** find grid poistion from indexes */
		Vecd GridPositionFromIndexes(Vecu grid_indexes);
		/** convert 1d vector index to mesh index. */
		Vecu transfer1DtoMeshIndex(Vecu mesh_size, size_t i);
		/** convert mesh index to 1d vector index. */
		size_t transferMeshIndexTo1D(Vecu mesh_size, Vecu mesh_index);

		/** allcate memories for the mesh data matrix*/
		virtual void AllocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() = 0;
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
		/** cell_spacing */
		Real cell_spacing_;
		/** number of cells by dimension */
		Vecu number_of_cells_;

		/** set the mesh lower bound including the buffer region. */
		void setMeshLowerBound(Vecd lower_bound, Real grid_spacing, size_t buffer_size);
		/** computing number of total lattices */
		Vecu calcNumberOfCells(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size);
		/** computing number of total lattices */
		Vecu getNumberOfGridPoints(Vecu number_of_cells) { return number_of_cells + Vecu(1); };
		Vecu getNumberOfCells(Vecu number_of_grid_points) { return number_of_grid_points - Vecu(1); };

		/** copy mesh properties from another mesh. */
		void copyMeshProperties(Mesh* another_mesh);
		/** shift position between cell and grid positions. */
		Vecd getGridPositionFromCellPosition(Vecd& cell_position) {
			return cell_position - Vecd(0.5 * cell_spacing_); 
		}
	public:
		/** Constructor using domain information. */
		Mesh(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0);
		/** Constructor using mesh information directly. */
		Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		virtual ~Mesh() {};

		/** accesss protected variables */
		Real getCellSpacing() { return cell_spacing_; };
		Vecu getNumberOfCells() { return number_of_cells_; };

		/** check whether a position well within in the mesh bounds */
		bool checkMeshBound(Vecd position);
		/** find cell indexes from point poistion */
		Vecu CellIndexesFromPosition(Vecd& position);
		/** Find cell position from indexes.
		  * It is the position shift in the upper-right direction half grid size */
		Vecd CellPositionFromIndexes(Vecu cell_indexes);

		/** output mesh data for Paraview visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) = 0;
		/** output mesh data for Tecplot visuallization */
		virtual void WriteMeshToPltFile(ofstream &output_file) = 0;
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
		/** displacment to the surface */
		Vecd n_;
		/** curvature */
		Real kappa_;

		/** Default constructor */
		LevelSetData() : phi_(0), n_(0), kappa_(1.0) {};
		/** Constructor for one data point. */
		LevelSetData(Real level_set, Vecd normal_direction);
		virtual ~LevelSetData() {};
	};

	/**
	 * @class MeshBackground
	 * @brief Background mesh gives a level set field.
	 */
	class MeshBackground : public Mesh
	{
		/** base mesh data */
		MeshDataMatrix<LevelSetData> mesh_background_data_;
	public:
		MeshBackground(Vecd lower_bound, Vecd upper_bound, 
			Real grid_spacing, size_t buffer_size = 0);
		virtual ~MeshBackground() { DeleteMeshDataMatrix(); };

		/** allocate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() override;

		/** initialize level set and displacement to surface
		  * for body region geometry */
		void InitializeLevelSetData(SPHBody &body);
		void ComputeCurvatureFromLevelSet(SPHBody &body);
		/** probe the mesh data */
		Vecd ProbeNormalDirection(Vecd Point);
		Real ProbeLevelSet(Vecd Point);
		Real ProbeCurvature(Vecd Point);

		/** output mesh data for Paraview visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) override;
		/** output mesh data for Tecplot visuallization */
		virtual void WriteMeshToPltFile(ofstream &output_file) override;
	};

	/**
	 * @class BaseDataPackage
	 * @brief Abstract base class for a data package 
	 * which is given by a small mesh patch.
	 */
	class BaseDataPackage : public BaseMesh
	{
	public:
		/** lower bound for the addresses as coordinate reference */
		Vecd data_lower_bound_;
		/** number of the addresses for each package data */
		Vecu number_of_addrs_;
		/** the size of data package. */
		size_t pkg_size_;
		/** If true, its data saved in memeory pool. */
		bool is_inner_pkg_;

		/** Constructor with package size information.  */
		BaseDataPackage(size_t pkg_size, size_t buffer_size) 
			: BaseMesh(Vecu(pkg_size), buffer_size), pkg_size_(pkg_size), 
			number_of_addrs_(Vecu(pkg_size + 2 * buffer_size)),
			is_inner_pkg_(false), data_lower_bound_(0) {};
		virtual ~BaseDataPackage() {};

		/** get the size of package*/
		size_t getDataPackageSize() { return pkg_size_; };
		/** initialize package mesh information. */
		void initializePackageGoemetry(Vecd& pkg_lower_bound, Real data_spacing);
		/** initialize the defaultly construted package. */
		virtual void initializeDataPackage(SPHBody* sph_body) = 0;
		/** Bi and tri-linear interpolation within the package. */
		template<class DataType>
		DataType ProbeDataPackage(MeshDataMatrix<DataType*> data, Vecd& position);
	};

	/**
	 * @class LevelSetDataPackage
	 * @brief Fixed memory level data located in a package.
	 */
	class LevelSetDataPackage : public BaseDataPackage
	{
	protected:
		/** Heavside function for computing volume fraction is a cell. */
		Real getHeaviside(Real phi);
		/** Aperture at a grid point for compueting patch area in a cell. */
		Vecd getAperture(Vecu grid_index);
	public:
		/** level set, signed distance. */
		MeshDataMatrix<Real> phi_;
		/** normal direction. */
		MeshDataMatrix<Vecd> n_;
		/** address of the mesh data for spatial operatings.
		  * only valid for inner packages.
		  * The extra layer of data is from the neighbor packages. */
		MeshDataMatrix<Real*> phi_addrs_;
		MeshDataMatrix<Vecd*> n_addrs_;
		/** If true, the packase is near to zero level set. */
		bool is_core_pkg_;

		/** */
		/** defualt constructor */
		LevelSetDataPackage() 
			: BaseDataPackage(4, 1), is_core_pkg_(false) { 
			AllocateMeshDataMatrix(); 
		}
		virtual ~LevelSetDataPackage() { DeleteMeshDataMatrix(); };;

		/** allocate memory for package data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memory for package data */
		virtual void DeleteMeshDataMatrix() override;

		/** initialize the defaultly constructed package. */
		virtual void initializeDataPackage(SPHBody* sph_body) override;
		/** initialize with uniform data*/
		void initializeWithUniformData(Real level_set, Vecd normal_direction);

		/** Get curvation at a grid point.*/
		Real getCurvature(Vecu grid_index);
		/** Get normal direction at a grid point.*/
		Vecd getNormalDirection(Vecu grid_index);
		/** Get area of a surface patch in a cell.*/
		Real getSurfaceInCell(Vecu cell_index);
		/** Get volume fraction in a cell.*/
		Real getVolumeInCell(Vecu cell_index);
		/** Indicate if a grid point near the surafce.*/
		bool isNearSurfaceGrid(Vecu grid_index);
		/** Indicate if a grid point within the narrow band.*/
		bool isNarrowBandGrid(Vecu grid_index);

		/** A Reinitialization step.*/
		bool stepReinitialization(Vecu grid_index);
		/** Reset levelset by computing distance to another levelset.*/
		bool redistanceToLevelset(Vecu grid_index, Real phi);
	};

	/**
	 * @class BaseMeshWithDataPackages
	 * @brief Abstract class for a mesh with data packages
	 */
	class BaseMeshWithDataPackages : public Mesh
	{
	protected:
		/** initialize a cell with package data. */
		virtual void initializeDataInACell(Vecu cell_index, Real dt) = 0;
		/** initialize the addresses in a data package. */
		virtual void initializeAdressesInACell(Vecu cell_index, Real dt) = 0;
		/** tag if a data package is inner package. */
		virtual void tagACellIsInnerPackage(Vecu cell_index, Real dt) = 0;
	public:
		/** Constructor using domain information. */
		BaseMeshWithDataPackages(Vecd lower_bound, 
			Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0)
			: Mesh(lower_bound, upper_bound, grid_spacing, buffer_size) {};
		/** Constructor using mesh information directly. */
		BaseMeshWithDataPackages(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
			: Mesh(mesh_lower_bound, number_of_cells, cell_spacing) {};
		virtual ~BaseMeshWithDataPackages() {};

		/** initialize the data packages with external information. */
		virtual void InitializeDataPackages() = 0;
	};

	/**
	  * @class BaseLevelSet
	  * @brief A abstract desrcibes a mesh with level set data packages.
	  */
	class BaseLevelSet : public BaseMeshWithDataPackages
	{
	protected:
		SPHBody* sph_body_;

	public:
		/** Constructor using domain information. */
		BaseLevelSet(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0);
		/** Constructor using mesh information directly. */
		BaseLevelSet(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		virtual ~BaseLevelSet() {};

		/** set the SPH body externally*/
		void setSPHBody(SPHBody* sph_body) { sph_body_ = sph_body; };
		/** probe the level set at a off-grid position. */
		virtual Real probeLevelSet(Vecd position) = 0;
		/** probe the normal direction at a off-grid position. */
		virtual Vecd probeNormalDirection(Vecd position) = 0;

	};

	/**
	 * @class MeshWithDataPackages
	 * @brief Abstract class fpr mesh with data packages
	 */
	template<class BaseMeshType = BaseLevelSet, class DataPackageType = LevelSetDataPackage>
	class MeshWithDataPackages : public BaseMeshType
	{
	protected:
		/** spacing of data in the data packages*/
		Real data_spacing_;
		/** the size of the data packages*/
		size_t pkg_size_;
		/** toal numer of data in all the packages. */
		Vecu number_of_data_;
		/** singular data packages. prodvied for far field condition. */
		StdVec<DataPackageType*> singular_data_pkgs_addrs;

		/** find the poistion data from its global index */
		Vecd getDataPositionFromIndex(Vecu global_data_index) 
		{
			Vecd data_position;
			for (int n = 0; n < data_position.size(); n++)
			{
				data_position[n] = this->mesh_lower_bound_[n] + 0.5* data_spacing_
					+ Real(global_data_index[n]) * data_spacing_;
			}
			return data_position;
		};
		/** find the value of data from its global index */
		template<class DataType, MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
		DataType getValueFromGlobalDataIndex(Vecu global_data_index);
		/** initialize the addresses in a data package for one varibale. */
		template<class DataType, MeshDataMatrix<DataType*> DataPackageType:: * MemPtrAddrss,
			MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
		void initializeOneVariableAdressesInACell(Vecu cell_index);

	public:
		/** memory pool for all packages in the mesh. */
		Mypool<DataPackageType> data_pakg_pool_;
		/** Address of mesh date packages. */
		MeshDataMatrix<DataPackageType*> data_pkg_addrs_;
		/** Inner data packages which is able to carry out spatial operations. */
		ConcurrentVector<DataPackageType*> inner_data_pkgs_;

		/** allocate memories for addresses of data packages. */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for addresses of data packages. */
		virtual void DeleteMeshDataMatrix() override;

		/** Constructor using domain information. */
		MeshWithDataPackages(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0)
			: BaseMeshType(lower_bound, upper_bound, grid_spacing, buffer_size)
		{
			pkg_size_ = DataPackageType().getDataPackageSize();
			data_spacing_ = grid_spacing / (Real)pkg_size_;
			number_of_data_ = this->number_of_cells_ * pkg_size_;
			AllocateMeshDataMatrix();
		}
		/** Constructor using mesh information directly. */
		MeshWithDataPackages(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
			: BaseMeshType(mesh_lower_bound, number_of_cells, cell_spacing)
		{
			pkg_size_ = DataPackageType().getDataPackageSize();
			data_spacing_ = this->grid_spacing_ / (Real)pkg_size_;
			number_of_data_ = this->number_of_cells_ * pkg_size_;
			AllocateMeshDataMatrix();
		};
		virtual ~MeshWithDataPackages() { DeleteMeshDataMatrix(); };

		/** probe a off-mesh value. */
		template<class DataType, MeshDataMatrix<DataType*> DataPackageType:: * MemPtrAddrss,
			MeshDataMatrix<DataType> DataPackageType:: * MemPtr>
			DataType probeMesh(Vecd& position);
	};
	
	/**
	 * @class LevelSet
	 * @brief Mesh with level set data as packages.
	 */
	class LevelSet
		: public MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>
	{
	protected:
		/** initialize level set in a cell. */
		virtual void initializeDataInACell(Vecu cell_index, Real dt) override;
		/** initialize the addresses in a data package. */
		virtual void initializeAdressesInACell(Vecu cell_index, Real dt) override;
		/** tag if a data package is inner package. */
		virtual void tagACellIsInnerPackage(Vecu cell_index, Real dt) override;
	public:
		/** Core packages which are near to zero level set. */
		ConcurrentVector<LevelSetDataPackage*> core_data_pkgs_;

		/** Constructor using domain information. */
		LevelSet(Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		/** Constructor using mesh information directly. */
		LevelSet(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0);
		/** Constructor using domain and sph body information. */
		LevelSet(SPHBody* sph_body, Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0);
		virtual ~LevelSet() {};

		/** initialize the initial level set field. */
		virtual void InitializeDataPackages() override;
		/** probe the level set at a off-grid position. */
		virtual Real probeLevelSet(Vecd position) override;
		virtual Vecd probeNormalDirection(Vecd position) override;

		/** output mesh data for Paraview visuallization */
		virtual void WriteMeshToVtuFile(ofstream& output_file) override {};
		/** output mesh data for Tecplot visuallization */
		virtual void WriteMeshToPltFile(ofstream& output_file);
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

			/** bulid the zero level mesh first.*/
			int middle_level = ((int)total_levels - 1) / 2;
			Real zero_level_cell_spacing = reference_cell_spacing * powern(2.0, middle_level);
			cell_spacing_levels_.push_back(zero_level_cell_spacing);
			MeshLevelType* zero_level_mesh
				= new MeshLevelType(lower_bound, upper_bound, zero_level_cell_spacing, buffer_size);
			mesh_levels_.push_back(zero_level_mesh);
			Vecu zero_level_number_of_cells = zero_level_mesh->getNumberOfCells();
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

