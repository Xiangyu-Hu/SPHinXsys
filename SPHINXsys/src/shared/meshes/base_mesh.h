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
#include <mutex>
#include <functional>
using namespace std::placeholders;

using namespace std;

namespace SPH
{
	/**
	 * @brief preclaimed classes.
	 */
	class Kernel;

	/** Functor for operation on the mesh. */
	typedef std::function<void(Vecu, Real)> MeshFunctor;
	/** Iterator on the mesh by looping index. sequential computing. */
	void MeshIterator(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);
	/** Iterator on the mesh by looping index. parallel computing. */
	void MeshIterator_parallel(Vecu index_begin, Vecu index_end, MeshFunctor& mesh_functor, Real dt = 0.0);

	/**
	 * @class BaseMesh
	 * @brief Base class for all meshes
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
		BaseMesh(Vecu number_of_grid_points)
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
		 *@param[in] number_of_grid_points(Vecu) number_of_grid_points in each direction.
		 *@param[in] i(size_t) 1D index
		 *@param[out] (Vecu) 2 or 3 D index
		 */
		Vecu transfer1DtoMeshIndex(Vecu number_of_grid_points, size_t i);
		/** convert mesh index to 1d vector index. */
		/**
		 *@brief This function convert mesh index to 1d vector index.
		 *@param[in] number_of_grid_points(Vecu) Mesh size in each direction.
         *@param[in] grid_index Mesh index in each direction.
         *@param[out] (size_t) 1D index.
         */
        size_t transferMeshIndexTo1D(Vecu number_of_grid_points, Vecu grid_index);
        /** generates Morton number for the cells.*/
        /**
         *@brief This function converts mesh index to 1d vector index.
         *@param[in] i Mesh index in x, y or z direction.
         *@param[out] (size_t) Morton Z number.
         */
        size_t MortonCode(const size_t &i);
        /** converts mesh index into a Morton order.
         * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
         * https://stackoverflow.com/questions/18529057/
         * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
         */
        /**
         *@brief This function converts mesh index into a Morton order.
         *@param[in] grid_index Mesh index in each direction.
         *@param[out] (size_t) 1D index.
         */
        size_t transferMeshIndexToMortonOrder(Vecu grid_index);
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
		size_t buffer_width_;	/**< buffer width to avoid bound check.*/
		Real cell_spacing_; 	/**< cell_spacing */
		Vecu number_of_cells_; 	/**< number of cells by dimension */
		/**
		 *@brief This function set the mesh lower bound including the buffer region.
		 *@param[in] lower_bound(Vecd) Input mesh lower bound.
		 *@param[in] grid_spacing(Real) Input grid spacing.
		 *@param[in] buffer_width(Real)  Buffersize to extened the mesh from the physical bound of a body.
		 */
		void setMeshLowerBound(Vecd lower_bound, Real grid_spacing, size_t buffer_width);
		/**
		 *@brief This function compute number of total cells
		 *@param[in] lower_bound(Vecd) Input mesh lower bound.
		 *@param[in] grid_spacing(Real) Input grid spacing.
		 *@param[in] buffer_width(Real)  Buffersize to extend the mesh from physical domain of a body or something.
		 */
		Vecu calcNumberOfCells(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_width);
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
			size_t buffer_width = 0 /**< Buffer size. */
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
		size_t MeshBufferSize() { return buffer_width_; };
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
		virtual void writeMeshToVtuFile(ofstream& output_file) {};
		/** output mesh data for Tecplot visualization */
		virtual void writeMeshToPltFile(ofstream& output_file) {};

		/** allcate memories for the mesh data matrix*/
		virtual void allocateMeshDataMatrix() {};
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() {};
	};

	/**
	 * @class MultilevelMesh
	 * @brief Multi level Meshes with multi resolution mesh data
	 */
	template<class BaseMeshType, class MeshLevelType>
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
			Real reference_cell_spacing, size_t total_levels, size_t buffer_width = 0)
			: BaseMeshType(lower_bound, upper_bound, reference_cell_spacing, buffer_width),
			total_levels_(total_levels) {

			/** build the zero level mesh first.*/
			int middle_level = ((int)total_levels - 1) / 2;
			Real zero_level_cell_spacing = reference_cell_spacing * powern(2.0, middle_level);
			cell_spacing_levels_.push_back(zero_level_cell_spacing);
			MeshLevelType* zero_level_mesh
				= new MeshLevelType(lower_bound, upper_bound, zero_level_cell_spacing, buffer_width);
			mesh_levels_.push_back(zero_level_mesh);
			Vecu zero_level_number_of_cells = zero_level_mesh->NumberOfCells();
			/** copy zero level mesh perperties to this. */
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
}
