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
* @brief 	This is the base classes of mesh, which is used to describe ordered and indexed
*			data sets.  Depending on application, there are different data 
* 			saved on the mesh. The intersection points of mesh lines are called 
*			grid points, the element enclosed by mesh lines (2D) or faces (3D) called 
*			cells. The mesh line or face are also called cell faces. Grid points are
*			also called cell corners.
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef BASE_MESH_H
#define BASE_MESH_H

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "my_memory_pool.h"

#include <fstream>
#include <algorithm>
#include <mutex>
#include <functional>
using namespace std::placeholders;

namespace SPH
{
	class Kernel;
	class ParticleAdaptation;

	/** Functor for operation on the mesh. */
	typedef std::function<void(const Vecu &, Real)> MeshFunctor;
	/** Iterator on the mesh by looping index. sequential computing. */
	void MeshIterator(const Vecu &index_begin, const Vecu &index_end, MeshFunctor &mesh_functor, Real dt = 0.0);
	/** Iterator on the mesh by looping index. parallel computing. */
	void MeshIterator_parallel(const Vecu &index_begin, const Vecu &index_end, MeshFunctor &mesh_functor, Real dt = 0.0);

	/**
	 * @class BaseMesh
	 * @brief Base class for all meshes which may be grid or cell based.
	 * The basic properties of the mesh, such as lower bound, grid spacing 
	 * and number of grid points may be determined by the derived class.
	 * Note that there is no mesh-based data defined here.
	 */
	class BaseMesh
	{
	protected:
		Vecd mesh_lower_bound_;		 /**< mesh lower bound as reference coordinate */
		Real grid_spacing_;			 /**< grid_spacing */
		Vecu number_of_grid_points_; /**< number of grid points by dimension */
	public:
		BaseMesh();
		BaseMesh(Vecu number_of_grid_points);
		BaseMesh(Vecd mesh_lower_bound, Real grid_spacing, Vecu number_of_grid_points);
		BaseMesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
		virtual ~BaseMesh(){};

		Vecd MeshLowerBound() { return mesh_lower_bound_; };
		Real GridSpacing() { return grid_spacing_; };
		Vecu NumberOfGridPoints() { return number_of_grid_points_; };
		Vecu NumberOfGridPointsFromNumberOfCells(const Vecu &number_of_cells) { return number_of_cells + Vecu(1); };
		Vecu NumberOfCellsFromNumberOfGridPoints(const Vecu &number_of_grid_points) { return number_of_grid_points - Vecu(1); };
		Vecd GridPositionFromCellPosition(const Vecd &cell_position) { return cell_position - Vecd(0.5 * grid_spacing_); };

		Vecu CellIndexFromPosition(const Vecd &position);
		Vecd CellPositionFromIndex(const Vecu &cell_index);
		/** Note that, the lower corner grid of a cell has the same index as the cell. */
		Vecd GridPositionFromIndex(const Vecu &grid_index);
		Vecu transfer1DtoMeshIndex(const Vecu &number_of_grid_points, size_t i);
		size_t transferMeshIndexTo1D(const Vecu &number_of_grid_points, const Vecu &grid_index);
		/** converts mesh index into a Morton order.
         * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
         * https://stackoverflow.com/questions/18529057/
         * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
         */
		size_t MortonCode(const size_t &i);
		/** This function converts mesh index into a Morton order. */
		size_t transferMeshIndexToMortonOrder(const Vecu &grid_index);
	};

	/**
	 * @class Mesh
	 * @brief Abstract base class for cell-based mesh 
	 * by introducing number of cells, buffer width and mesh-based data in its derived classes.
	 */
	class Mesh : public BaseMesh
	{
	protected:
		size_t buffer_width_;  /**< buffer width to avoid bound check.*/
		Vecu number_of_cells_; /**< number of cells by dimension */

		void copyMeshProperties(Mesh *another_mesh);
		/** allocate memories for the mesh data matrix*/
		virtual void allocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() = 0;

	public:
		Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
		Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real grid_spacing);
		virtual ~Mesh(){};

		Vecu NumberOfCells() { return number_of_cells_; };
		size_t MeshBufferSize() { return buffer_width_; };
	};

	/**
	 * @class BaseMeshField
	 * @brief Abstract base class for the field data saved on a mesh.
	 */
	class BaseMeshField
	{
	protected:
		std::string name_;

	public:
		BaseMeshField(std::string name) : name_(name){};
		virtual ~BaseMeshField(){};

		std::string Name() { return name_; };
		/** output mesh data for Tecplot visualization */
		virtual void writeMeshFieldToPlt(std::ofstream &output_file) = 0;
	};

	/**
	 * @class MultilevelMesh
	 * @brief Multi level Meshes with successively double the resolutions
	 */
	template <class MeshFieldType, class MeshLevelType>
	class MultilevelMesh : public MeshFieldType
	{
	protected:
		size_t total_levels_; /**< level 0 is the coarsest */
		StdVec<MeshLevelType *> mesh_levels_;

	public:
		/**template parameter pack is used with rvalue reference and perfect forwarding to keep 
		 * the type of arguments when called by another function with template parameter pack too. */
		template <typename... Args>
		MultilevelMesh(BoundingBox tentative_bounds, Real reference_spacing, size_t total_levels,
					   Real maximum_spacing_ratio, Args &&...args)
			: MeshFieldType(std::forward<Args>(args)...), total_levels_(total_levels)
		{
			Real zero_level_spacing = reference_spacing * maximum_spacing_ratio;
			for (size_t level = 0; level != total_levels_; ++level)
			{
				Real spacing_level = zero_level_spacing * powerN(0.5, (int)level);
				/** all mesh levels aligned at the lower bound of tentative_bounds */
				MeshLevelType *mesh_level =
					new MeshLevelType(tentative_bounds, spacing_level, std::forward<Args>(args)...);
				mesh_levels_.push_back(mesh_level);
			}
		};

		virtual ~MultilevelMesh()
		{
			for (size_t l = 0; l != total_levels_; ++l)
				mesh_levels_[l]->~MeshLevelType();
		};

		StdVec<MeshLevelType *> getMeshLevels() { return mesh_levels_; };

		void writeMeshFieldToPlt(std::ofstream &output_file) override
		{
			for (size_t l = 0; l != total_levels_; ++l)
			{
				mesh_levels_[l]->writeMeshFieldToPlt(output_file);
			}
		}
	};
}
#endif //BASE_MESH_H