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


#ifndef BASE_MESH_H
#define BASE_MESH_H



#include "base_data_package.h"
#include "sph_data_conainers.h"
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
		BaseMesh();
		BaseMesh(Vecu number_of_grid_points);
		virtual ~BaseMesh() {};

		Vecd MeshLowerBound() { return mesh_lower_bound_; };
		Real GridSpacing() { return grid_spacing_; };
		Vecu NumberOfGridPoints() { return number_of_grid_points_; };

		Vecu GridIndexFromPosition(const Vecd& position);
		Vecd GridPositionFromIndex(Vecu grid_index);
		Vecu transfer1DtoMeshIndex(Vecu number_of_grid_points, size_t i);
        size_t transferMeshIndexTo1D(Vecu number_of_grid_points, Vecu grid_index);
        /** converts mesh index into a Morton order.
         * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
         * https://stackoverflow.com/questions/18529057/
         * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
         */
        size_t MortonCode(const size_t &i);
        /** This function converts mesh index into a Morton order. */
        size_t transferMeshIndexToMortonOrder(Vecu grid_index);
    };

	/**
	 * @class Mesh
	 * @brief Abstract base class for cell-based mesh properties.
	 * The mesh is proposed for several functions.
	 * First, it is used in cell linked list for neighbor search.
	 * Second, it is used for background maps such as level sets.
	 * This class is the counterpart of the class particles.
	 */
	class Mesh : public BaseMesh
	{
	protected:
		std::string name_;
		size_t buffer_width_;	/**< buffer width to avoid bound check.*/
		Vecu number_of_cells_; 	/**< number of cells by dimension */
		Vecu NumberOfGridPointsFromNumberOfCells(Vecu number_of_cells) { return number_of_cells + Vecu(1); };
		Vecu NumberOfCellsFromNumberOfGridPoints(Vecu number_of_grid_points) { return number_of_grid_points - Vecu(1); };
		void copyMeshProperties(Mesh* another_mesh);
		Vecd GridPositionFromCellPosition(Vecd& cell_position)
		{
			return cell_position - Vecd(0.5 * grid_spacing_);
		};
		void initializeWithBoundingBox(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
	public:
		Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
		Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real grid_spacing);
		virtual ~Mesh() {};

		std::string Name() { return name_; };
		Vecu NumberOfCells() { return number_of_cells_; };
		size_t MeshBufferSize() { return buffer_width_; };
		/** This function check whether a position well within in the mesh bounds */
		bool isWithinMeshBound(Vecd position);
		Vecu CellIndexFromPosition(Vecd& position);
		Vecd CellPositionFromIndex(Vecu cell_index);

		/** output mesh data for Tecplot visualization */
		virtual void writeMeshToPltFile(std::ofstream& output_file) {};

		/** allocate memories for the mesh data matrix*/
		virtual void allocateMeshDataMatrix() {};
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() {};
	};

	/**
	 * @class MultilevelMesh
	 * @brief Multi level Meshes with successively double the resolutions
	 */
	template<class MeshCompositionType, class BaseMeshType, class MeshLevelType>
	class MultilevelMesh : public BaseMeshType
	{
	protected:
		size_t total_levels_; /**< level 0 is the coarsest */
		StdVec<MeshLevelType*> mesh_levels_;
	public:
		MultilevelMesh(MeshCompositionType& mesh_composition, ParticleAdaptation& particle_adaptation, 
			BoundingBox tentative_bounds, Real reference_spacing, 
			size_t total_levels, Real maximum_spacing_ratio, size_t buffer_width) : 
			BaseMeshType(mesh_composition, particle_adaptation, tentative_bounds, 
			reference_spacing, buffer_width), total_levels_(total_levels) 
		{
			Real zero_level_spacing = reference_spacing * maximum_spacing_ratio;
			for (size_t level = 0; level != total_levels_; ++level) {
				Real spacing_level = zero_level_spacing * powerN(0.5, (int)level);
				/** all mesh levels aligned at the lower bound of tentative_bounds */
				MeshLevelType* mesh_level = 
					new MeshLevelType(mesh_composition, particle_adaptation, tentative_bounds, spacing_level, buffer_width);
				mesh_levels_.push_back(mesh_level);
			}
		};

		virtual ~MultilevelMesh() 
		{ 
			for (size_t l = 0; l != total_levels_; ++l) mesh_levels_[l]->~MeshLevelType();
		};

		StdVec<MeshLevelType*> getMeshLevels() { return mesh_levels_; };

		void writeMeshToPltFile(std::ofstream& output_file) override
		{
			for (size_t l = 0; l != total_levels_; ++l) {
				mesh_levels_[l]->writeMeshToPltFile(output_file);
			}
		}
	};
}
#endif //BASE_MESH_H