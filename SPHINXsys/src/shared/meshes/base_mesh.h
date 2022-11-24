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
 * @file 	base_mesh.h
 * @brief 	This is the base classes of mesh, which is used to describe ordered and indexed
 *			data sets.  Depending on application, there are different data
 * 			saved on the mesh. The intersection points of mesh lines are called
 *			grid points, the element enclosed by mesh lines (2D) or faces (3D) called
 *			cells. The mesh line or face are also called cell faces. Grid points are
 *			also called cell corners.
 * @author	Chi Zhang and Xiangyu Hu
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
	/**
	 * @class BaseMesh
	 * @brief Base class for all structured meshes which may be grid or cell based.
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
		explicit BaseMesh(Vecu number_of_grid_points);
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

		/** iteration with void (non_value_returning) function. 2D case. */
		template <int lower0, int upper0,
				  int lower1, int upper1, typename FunctionOnEach>
		inline void for_each2d(const FunctionOnEach &function)
		{
			for (int l = lower0; l != upper0; ++l)
				for (int m = lower1; m != upper1; ++m)
				{
					function(l, m);
				}
		};
		template <int lower, int upper, typename FunctionOnEach>
		inline void for_each2d(const FunctionOnEach &function)
		{
			for_each2d<lower, upper, lower, upper, FunctionOnEach>(function);
		};

		/** iteration with void (non_value_returning) function. 2D case. */
		template <int lower0, int upper0,
				  int lower1, int upper1, typename CheckOnEach>
		inline Vec2i find_if2d(const CheckOnEach &function)
		{
			for (int l = lower0; l != upper0; ++l)
				for (int m = lower1; m != upper1; ++m)
				{
					if (function(l, m))
						return Vec2i(l, m);
				}
			return Vec2i(upper0, upper1);
		};
		template <int lower, int upper, typename CheckOnEach>
		inline Vec2i find_if2d(const CheckOnEach &function)
		{
			return for_each2d<lower, upper, lower, upper, CheckOnEach>(function);
		};

		/** iteration with void (non_value_returning) function. 3D case. */
		template <int lower0, int upper0,
				  int lower1, int upper1,
				  int lower2, int upper2, typename FunctionOnEach>
		inline void for_each3d(const FunctionOnEach &function)
		{
			for (int l = lower0; l != upper0; ++l)
				for (int m = lower1; m != upper1; ++m)
					for (int n = lower2; n != upper2; ++n)
					{
						function(l, m, n);
					}
		};
		template <int lower, int upper, typename FunctionOnEach>
		inline void for_each3d(const FunctionOnEach &function)
		{
			for_each3d<lower, upper, lower, upper, lower, upper, FunctionOnEach>(function);
		};

		/** iteration with void (non_value_returning) function. 3D case. */
		template <int lower0, int upper0,
				  int lower1, int upper1,
				  int lower2, int upper2, typename CheckOnEach>
		inline Vec3i find_if3d(const CheckOnEach &function)
		{
			for (int l = lower0; l != upper0; ++l)
				for (int m = lower1; m != upper1; ++m)
					for (int n = lower2; n != upper2; ++n)
					{
						if (function(l, m, n))
							return Vec3i(l, m, n);
					}
			return Vec3i(upper0, upper1, upper2);
		};
		template <int lower, int upper, typename CheckOnEach>
		inline Vec3i find_if3d(const CheckOnEach &function)
		{
			return find_if3d<lower, upper, lower, upper, lower, upper, CheckOnEach>(function);
		};

		Vecu transfer1DtoMeshIndex(const Vecu &number_of_mesh_indexes, size_t i);
		size_t transferMeshIndexTo1D(const Vecu &number_of_mesh_indexes, const Vecu &mesh_index);
		/** converts mesh index into a Morton order.
		 * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
		 * https://stackoverflow.com/questions/18529057/
		 * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
		 */
		size_t MortonCode(const size_t &i);
		/** This function converts mesh index into a Morton order. */
		size_t transferMeshIndexToMortonOrder(const Vecu &mesh_index);
	};

	/**
	 * @class Mesh
	 * @brief Abstract base class for cell-based mesh
	 * by introducing number of cells, buffer width and mesh-based data in its derived classes.
	 * Note that we identify the difference between grid spacing and data spacing.
	 * The latter is different from grid spacing when MeshWithDataPackage is considered.
	 */
	class Mesh : public BaseMesh
	{
	protected:
		size_t buffer_width_;  /**< buffer width to avoid bound check.*/
		Vecu number_of_cells_; /**< number of cells by dimension */

		void copyMeshProperties(Mesh *another_mesh);

	public:
		Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
		Mesh(Vecd mesh_lower_bound, Vecu number_of_cells, Real grid_spacing);
		virtual ~Mesh(){};

		Vecu NumberOfCells() { return number_of_cells_; };
		size_t MeshBufferSize() { return buffer_width_; };
		virtual Real DataSpacing() { return grid_spacing_; };
	};

	/**
	 * @class BaseMeshField
	 * @brief Abstract base class for the geometric or physics field.
	 */
	class BaseMeshField
	{
	protected:
		std::string name_;

	public:
		explicit BaseMeshField(const std::string &name) : name_(name){};
		virtual ~BaseMeshField(){};

		std::string Name() { return name_; };
		/** output mesh data for Tecplot visualization */
		virtual void writeMeshFieldToPlt(std::ofstream &output_file) = 0;
	};

	/**
	 * @class RefinedMesh
	 * @brief Abstract base class derived from the coarse mesh but has double resolution.
	 * Currently, the design is simple but can be extending for more inter-mesh operations.
	 */
	template <class CoarseMeshType>
	class RefinedMesh : public CoarseMeshType
	{
	public:
		template <typename... Args>
		RefinedMesh(BoundingBox tentative_bounds, CoarseMeshType &coarse_mesh, Args &&...args)
			: CoarseMeshType(tentative_bounds, 0.5 * coarse_mesh.DataSpacing(), std::forward<Args>(args)...),
			  coarse_mesh_(coarse_mesh){};
		virtual ~RefinedMesh(){};

	protected:
		CoarseMeshType &coarse_mesh_;
	};

	/**
	 * @class MultilevelMesh
	 * @brief Multi-level Meshes with successively double the resolution
	 */
	template <class MeshFieldType, class CoarsestMeshType, class RefinedMeshType>
	class MultilevelMesh : public MeshFieldType
	{
	private:
		UniquePtrKeepers<CoarsestMeshType> mesh_level_ptr_vector_keeper_;

	protected:
		size_t total_levels_; /**< level 0 is the coarsest */
		StdVec<CoarsestMeshType *> mesh_levels_;

	public:
		/**template parameter pack is used with rvalue reference and perfect forwarding to keep
		 * the type of arguments when called by another function with template parameter pack too. */
		template <typename... Args>
		MultilevelMesh(BoundingBox tentative_bounds, Real reference_spacing, size_t total_levels, Args &&...args)
			: MeshFieldType(std::forward<Args>(args)...), total_levels_(total_levels)
		{
			Real coarsest_spacing = reference_spacing;
			mesh_levels_.push_back(
				mesh_level_ptr_vector_keeper_
					.template createPtr<CoarsestMeshType>(tentative_bounds, coarsest_spacing, std::forward<Args>(args)...));

			for (size_t level = 1; level != total_levels_; ++level)
			{
				/** all mesh levels aligned at the lower bound of tentative_bounds */
				mesh_levels_.push_back(
					mesh_level_ptr_vector_keeper_
						.template createPtr<RefinedMeshType>(tentative_bounds, *mesh_levels_.back(), std::forward<Args>(args)...));
			}
		};

		virtual ~MultilevelMesh(){};

		StdVec<CoarsestMeshType *> getMeshLevels() { return mesh_levels_; };

		void writeMeshFieldToPlt(std::ofstream &output_file) override
		{
			for (size_t l = 0; l != total_levels_; ++l)
			{
				mesh_levels_[l]->writeMeshFieldToPlt(output_file);
			}
		}
	};
}
#endif // BASE_MESH_H