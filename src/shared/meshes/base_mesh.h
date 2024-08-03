/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
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
#include "my_memory_pool.h"
#include "sph_data_containers.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <mutex>
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
class Mesh
{
  protected:
    Vecd mesh_lower_bound_;  /**< mesh lower bound as reference coordinate */
    Real grid_spacing_;      /**< grid_spacing */
    Arrayi all_grid_points_; /**< number of grid points by dimension */
  public:
    Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
    ~Mesh(){};

    Vecd MeshLowerBound() { return mesh_lower_bound_; };
    Real GridSpacing() { return grid_spacing_; };
    Arrayi AllGridPoints() { return all_grid_points_; };
    Arrayi AllCells() { return all_grid_points_ - Arrayi::Ones(); };
    size_t NumberOfGridPoints() { return transferMeshIndexTo1D(all_grid_points_, all_grid_points_); };
    size_t NumberOfCells() { return transferMeshIndexTo1D(AllCells(), AllCells()); };
    /** Given the grid point number, return the cell number. */
    Arrayi AllCellsFromAllGridPoints(const Arrayi &all_grid_points) { return all_grid_points - Arrayi::Ones(); };
    /** Given the cell position, return the grid position. */
    Vecd GridPositionFromCellPosition(const Vecd &cell_position) { return cell_position - 0.5 * grid_spacing_ * Vecd::Ones(); };
    /** Given the position, return the cell index. */
    Arrayi CellIndexFromPosition(const Vecd &position);
    /** Given the cell index, return the cell position. */
    Vecd CellPositionFromIndex(const Arrayi &cell_index);
    /** Given the cell index, return the position of its lower corner. */
    Vecd CellLowerCorner(const Arrayi &cell_index);
    /** Given the index, return the grid position. */
    Vecd GridPositionFromIndex(const Arrayi &grid_index);
    /** Transfer 1D int to mesh index.  */
    Arrayi transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i);
    /** Transfer mesh index to 1D int.  */
    size_t transferMeshIndexTo1D(const Arrayi &mesh_size, const Arrayi &mesh_index);
    /** converts mesh index into a Morton order.
     * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
     * https://stackoverflow.com/questions/18529057/
     * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
     */
    size_t MortonCode(const size_t &i);
    /** Converts mesh index into a Morton order. */
    size_t transferMeshIndexToMortonOrder(const Arrayi &mesh_index);
};

/**
 * @class BaseMeshField
 * @brief Abstract base class for the geometric or physics field.
 */
class BaseMeshField
{
  protected:
    std::string name_{};

  public:
    explicit BaseMeshField(const std::string &name) : name_(name){};
    virtual ~BaseMeshField(){};
    /** Return the mesh field name. */
    std::string Name() { return name_; };
    /** output mesh data for Tecplot visualization */
    virtual void writeMeshFieldToPlt(std::ofstream &output_file) = 0;
};

/**
 * @class 	RefinedMesh
 * @brief 	Abstract base class derived from the coarse mesh but has double resolution.
 */
template <class CoarseMeshType>
class RefinedMesh;

/**
 * @class 	MultilevelMesh
 * @brief 	Multi-level Meshes with successively double the resolution
 */
template <class MeshFieldType, class CoarsestMeshType>
class MultilevelMesh : public MeshFieldType
{
  public:
    /**template parameter pack is used with rvalue reference and perfect forwarding to keep
     * the type of arguments when called by another function with template parameter pack too. */
    template <typename... Args>
    MultilevelMesh(BoundingBox tentative_bounds, Real reference_spacing, size_t total_levels, Args &&...args)
        : MeshFieldType(std::forward<Args>(args)...), total_levels_(total_levels)
    {
        mesh_levels_.push_back(
            mesh_level_ptr_vector_keeper_
                .template createPtr<CoarsestMeshType>(tentative_bounds, reference_spacing, std::forward<Args>(args)...));

        for (size_t level = 1; level != total_levels_; ++level)
        {
            /** all mesh levels aligned at the lower bound of tentative_bounds */
            mesh_levels_.push_back(
                mesh_level_ptr_vector_keeper_
                    .template createPtr<RefinedMesh<CoarsestMeshType>>(tentative_bounds, *mesh_levels_.back(), std::forward<Args>(args)...));
        }
    };
    virtual ~MultilevelMesh(){};

  private:
    UniquePtrsKeeper<CoarsestMeshType> mesh_level_ptr_vector_keeper_;

  protected:
    size_t total_levels_;                    /**< level 0 is the coarsest */
    StdVec<CoarsestMeshType *> mesh_levels_; /**< Mesh in different coarse level. */

  public:
    /** Return the mesh at different level. */
    StdVec<CoarsestMeshType *> getMeshLevels() { return mesh_levels_; };
    /** Write mesh data to file. */
    void writeMeshFieldToPlt(std::ofstream &output_file) override
    {
        for (size_t l = 0; l != total_levels_; ++l)
        {
            mesh_levels_[l]->writeMeshFieldToPlt(output_file);
        }
    }
};
} // namespace SPH
#endif // BASE_MESH_H