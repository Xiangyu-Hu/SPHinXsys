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
#include "sphinxsys_containers.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <mutex>
using namespace std::placeholders;

namespace SPH
{
/**
 * @class Mesh
 * @brief Base class for all structured meshes which may be grid or cell based.
 * The basic properties of the mesh, such as lower bound, grid spacing
 * and number of grid points may be determined by the derived class.
 * Note that there is no mesh-based data defined here.
 */
class Mesh
{
  public:
    Mesh(BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width);
    Mesh(Vecd mesh_lower_bound, Real grid_spacing, Arrayi all_grid_points);
    ~Mesh(){};

    Vecd MeshLowerBound() const { return mesh_lower_bound_; };
    Real GridSpacing() const { return grid_spacing_; };
    Arrayi AllGridPoints() const { return all_grid_points_; };
    Arrayi AllCells() const { return all_cells_; };
    size_t NumberOfGridPoints() const { return transferMeshIndexTo1D(all_grid_points_, all_grid_points_); };
    size_t NumberOfCells() const { return transferMeshIndexTo1D(all_cells_, all_cells_); };

    Arrayi CellIndexFromPosition(const Vecd &position) const
    {
        return floor((position - mesh_lower_bound_).array() / grid_spacing_)
            .cast<int>()
            .max(Arrayi::Zero())
            .min(all_grid_points_ - 2 * Arrayi::Ones());
    };

    size_t LinearCellIndexFromPosition(const Vecd &position) const
    {
        return transferMeshIndexTo1D(all_cells_, CellIndexFromPosition(position));
    };

    size_t LinearCellIndexFromCellIndex(const Arrayi &cell_index) const
    {
        return transferMeshIndexTo1D(all_cells_, cell_index);
    };

    Vecd CellPositionFromIndex(const Arrayi &cell_index) const;
    Vecd GridPositionFromIndex(const Arrayi &grid_index) const;
    Vecd CellLowerCornerPosition(const Arrayi &cell_index) const;
    //----------------------------------------------------------------------
    // Transferring between 1D mesh indexes.
    // Here, mesh size can be either AllGridPoints or AllCells.
    //----------------------------------------------------------------------
    Arrayi transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i) const;

    size_t transferMeshIndexTo1D(const Array2i &mesh_size, const Array2i &mesh_index) const
    {
        return mesh_index[0] * mesh_size[1] + mesh_index[1];
    };

    size_t transferMeshIndexTo1D(const Array3i &mesh_size, const Array3i &mesh_index) const
    {
        return mesh_index[0] * mesh_size[1] * mesh_size[2] +
               mesh_index[1] * mesh_size[2] +
               mesh_index[2];
    };
    /** converts mesh index into a Morton order.
     * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
     * https://stackoverflow.com/questions/18529057/
     * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
     */
    size_t transferMeshIndexToMortonOrder(const Array2i &mesh_index) const
    {
        return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
    };

    size_t transferMeshIndexToMortonOrder(const Array3i &mesh_index) const
    {
        return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1) | (MortonCode(mesh_index[2]) << 2);
    };

  protected:
    Vecd mesh_lower_bound_;  /**< mesh lower bound as reference coordinate */
    Real grid_spacing_;      /**< grid_spacing */
    size_t buffer_width_;    /**< buffer width to avoid bound check.*/
    Arrayi all_grid_points_; /**< number of grid points by dimension */
    Arrayi all_cells_;       /**< number of cells by dimension */

    size_t MortonCode(const size_t &i) const
    {
        size_t x = i;
        x &= 0x3ff;
        x = (x | x << 16) & 0x30000ff;
        x = (x | x << 8) & 0x300f00f;
        x = (x | x << 4) & 0x30c30c3;
        x = (x | x << 2) & 0x9249249;
        return x;
    };
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