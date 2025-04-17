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
 * @file    mesh_with_data_packages.h
 * @brief   This class is designed to save memory and increase computational efficiency on mesh.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_H
#define MESH_WITH_DATA_PACKAGES_H

#include "base_mesh.h"
#include "my_memory_pool.h"
#include "sphinxsys_variable.h"
#include "tbb/parallel_sort.h"
#include "mesh_iterators.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <mutex>
using namespace std::placeholders;

namespace SPH
{
/**
 * @class MeshWithGridDataPackages
 * @brief Abstract class for mesh with grid-based data packages.
 * @details The idea is to save sparse data on a cell-based mesh.
 * We say sparse data, it means that only in some inner mesh cells there are no trivial data.
 * A typical example is a level-set field which only has meaningful values near the interface,
 * while the latter is in the inner region of a mesh.
 * In this class, only some inner mesh cells are filled with data packages.
 * Each data package is again a mesh, but grid based, where two sets of data are saved on its grid points.
 * One is the field data of matrices with pkg_size, the other is corresponding address data of matrices with pkg_addrs_size.
 * For two neighboring data packages, they share the data in the buffer which is in the overlap region.
 * The filling of field data is achieved first by the data matrices by the function initializeDataInACell
 * and then the address matrix by the function initializeAddressesInACell.
 * All these data packages are indexed by a concurrent vector inner_data_pkgs_.
 * Note that a data package should be not near the mesh bound, otherwise one will encounter the error "out of range".
 */
template <int PKG_SIZE>
class MeshWithGridDataPackages : public Mesh
{
  private:
    DataContainerUniquePtrAssemble<MeshVariable> mesh_variable_ptrs_;

  public:
    template <typename... Args>
    explicit MeshWithGridDataPackages(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size)
        : Mesh(tentative_bounds, pkg_size * data_spacing, buffer_size),
          global_mesh_(mesh_lower_bound_ + 0.5 * data_spacing * Vecd::Ones(), data_spacing, all_cells_ * pkg_size),
          data_spacing_(data_spacing)
    {
        allocateIndexDataMatrix();
    };
    virtual ~MeshWithGridDataPackages()
    {
        deleteIndexDataMatrix();
        delete[] cell_neighborhood_;
        delete[] meta_data_cell_;
    };
    /** spacing between the data, which is 1/ pkg_size of this grid spacing */
    Real DataSpacing() { return data_spacing_; };
    Real GridSpacing() { return grid_spacing_; };
    size_t BufferWidth() { return buffer_width_; };

  public:
    ConcurrentVec<std::pair<size_t, int>> occupied_data_pkgs_; /**< (size_t)sort_index, (int)core1/inner0. */
    CellNeighborhood *cell_neighborhood_;                  /**< 3*3(*3) array to store indicies of neighborhood cells. */
    std::pair<Arrayi, int> *meta_data_cell_;          /**< metadata for each occupied cell: (arrayi)cell index, (int)core1/inner0. */
    Mesh global_mesh_;                            /**< the mesh for the locations of all possible data points. */
    size_t num_grid_pkgs_ = 2;                        /**< the number of all distinct packages, initially only 2 singular packages. */

  protected:
    MeshVariableAssemble all_mesh_variables_;         /**< all mesh variables on this mesh. */
    static constexpr int pkg_size = PKG_SIZE;         /**< the size of the data package matrix*/
    const Real data_spacing_;                         /**< spacing of data in the data packages*/
    using MetaData = std::pair<int, size_t>;          /**< stores the metadata for each cell: (int)singular0/inner1/core2, (size_t)package data index*/
    MeshDataMatrix<size_t> index_data_mesh_;         /**< metadata for all cells. */
    using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
    template <typename DataType>
    using PackageData = PackageDataMatrix<DataType, pkg_size>;
    /** Matrix data for temporary usage. */
    template <typename DataType>
    using PackageTemporaryData = PackageDataMatrix<DataType, pkg_size + 1>;

    void allocateIndexDataMatrix(); /**< allocate memories for metadata of data packages. */
    void deleteIndexDataMatrix();   /**< delete memories for metadata of data packages. */

    /** resize all mesh variable data field with `num_grid_pkgs_` size(initially only singular data) */
    struct ResizeMeshVariableData
    {
        template <typename DataType>
        void operator()(DataContainerAddressKeeper<MeshVariable<DataType>> &all_mesh_variables_,
                        const size_t num_grid_pkgs_)
        {
            for (size_t l = 0; l != all_mesh_variables_.size(); ++l)
            {
                MeshVariable<DataType> *variable = all_mesh_variables_[l];
                variable->allocateAllMeshVariableData(num_grid_pkgs_);
            }
        }
    };
    OperationOnDataAssemble<MeshVariableAssemble, ResizeMeshVariableData> resize_mesh_variable_data_{};

    /** probe by applying bi and tri-linear interpolation within the package. */
    template <class DataType>
    DataType probeDataPackage(MeshVariable<DataType> &mesh_variable, size_t package_index, const Arrayi &cell_index, const Vecd &position);
    bool isSingularDataPackage(const Arrayi &cell_index);

    /** return the position of the lower bound data in a cell. */
    Vecd DataLowerBoundInCell(const Arrayi &cell_index)
    {
        return CellLowerCornerPosition(cell_index) + 0.5 * data_spacing_ * Vecd::Ones();
    }

    /** return the grid index from its position and the index of the cell it belongs to. */
    Arrayi DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position)
    {
        return floor((position - DataLowerBoundInCell(cell_index)).array() / data_spacing_)
            .template cast<int>()
            .max(Arrayi::Zero())
            .min((pkg_size - 1) * Arrayi::Ones());
    }

  public:
    std::pair<size_t, Arrayi> NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour);
    bool isInnerDataPackage(const Arrayi &cell_index);
    int DataPackageSize() const { return pkg_size; }
    void assignDataPackageIndex(const Arrayi &cell_index, const size_t package_index);
    size_t PackageIndexFromCellIndex(const Arrayi &cell_index);
    /** obtain averaged value at a corner of a data cell */
    template <typename DataType>
    DataType CornerAverage(MeshVariable<DataType> &mesh_variable,
                           Arrayi addrs_index, Arrayi corner_direction,
                           CellNeighborhood &neighborhood);
    /** This function find the value of data from its index from global mesh. */
    template <typename DataType>
    DataType DataValueFromGlobalIndex(MeshVariable<DataType> &mesh_variable,
                                      const Arrayi &global_grid_index);
    /** assign value to data package according to the position of data */
    template <typename DataType, typename FunctionByPosition>
    void assignByPosition(MeshVariable<DataType> &mesh_variable,
                          const Arrayi &cell_index,
                          const FunctionByPosition &function_by_position);
    /** This function probe a mesh value */
    template <class DataType>
    DataType probeMesh(MeshVariable<DataType> &mesh_variable, const Vecd &position);
    /** return the position of data from its local grid index and the index of the cell it belongs to. */
    Vecd DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index)
    {
        return DataLowerBoundInCell(cell_index) + data_index.cast<Real>().matrix() * data_spacing_;
    }
    /** void (non_value_returning) function iterate on all data points by value. */
    template <typename FunctionOnData>
    void for_each_cell_data(const FunctionOnData &function);
    void resizeMeshVariableData()
    {
        resize_mesh_variable_data_(all_mesh_variables_, num_grid_pkgs_);
    }

    template <typename DataType>
    MeshVariable<DataType> *getMeshVariable(const std::string &variable_name)
    {
        return findVariableByName<DataType>(all_mesh_variables_, variable_name);
    }

    /** compute gradient transform within data package at `package_index` */
    template <typename InDataType, typename OutDataType>
    void computeGradient(MeshVariable<InDataType> &in_variable,
                         MeshVariable<OutDataType> &out_variable,
                         const size_t package_index);

    void registerOccupied(std::pair<size_t, int> &occupied)
    {
        occupied_data_pkgs_.push_back(occupied);
    }

    Arrayi CellIndexFromPositionOnGlobalMesh(const Vecd &position)
    {
        return global_mesh_.CellIndexFromPosition(position);
    }

    Vecd GridPositionFromIndexOnGlobalMesh(const Arrayi &cell_index)
    {
        return global_mesh_.GridPositionFromIndex(cell_index);
    }

    bool isCoreDataPackage(const Arrayi &cell_index)
    {
        size_t package_index = PackageIndexFromCellIndex(cell_index);
        return meta_data_cell_[package_index].second == 1;
    }

    void organizeOccupiedPackages(){
        parallel_sort(occupied_data_pkgs_.begin(), occupied_data_pkgs_.end(),
                    [](const std::pair<size_t, int>& a, const std::pair<size_t, int>& b)
                    {
                        return a.first < b.first; 
                    });
        num_grid_pkgs_ = occupied_data_pkgs_.size() + 2;
        cell_neighborhood_ = new CellNeighborhood[num_grid_pkgs_];
        meta_data_cell_ = new std::pair<Arrayi, int>[num_grid_pkgs_];
    }

    bool isWithinCorePackage(Vecd position)
    {
        Arrayi cell_index = CellIndexFromPosition(position);
        return isCoreDataPackage(cell_index);
    }

    template <typename DataType>
    MeshVariable<DataType> *registerMeshVariable(const std::string &variable_name)
    {
        MeshVariable<DataType> *variable =
            findVariableByName<DataType>(all_mesh_variables_, variable_name);
        if (variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            size_t new_variable_index = std::get<type_index>(all_mesh_variables_).size();
            return addVariableToAssemble<DataType>(all_mesh_variables_, mesh_variable_ptrs_,
                                                   variable_name, new_variable_index);
        }
        return variable;
    };
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
