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
#include "sphinxsys_variable.h"
#include "tbb/parallel_sort.h"
#include "execution_policy.h"

#include <algorithm>
#include <fstream>
#include <functional>
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
  public:
    template <typename... Args>
    explicit MeshWithGridDataPackages(BoundingBox tentative_bounds,
                                      Real data_spacing,
                                      size_t buffer_size)
        : Mesh(tentative_bounds, pkg_size * data_spacing, buffer_size),
          global_mesh_(mesh_lower_bound_ + 0.5 * data_spacing * Vecd::Ones(), data_spacing, all_cells_ * pkg_size),
          data_spacing_(data_spacing),
          index_handler_("index_handler", IndexHandler{data_spacing_, all_cells_, *static_cast<Mesh*>(this)})
    {
        allocateIndexDataMatrix(all_cells_);
    };
    virtual ~MeshWithGridDataPackages(){};

    /** spacing between the data, which is 1/ pkg_size of this grid spacing */
    Real DataSpacing() { return data_spacing_; };
    Real GridSpacing() { return grid_spacing_; };
    size_t BufferWidth() { return buffer_width_; };
    Mesh global_mesh_;                                         /**< the mesh for the locations of all possible data points. */
    size_t num_grid_pkgs_ = 2;                                 /**< the number of all distinct packages, initially only 2 singular packages. */
    DiscreteVariable<std::pair<Arrayi, int>> meta_data_cell_{"meta_data_cell", 2};          /**< metadata for each occupied cell: (arrayi)cell index, (int)core1/inner0. */
    DiscreteVariable<CellNeighborhood> cell_neighborhood_{"mesh_cell_neighborhood", 2};     /**< 3*3(*3) array to store indicies of neighborhood cells. */
    DiscreteVariable<size_t> cell_package_index_{"cell_package_index", 2};                  /**< the package index for each cell in a 1-d array. */
    ConcurrentVec<std::pair<size_t, int>> occupied_data_pkgs_; /**< (size_t)sort_index, (int)core1/inner0. */

  private:
    DataContainerUniquePtrAssemble<MeshVariable> mesh_variable_ptrs_;
    MeshVariableAssemble all_mesh_variables_;         /**< all mesh variables on this mesh. */
    static constexpr int pkg_size = PKG_SIZE;         /**< the size of the data package matrix. */
    const Real data_spacing_;                         /**< spacing of data in the data packages. */

    /**< allocate memories for metadata of data packages. */
    void allocateIndexDataMatrix(Array2i mesh_size){ cell_package_index_.reallocateData(par, all_cells_[0] * all_cells_[1]); };
    void allocateIndexDataMatrix(Array3i mesh_size){ cell_package_index_.reallocateData(par, all_cells_[0] * all_cells_[1] * all_cells_[2]); };

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
                variable->reallocateData(par, num_grid_pkgs_);
            }
        }
    };
    OperationOnDataAssemble<MeshVariableAssemble, ResizeMeshVariableData>
        resize_mesh_variable_data_{};

    struct SyncMeshVariableData
    {
        template <typename DataType, class ExecutionPolicy>
        void operator()(DataContainerAddressKeeper<MeshVariable<DataType>> &all_mesh_variables_,
                        ExecutionPolicy &ex_policy)
        {
            for(size_t l = 0; l != all_mesh_variables_.size(); l++)
            {
                MeshVariable<DataType> *variable = all_mesh_variables_[l];
                variable->prepareForOutput(ex_policy);
            }
        }
    };
    OperationOnDataAssemble<MeshVariableAssemble, SyncMeshVariableData>
        sync_mesh_variable_data_{};
  public:
    /** wrapper for all index exchange related functions. */
    struct IndexHandler
    {
        Real data_spacing_;
        Arrayi all_cells_;
        Mesh mesh_;

        Arrayi CellIndexFromPosition(const Vecd &position) const
        {
            return mesh_.CellIndexFromPosition(position);
        }

        /** return the position of the lower bound data in a cell. */
        Vecd DataLowerBoundInCell(const Arrayi &cell_index)
        {
            return mesh_.CellLowerCornerPosition(cell_index) + 0.5 * data_spacing_ * Vecd::Ones();
        }

        /** return the grid index from its position and the index of the cell it belongs to. */
        Arrayi DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position)
        {
            return floor((position - DataLowerBoundInCell(cell_index)).array() / data_spacing_)
                .template cast<int>()
                .max(Arrayi::Zero())
                .min((pkg_size - 1) * Arrayi::Ones());
        }

        /** return the position of data from its local grid index and the index of the cell it belongs to. */
        Vecd DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index)
        {
            return DataLowerBoundInCell(cell_index) + data_index.cast<Real>().matrix() * data_spacing_;
        }

        /** return the package index in the data array from the cell index it belongs to. */
        size_t PackageIndexFromCellIndex(size_t *cell_package_index, const Arrayi &cell_index)
        {
            size_t index_1d = mesh_.transferMeshIndexTo1D(all_cells_, cell_index);
            return cell_package_index[index_1d];
        }
    };

    SingularVariable<IndexHandler> index_handler_;
  public:
    void resizeMeshVariableData() { resize_mesh_variable_data_(all_mesh_variables_, num_grid_pkgs_); };
    template <class ExecutionPolicy>
    void syncMeshVariableData(ExecutionPolicy &ex_policy) { sync_mesh_variable_data_(all_mesh_variables_, ex_policy); };

    template <typename DataType>
    MeshVariable<DataType> *registerMeshVariable(const std::string &variable_name)
    {
        MeshVariable<DataType> *variable =
            findVariableByName<DataType, MeshVariable>(all_mesh_variables_, variable_name);
        if (variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            size_t new_variable_index = std::get<type_index>(all_mesh_variables_).size();
            return addVariableToAssemble<DataType, MeshVariable>(all_mesh_variables_, mesh_variable_ptrs_,
                                                   variable_name, new_variable_index);
        }
        return variable;
    }

    /** return the mesh variable according to the name registered */
    template <typename DataType>
    MeshVariable<DataType> *getMeshVariable(const std::string &variable_name)
    {
        return findVariableByName<DataType, MeshVariable>(all_mesh_variables_, variable_name);
    }

    void registerOccupied(size_t sort_index, int type)
    {
        occupied_data_pkgs_.push_back(std::make_pair(sort_index, type));
    }

    void organizeOccupiedPackages()
    {
        parallel_sort(occupied_data_pkgs_.begin(), occupied_data_pkgs_.end(),
                      [](const std::pair<size_t, int>& a, const std::pair<size_t, int>& b)
                      {
                          return a.first < b.first; 
                      });
        num_grid_pkgs_ = occupied_data_pkgs_.size() + 2;
        cell_neighborhood_.reallocateData(par, num_grid_pkgs_);
        meta_data_cell_.reallocateData(par, num_grid_pkgs_);
    }

    bool isInnerDataPackage(const Arrayi &cell_index)
    {
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        /**
         * NOTE currently this func is only used in non-device mode;
         *      use the `DelegatedDataField` version when needed.
         */
        return cell_package_index_.Data()[index_1d] > 1;
    }
    bool isWithinCorePackage(size_t *cell_package_index,
                             std::pair<Arrayi, int> *meta_data_cell,
                             Vecd position)
    {
        Arrayi cell_index = CellIndexFromPosition(position);
        size_t package_index = PackageIndexFromCellIndex(cell_package_index, cell_index);
        return meta_data_cell[package_index].second == 1;
    }

    //TODO eliminate the function with the index handler
    /** return the package index in the data array from the cell index it belongs to. */
    size_t PackageIndexFromCellIndex(size_t *cell_package_index, const Arrayi &cell_index)
    {
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        return cell_package_index[index_1d];
    }
    void assignDataPackageIndex(const Arrayi &cell_index, const size_t package_index){
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        /**
         * NOTE currently the `cell_package_index_` is only assigned in the host;
         *      use the `DelegatedDataField` version when needed.
         */
        cell_package_index_.Data()[index_1d] = package_index;
    }
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
