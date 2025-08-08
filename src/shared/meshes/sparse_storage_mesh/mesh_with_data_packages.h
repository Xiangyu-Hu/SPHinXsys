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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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

#include "execution_policy.h"
#include "grid_data_package_type.hpp"
#include "sphinxsys_variable.h"

#include "tbb/parallel_sort.h"

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
 * Each data package is again a mesh, but grid based with pkg_size grids on each dimension.
 * The operation on field data is achieved by mesh dynamics.
 * Note that a data package should be not near the mesh bound, otherwise one will encounter the error "out of range".
 */
template <size_t PKG_SIZE>
class MeshWithGridDataPackages : public Mesh
{
  public:
    template <class DataType>
    using MeshVariableData = PackageDataMatrix<DataType, PKG_SIZE>;
    template <typename DataType>
    using MeshVariable = DiscreteVariable<MeshVariableData<DataType>>;
    template <typename DataType>
    using BKGMeshVariable = DiscreteVariable<DataType>;
    static constexpr int pkg_size = PKG_SIZE; /**< the size of the data package matrix. */

  protected:
    typedef DataContainerAddressAssemble<MeshVariable> MeshVariableAssemble;
    DataContainerUniquePtrAssemble<MeshVariable> mesh_variable_ptrs_;
    MeshVariableAssemble all_mesh_variables_;     /**< all mesh variables on this mesh. */
    MeshVariableAssemble mesh_variable_to_write_; /**< mesh variables to write, which are not empty. */
    typedef DataContainerAddressAssemble<DiscreteVariable> BKGMeshVariableAssemble;
    DataContainerUniquePtrAssemble<DiscreteVariable> bkg_mesh_variable_ptrs_;
    BKGMeshVariableAssemble all_bkg_mesh_variables_;     /**< all discrete variables on this mesh. */
    BKGMeshVariableAssemble bkg_mesh_variable_to_write_; /**< discrete variables to write, which are not empty. */

  public:
    template <typename... Args>
    explicit MeshWithGridDataPackages(
        BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, UnsignedInt num_singular_pkgs = 2)
        : Mesh(tentative_bounds, pkg_size * data_spacing, buffer_size),
          global_mesh_(mesh_lower_bound_ + 0.5 * data_spacing * Vecd::Ones(), data_spacing, all_cells_ * pkg_size),
          num_singular_pkgs_(num_singular_pkgs), num_grid_pkgs_(num_singular_pkgs),
          dv_pkg_cell_info_("PackageCellInfo", num_singular_pkgs_),
          cell_neighborhood_("CellNeighborhood", num_singular_pkgs_),
          cell_pkg_index_(*registerBKGMeshVariable<UnsignedInt>("CellPackageIndex")),
          data_spacing_(data_spacing),
          index_handler_("index_handler", IndexHandler{data_spacing_, all_cells_, *static_cast<Mesh *>(this)}){};
    virtual ~MeshWithGridDataPackages() {};

    /** spacing between the data, which is 1/ pkg_size of this grid spacing */
    Real DataSpacing() { return data_spacing_; };
    Real GridSpacing() { return grid_spacing_; };
    size_t BufferWidth() { return buffer_width_; };
    int DataPackageSize() { return pkg_size; };
    Mesh global_mesh_;                                         /**< singular packages used for far field. */
    UnsignedInt num_singular_pkgs_;                            /**< the number of all packages, initially only singular packages. */
    UnsignedInt num_grid_pkgs_;                                /**< the number of all packages, initially only with singular packages. */
    DiscreteVariable<std::pair<Arrayi, int>> dv_pkg_cell_info_;   /**< metadata for each occupied cell: (arrayi)cell index, (int)core1/inner0. */
    DiscreteVariable<CellNeighborhood> cell_neighborhood_;     /**< 3*3(*3) array to store indicies of neighborhood cells. */
    BKGMeshVariable<UnsignedInt> &cell_pkg_index_;             /**< the package index for each cell in a 1-d array. */
    ConcurrentVec<std::pair<size_t, int>> occupied_data_pkgs_; /**< (size_t)sort_index, (int)core1/inner0. */
    UnsignedInt NumSingularPackages() const { return num_singular_pkgs_; };

    UnsignedInt NumGridPackages()
    {
        return checkOrganized("NumGridPackages", num_grid_pkgs_);
    };

    DiscreteVariable<std::pair<Arrayi, int>> &dvPkgCellInfo()
    {
        return checkOrganized("dvPkgCellInfo", dv_pkg_cell_info_);
    };

    DiscreteVariable<CellNeighborhood> &getCellNeighborhood()
    {
        return checkOrganized("getCellNeighborhood", cell_neighborhood_);
    };

    DiscreteVariable<size_t> &getCellPackageIndex()
    {
        return checkOrganized("getCellPackageIndex", cell_pkg_index_);
    };

    template <typename DataType>
    void addMeshVariableToWrite(const std::string &variable_name)
    {
        addVariableToList<MeshVariable, DataType>(mesh_variable_to_write_, all_mesh_variables_, variable_name);
    };

    template <typename DataType>
    void addBKGMeshVariableToWrite(const std::string &variable_name)
    {
        addVariableToList<DiscreteVariable, DataType>(bkg_mesh_variable_to_write_, all_bkg_mesh_variables_, variable_name);
    };

    void writeMeshVariableToPlt(std::ofstream &output_file);
    void writeBKGMeshVariableToPlt(std::ofstream &output_file);

  protected:
    /** Generalized mesh data type */
    const Real data_spacing_;   /**< spacing of data in the data packages. */
    bool is_organized_ = false; /**< whether the data packages are organized. */

    template <typename T>
    T &checkOrganized(std::string func_name, T &value)
    {
        if (!is_organized_)
        {
            std::cout << "\n Error: the mesh is not organized! (called from " << func_name << ")" << std::endl;
            exit(1);
        }
        return value;
    };

    template <template <typename> typename ContainerType, typename DataType>
    void addVariableToList(DataContainerAddressAssemble<ContainerType> &variable_set,
                           DataContainerAddressAssemble<ContainerType> &all_variable_set,
                           const std::string &variable_name)
    {
        ContainerType<DataType> *variable =
            findVariableByName<DataType, ContainerType>(all_variable_set, variable_name);

        if (variable == nullptr)
        {
            std::cout << "\n Error: the mesh variable '" << variable_name << "' is  not exist!" << std::endl;
            exit(1);
        }

        ContainerType<DataType> *listed_variable =
            findVariableByName<DataType, ContainerType>(variable_set, variable_name);
        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
        }
    };

    struct SyncMeshVariableData
    {
        template <typename DataType, class ExecutionPolicy>
        void operator()(DataContainerAddressKeeper<MeshVariable<DataType>> &all_mesh_variables_,
                        ExecutionPolicy &ex_policy)
        {
            for (size_t l = 0; l != all_mesh_variables_.size(); l++)
            {
                MeshVariable<DataType> *variable = all_mesh_variables_[l];
                variable->prepareForOutput(ex_policy);
            }
        }
    };
    OperationOnDataAssemble<MeshVariableAssemble, SyncMeshVariableData> sync_mesh_variable_data_{};

    void fillFarFieldCellNeighborhood(CellNeighborhood *neighbor)
    {
        for (size_t i = 0; i != num_singular_pkgs_; i++)
        {
            mesh_for_each(
                -Arrayi::Ones(), Arrayi::Ones() * 2,
                [&](const Arrayi &index)
                {
                    neighbor[i](index + Arrayi::Ones()) = i;
                });
        }
    };

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
    template <class ExecutionPolicy>
    void syncMeshVariableData(ExecutionPolicy &ex_policy) { sync_mesh_variable_data_(all_mesh_variables_, ex_policy); };

    template <template <typename> typename ContainerType, typename DataType, typename... Args>
    ContainerType<DataType> *registerVariable(DataContainerAddressAssemble<ContainerType> &all_variable_set,
                                              DataContainerUniquePtrAssemble<ContainerType> &all_variable_ptrs_,
                                              const std::string &variable_name, Args &&...args)
    {
        ContainerType<DataType> *variable =
            findVariableByName<DataType, ContainerType>(all_variable_set, variable_name);
        if (variable == nullptr)
        {
            return addVariableToAssemble<DataType, ContainerType>(
                all_variable_set, all_variable_ptrs_, variable_name, std::forward<Args>(args)...);
        }
        return variable;
    }

    template <typename DataType>
    MeshVariable<DataType> *registerMeshVariable(const std::string &variable_name)
    {
        if (!is_organized_)
        {
            std::cout << "\n Error: the mesh variable '" << variable_name
                      << "' is registered before the data packages are organized!" << std::endl;
            exit(1);
        }
        return registerVariable<MeshVariable, DataType>(
            all_mesh_variables_, mesh_variable_ptrs_, variable_name, num_grid_pkgs_);
    }

    template <typename DataType, typename... Args>
    BKGMeshVariable<DataType> *registerBKGMeshVariable(const std::string &variable_name, Args &&...args)
    {
        return registerVariable<BKGMeshVariable, DataType>(
            all_bkg_mesh_variables_, bkg_mesh_variable_ptrs_, variable_name, AllCells().prod(),
            std::forward<Args>(args)...);
    }

    /** return the mesh variable according to the name registered */
    template <typename DataType>
    MeshVariable<DataType> *getMeshVariable(const std::string &variable_name)
    {
        return findVariableByName<DataType, MeshVariable>(all_mesh_variables_, variable_name);
    }

    template <typename DataType>
    BKGMeshVariable<DataType> *getBKGMeshVariable(const std::string &variable_name)
    {
        return findVariableByName<DataType, BKGMeshVariable>(all_bkg_mesh_variables_, variable_name);
    }

    void registerOccupied(size_t sort_index, int type)
    {
        occupied_data_pkgs_.push_back(std::make_pair(sort_index, type));
    }

    void organizeOccupiedPackages()
    {
        parallel_sort(
            occupied_data_pkgs_.begin(), occupied_data_pkgs_.end(),
            [](const std::pair<size_t, int> &a, const std::pair<size_t, int> &b)
            {
                return a.first < b.first;
            });
        num_grid_pkgs_ = occupied_data_pkgs_.size() + num_singular_pkgs_;
        cell_neighborhood_.reallocateData(par, num_grid_pkgs_);
        fillFarFieldCellNeighborhood(cell_neighborhood_.Data());
        dv_pkg_cell_info_.reallocateData(par, num_grid_pkgs_);
        is_organized_ = true;
    }

    bool isInnerDataPackage(const Arrayi &cell_index)
    {
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        /**
         * NOTE currently this func is only used in non-device mode;
         *      use the `DelegatedData` version when needed.
         */
        return cell_pkg_index_.Data()[index_1d] > 1;
    }

    bool isWithinCorePackage(size_t *cell_package_index,
                             std::pair<Arrayi, int> *meta_data_cell,
                             Vecd position)
    {
        Arrayi cell_index = CellIndexFromPosition(position);
        size_t package_index = PackageIndexFromCellIndex(cell_package_index, cell_index);
        return meta_data_cell[package_index].second == 1;
    }

    /** return the package index in the data array from the cell index it belongs to. */
    size_t PackageIndexFromCellIndex(size_t *cell_package_index, const Arrayi &cell_index)
    {
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        return cell_package_index[index_1d];
    }
    void assignDataPackageIndex(const Arrayi &cell_index, const size_t package_index)
    {
        size_t index_1d = transferMeshIndexTo1D(all_cells_, cell_index);
        /**
         * NOTE currently the `cell_pkg_index_` is only assigned in the host;
         *      use the `DelegatedData` version when needed.
         */
        cell_pkg_index_.Data()[index_1d] = package_index;
    }
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
