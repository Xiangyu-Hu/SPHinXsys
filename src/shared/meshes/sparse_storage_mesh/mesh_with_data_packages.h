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
template <UnsignedInt PKG_SIZE>
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
    MeshVariableAssemble all_mesh_variables_;      /**< all mesh variables on this mesh. */
    MeshVariableAssemble mesh_variables_to_write_;  /**< mesh variables to write, which are not empty. */
    MeshVariableAssemble mesh_variables_to_probe_; /**< mesh variables to probe. */
    typedef DataContainerAddressAssemble<DiscreteVariable> BKGMeshVariableAssemble;
    DataContainerUniquePtrAssemble<DiscreteVariable> bkg_mesh_variable_ptrs_;
    BKGMeshVariableAssemble all_bkg_mesh_variables_;     /**< all discrete variables on this mesh. */
    BKGMeshVariableAssemble bkg_mesh_variables_to_write_; /**< discrete variables to write, which are not empty. */

  public:
    MeshWithGridDataPackages(BoundingBox tentative_bounds, Real data_spacing,
                             UnsignedInt buffer_size, UnsignedInt num_singular_pkgs = 2);
    virtual ~MeshWithGridDataPackages() {};

    /** spacing between the data, which is 1/ pkg_size of this grid spacing */
    Real DataSpacing() { return data_spacing_; };
    Real GridSpacing() { return grid_spacing_; };
    UnsignedInt BufferWidth() { return buffer_width_; };
    int DataPackageSize() { return pkg_size; };
    UnsignedInt NumSingularPackages() const { return num_singular_pkgs_; };
    UnsignedInt NumGridPackages();
    DiscreteVariable<std::pair<Arrayi, int>> &dvPkgCellInfo();
    DiscreteVariable<CellNeighborhood> &getCellNeighborhood();
    DiscreteVariable<UnsignedInt> &getCellPackageIndex();
    ConcurrentVec<std::pair<UnsignedInt, int>> &getOccupiedDataPackages();
    template <typename DataType>
    void addMeshVariableToWrite(const std::string &variable_name);
    template <typename DataType>
    void addMeshVariableToProbe(const std::string &variable_name);
    template <typename DataType>
    void addBKGMeshVariableToWrite(const std::string &variable_name);
    void writeMeshVariableToPlt(std::ofstream &output_file);
    void writeBKGMeshVariableToPlt(std::ofstream &output_file);

  protected:
    Mesh global_mesh_;                                              /**< the global mesh with the size of data spacing. */
    UnsignedInt num_singular_pkgs_;                                 /**< the number of all packages, initially only singular packages. */
    UnsignedInt num_grid_pkgs_;                                     /**< the number of all packages, initially only with singular packages. */
    DiscreteVariable<std::pair<Arrayi, int>> dv_pkg_cell_info_;     /**< metadata for each occupied cell: (arrayi)cell index, (int)core1/inner0. */
    DiscreteVariable<CellNeighborhood> cell_neighborhood_;          /**< 3*3(*3) array to store indicies of neighborhood cells. */
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;              /**< the package index for each cell in a 1-d array. */
    ConcurrentVec<std::pair<UnsignedInt, int>> occupied_data_pkgs_; /**< (UnsignedInt)sort_index, (int)core1/inner0. */
    const Real data_spacing_;                                       /**< spacing of data in the data packages. */
    bool is_organized_ = false;                                     /**< whether the data packages are organized. */

    template <typename T>
    T &checkOrganized(std::string func_name, T &value);

    template <template <typename> typename ContainerType, typename DataType>
    void addVariableToList(DataContainerAddressAssemble<ContainerType> &variable_set,
                           DataContainerAddressAssemble<ContainerType> &all_variable_set,
                           const std::string &variable_name);

    template <template <typename> class MeshVariableType>
    struct SyncMeshVariableData
    {
        template <typename DataType, class ExecutionPolicy>
        void operator()(DataContainerAddressKeeper<MeshVariableType<DataType>> &mesh_variables,
                        ExecutionPolicy &ex_policy);
    };
    OperationOnDataAssemble<MeshVariableAssemble, SyncMeshVariableData<MeshVariable>> sync_mesh_variable_data_{};
    OperationOnDataAssemble<BKGMeshVariableAssemble, SyncMeshVariableData<BKGMeshVariable>> sync_bkg_mesh_variable_data_{};

  public:
    /** wrapper for all index exchange related functions. */
    struct IndexHandler
    {
        Real data_spacing_;
        Arrayi all_cells_;
        Mesh mesh_;

        Arrayi CellIndexFromPosition(const Vecd &position) const;
        Vecd DataLowerBoundInCell(const Arrayi &cell_index);
        Arrayi DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position);
        Vecd DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index);
        UnsignedInt PackageIndexFromCellIndex(UnsignedInt *cell_package_index, const Arrayi &cell_index);
    };
    SingularVariable<IndexHandler> index_handler_;

  public:
    template <class ExecutionPolicy>
    void syncMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncBKGMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncMeshVariablesToProbe(ExecutionPolicy &ex_policy);
    template <template <typename> typename ContainerType, typename DataType, typename... Args>
    ContainerType<DataType> *registerVariable(DataContainerAddressAssemble<ContainerType> &all_variable_set,
                                              DataContainerUniquePtrAssemble<ContainerType> &all_variable_ptrs_,
                                              const std::string &variable_name, Args &&...args);
    template <typename DataType>
    MeshVariable<DataType> *registerMeshVariable(const std::string &variable_name);
    template <typename DataType, typename... Args>
    BKGMeshVariable<DataType> *registerBKGMeshVariable(const std::string &variable_name, Args &&...args);
    template <typename DataType>
    MeshVariable<DataType> *getMeshVariable(const std::string &variable_name);
    template <typename DataType>
    BKGMeshVariable<DataType> *getBKGMeshVariable(const std::string &variable_name);

    void registerOccupied(UnsignedInt sort_index, int type);
    void organizeOccupiedPackages();
    bool isInnerDataPackage(const Arrayi &cell_index);
    bool isWithinCorePackage(UnsignedInt *cell_package_index, std::pair<Arrayi, int> *meta_data_cell, Vecd position);
    Arrayi boundCellIndex(const Arrayi &input) const;
    UnsignedInt PackageIndexFromCellIndex(UnsignedInt *cell_package_index, const Arrayi &cell_index);
    void assignDataPackageIndex(const Arrayi &cell_index, const UnsignedInt package_index);
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
