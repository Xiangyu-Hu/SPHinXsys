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

#include "base_mesh.hpp"

#include "execution_policy.h"
#include "data_package_type.hpp"

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
 * Each data package is again a mesh, but grid based with PKG_SIZE grids on each dimension.
 * The operation on field data is achieved by mesh dynamics.
 * Note that a data package should be not near the mesh bound, otherwise one will encounter the error "out of range".
 */
template <int PKG_SIZE>
class MeshWithGridDataPackages
{
  public:
    template <class DataType>
    using MeshVariableData = PackageData<DataType, PKG_SIZE>;
    template <typename DataType>
    using MeshVariable = DiscreteVariable<MeshVariableData<DataType>>;
    template <typename DataType>
    using BKGMeshVariable = DiscreteVariable<DataType>;
    template <typename DataType>
    using MetaVariable = DiscreteVariable<DataType>;

    /** wrapper for all index exchange related functions. */
    class IndexHandler : public Mesh
    {
        Real data_spacing_;

      public:
        IndexHandler(BoundingBoxd tentative_bounds, Real grid_spacing,
                     UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset,
                     Real data_spacing);
        Vecd DataLowerBoundInCell(const Arrayi &cell_index) const;
        Arrayi DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position) const;
        Vecd DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index) const;
        UnsignedInt PackageIndexFromCellIndex(UnsignedInt *cell_package_index, const Arrayi &cell_index) const;
        bool isWithinCorePackage(UnsignedInt *cell_package_index, int *pkg_type, const Vecd &position);
        Real DataSpacing() const { return data_spacing_; };
    };
    typedef DataContainerAddressAssemble<MeshVariable> MeshVariableAssemble;
    typedef DataContainerAddressAssemble<BKGMeshVariable> BKGMeshVariableAssemble;
    typedef DataContainerAddressAssemble<MetaVariable> MetaVariableAssemble;

  protected:
    DataContainerUniquePtrAssemble<MeshVariable> mesh_variable_ptrs_;
    MeshVariableAssemble all_mesh_variables_;      /**< all mesh variables on this mesh. */
    MeshVariableAssemble mesh_variables_to_write_; /**< mesh variables to write, which are not empty. */
    MeshVariableAssemble mesh_variables_to_probe_; /**< mesh variables to probe. */
    MeshVariableAssemble evolving_mesh_variables_;

    DataContainerUniquePtrAssemble<BKGMeshVariable> bkg_mesh_variable_ptrs_;
    BKGMeshVariableAssemble all_bkg_mesh_variables_;      /**< all discrete variables on this mesh. */
    BKGMeshVariableAssemble bkg_mesh_variables_to_write_; /**< discrete variables to write, which are not empty. */

    DataContainerUniquePtrAssemble<MetaVariable> meta_variable_ptrs_;
    MetaVariableAssemble all_meta_variables_;
    MetaVariableAssemble evolving_meta_variables_;

    UniquePtrsKeeper<Entity> unique_variable_ptrs_;

  public:
    MeshWithGridDataPackages(BoundingBoxd tentative_bounds, Real data_spacing,
                             UnsignedInt buffer_size, UnsignedInt num_singular_pkgs = 2);
    virtual ~MeshWithGridDataPackages() {};

    static constexpr int DataPackageSize() { return PKG_SIZE; };
    UnsignedInt NumSingularPackages() const { return num_singular_pkgs_; };
    UnsignedInt PackageBound() const { return pkgs_bound_; };
    SingularVariable<UnsignedInt> &svNumGridPackages();
    BKGMeshVariable<UnsignedInt> &getCellPackageIndex() { return *bmv_cell_pkg_index_; };
    ConcurrentVec<std::pair<UnsignedInt, int>> &getOccupiedDataPackages() { return occupied_data_pkgs_; };
    MetaVariable<CellNeighborhood> &getCellNeighborhood();
    MetaVariable<UnsignedInt> &getPackage1DCellIndex();
    MetaVariable<int> &getPackageType();
    MetaVariableAssemble &getEvolvingMetaVariables() { return evolving_meta_variables_; };
    MeshVariableAssemble &getEvolvingMeshVariables() { return evolving_mesh_variables_; };

  protected:
    IndexHandler index_handler_;
    UnsignedInt num_singular_pkgs_;                  /**< the number of all packages, initially only singular packages. */
    SingularVariable<UnsignedInt> sv_num_grid_pkgs_; /**< the number of all packages, initially only with singular packages. */
    UnsignedInt pkgs_bound_;
    MetaVariable<UnsignedInt> *dv_pkg_1d_cell_index_;               /**< metadata for data pckages: cell index. */
    MetaVariable<int> *dv_pkg_type_;                                /**< metadata for data pckages: (int)core1/inner0. */
    MetaVariable<CellNeighborhood> *cell_neighborhood_;             /**< 3*3(*3) array to store indicies of neighborhood cells. */
    BKGMeshVariable<UnsignedInt> *bmv_cell_pkg_index_;              /**< the package index for each cell in a 1-d array. */
    ConcurrentVec<std::pair<UnsignedInt, int>> occupied_data_pkgs_; /**< (UnsignedInt)sort_index, (int)core1/inner0. */
    bool is_organized_ = false;                                     /**< whether the data packages are organized. */
    Mesh global_mesh_;                                              /**< the global mesh with the size of data spacing. */

    template <typename T>
    T &checkOrganized(std::string func_name, T &value);

    OperationOnDataAssemble<MeshVariableAssemble, PrepareVariablesToWrite<MeshVariable>> sync_mesh_variable_data_{};
    OperationOnDataAssemble<BKGMeshVariableAssemble, PrepareVariablesToWrite<BKGMeshVariable>> sync_bkg_mesh_variable_data_{};

    template <typename DataType>
    DataType DataValueFromGlobalIndex(PackageData<DataType, PKG_SIZE> *pkg_data,
                                      const Arrayi &global_grid_index, UnsignedInt *cell_package_index);

  public:
    IndexHandler &getIndexHandler() { return index_handler_; };
    template <class ExecutionPolicy>
    void syncMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncBKGMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncMeshVariablesToProbe(ExecutionPolicy &ex_policy);
    template <typename DataType>
    MeshVariable<DataType> *registerMeshVariable(const std::string &variable_name);
    template <typename DataType, typename... Args>
    BKGMeshVariable<DataType> *registerBKGMeshVariable(const std::string &variable_name, Args &&...args);
    template <typename DataType, typename... Args>
    MetaVariable<DataType> *registerMetaVariable(const std::string &variable_name, Args &&...args);
    template <typename DataType>
    MeshVariable<DataType> *getMeshVariable(const std::string &variable_name);
    template <typename DataType>
    BKGMeshVariable<DataType> *getBKGMeshVariable(const std::string &variable_name);
    template <typename DataType>
    MetaVariable<DataType> *getMetaVariable(const std::string &variable_name);
    template <typename DataType>
    void addMeshVariableToWrite(const std::string &variable_name);
    template <typename DataType>
    void addMeshVariableToProbe(const std::string &variable_name);
    template <typename DataType>
    void addBKGMeshVariableToWrite(const std::string &variable_name);
    template <typename DataType>
    void addEvolvingMetaVariable(const std::string &variable_name);
    void writeMeshVariableToPlt(std::ofstream &output_file);
    void writeBKGMeshVariableToPlt(std::ofstream &output_file);

    void organizeOccupiedPackages();
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
