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
 * @file    spares_mesh_field.h
 * @brief   This class is designed to save memory and increase computational efficiency on mesh.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_H
#define MESH_WITH_DATA_PACKAGES_H

#include "base_mesh.hpp"

#include "data_package_type.hpp"
#include "execution_policy.h"

namespace SPH
{
template <int PKG_SIZE>
class PackageMesh : public Mesh
{
    Real data_spacing_;

  public:
    PackageMesh() {};
    PackageMesh(BoundingBoxd tentative_bounds, Real grid_spacing,
                UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset);
    Vecd DataLowerBoundInCell(const Arrayi &cell_index) const;
    Arrayi DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position) const;
    Vecd DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index) const;
    UnsignedInt PackageIndexFromCellIndex(UnsignedInt *cell_pkg_index, const Arrayi &cell_index) const;
    bool isWithinPackageType(int target_type, UnsignedInt *cell_pkg_index, int *pkg_type, const Vecd &position);
    Real DataSpacing() const { return data_spacing_; };
    Mesh GlobalMesh() const;
    template <typename DataType>
    DataType ValueByGlobalMesh(PackageData<DataType, PKG_SIZE> *pkg_data, const Arrayi &global_grid_index,
                               UnsignedInt *cell_package_index) const;
};

/**
 * @class SparseMeshField
 * @brief Class for mesh field using grid-based data packages.
 * @details The idea is to save sparse data on a cell-based mesh.
 * We say sparse data, it means that only in some inner mesh cells there are no trivial data.
 * A typical example is a level-set field which only has meaningful values near the interface,
 * while the latter is in the inner region of a mesh.
 * In this class, only some inner mesh cells are filled with data packages.
 * Each data package is again grid based with PKG_SIZE grids on each dimension.
 * The operation on field data is achieved by mesh dynamics.
 * Note that a data package should be not near the mesh bound,
 * otherwise one will encounter the error "out of range".
 */
template <int PKG_SIZE>
class SparseMeshField : public MultiResolutionMeshField<PackageMesh<PKG_SIZE>>
{
  public:
    template <class DataType>
    using PackageVariableData = PackageData<DataType, PKG_SIZE>;
    template <typename DataType>
    using PackageVariable = DiscreteVariable<PackageVariableData<DataType>>;
    template <typename DataType>
    using MetaVariable = DiscreteVariable<DataType>;
    using IndexHandler = PackageMesh<PKG_SIZE>;
    typedef DataContainerAddressAssemble<PackageVariable> PackageVariableAssemble;
    typedef DataContainerAddressAssemble<MetaVariable> MetaVariableAssemble;

  public:
    SparseMeshField(const std::string &name, UnsignedInt resolution_levels,
                    BoundingBoxd tentative_bounds, Real Reference_grid_spacing,
                    UnsignedInt buffer_width, UnsignedInt num_boundary_pkgs);
    virtual ~SparseMeshField() {};
    static constexpr int PackageDataSize() { return PKG_SIZE; };
    UnsignedInt NumBoundaryPackages() const { return num_boundary_pkgs_; };
    UnsignedInt PackageBound() const { return pkgs_bound_; };
    StdVec<UnsignedInt> &getNumPackageOffsets() { return num_pkgs_offsets_; };
    CellVariable<UnsignedInt> &getCellPackageIndex() { return *mcv_cell_pkg_index_; };
    ConcurrentVec<std::pair<UnsignedInt, int>> &getOccupiedPackageDatas() { return occupied_data_pkgs_; };
    MetaVariable<CellNeighborhood> &getCellNeighborhood();
    MetaVariable<UnsignedInt> &getPackage1DCellIndex();
    MetaVariable<int> &getPackageType();
    MetaVariableAssemble &getEvolvingMetaVariables() { return evolving_meta_variables_; };
    PackageVariableAssemble &getEvolvingPackageVariables() { return evolving_pkg_variables_; };
    void writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence = 0) override;
    template <class ExecutionPolicy>
    void syncPackageVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncPackageVariablesToProbe(ExecutionPolicy &ex_policy);
    template <typename DataType>
    PackageVariable<DataType> *registerPackageVariable(const std::string &name);
    template <typename DataType, typename... Args>
    MetaVariable<DataType> *registerMetaVariable(const std::string &name, Args &&...args);
    template <typename DataType>
    PackageVariable<DataType> *getPackageVariable(const std::string &name);
    template <typename DataType>
    MetaVariable<DataType> *getMetaVariable(const std::string &name);
    template <typename DataType>
    void addPackageVariableToWrite(const std::string &name);
    template <typename DataType>
    void addPackageVariableToProbe(const std::string &name);
    template <typename DataType>
    void addEvolvingMetaVariable(const std::string &name);
    template <class DiscreteVariableType, class BoundaryDataFunction>
    void setBoundaryData(DiscreteVariableType *variable, UnsignedInt resolution_level,
                         const BoundaryDataFunction &boundary_data_function);
    void organizeOccupiedPackages(UnsignedInt resolution_level);

    template <typename DataType>
    class ProbeMesh
    {
      public:
        template <class ExecutionPolicy>
        ProbeMesh(const ExecutionPolicy &ex_policy, SparseMeshField<PKG_SIZE> &encloser,
                  const std::string &variable_name);
        DataType operator()(IndexHandler &index_handler, const Vecd &position);

      protected:
        PackageData<DataType, PKG_SIZE> *pkg_data_;
        IndexHandler *index_handler_;
        UnsignedInt resolution_levels_;
        UnsignedInt *cell_pkg_index_;
        int *pkg_type_;
        CellNeighborhood *cell_neighborhood_;

        UnsignedInt locateResolutionLevelByPackageType(int target_type, const Vecd &position);
        DataType probeInResolutionLevel(UnsignedInt level, const Vecd &position);
        DataType probeBetweenResolutionLevels(UnsignedInt coarser_level, Real coarser_weight, const Vecd &position);
        DataType probePackageData(const IndexHandler &index_handler, UnsignedInt package_index,
                                  const Arrayi &cell_index, const Vecd &position);
    };

  protected:
    DataContainerUniquePtrAssemble<PackageVariable> pkg_variable_ptrs_;
    PackageVariableAssemble all_pkg_variables_;      /**< all mesh variables on this mesh. */
    PackageVariableAssemble pkg_variables_to_write_; /**< mesh variables to write, which are not empty. */
    PackageVariableAssemble pkh_variables_to_probe_; /**< mesh variables to probe. */
    PackageVariableAssemble evolving_pkg_variables_;
    DataContainerUniquePtrAssemble<MetaVariable> meta_variable_ptrs_;
    MetaVariableAssemble all_meta_variables_;
    MetaVariableAssemble evolving_meta_variables_;
    UniquePtrsKeeper<Entity> unique_variable_ptrs_;
    UnsignedInt num_boundary_pkgs_;        /**< the number boundary packages for each resolution level. */
    StdVec<UnsignedInt> num_pkgs_offsets_; /**< save stating and ending indexes of non-boundary packages. */
    UnsignedInt pkgs_bound_;
    MetaVariable<UnsignedInt> *dv_pkg_1d_cell_index_;               /**< metadata for data packages: cell index. */
    MetaVariable<int> *dv_pkg_type_;                                /**< metadata for data packages: (int)core1/inner0. */
    MetaVariable<CellNeighborhood> *cell_neighborhood_;             /**< 3*3(*3) array to store indicies of neighborhood cells. */
    CellVariable<UnsignedInt> *mcv_cell_pkg_index_;                 /**< the package index for each cell in as 1-d array. */
    ConcurrentVec<std::pair<UnsignedInt, int>> occupied_data_pkgs_; /**< (UnsignedInt)sort_index, (int)core1/inner0. */
    bool is_organized_ = false;                                     /**< whether the data packages are organized. */

    template <typename T>
    T &checkOrganized(std::string func_name, T &value);
    OperationOnDataAssemble<PackageVariableAssemble, PrepareVariablesToWrite<PackageVariable>> sync_mesh_variable_data_{};
    void writePackageVariablesToPltByMesh(UnsignedInt resolution_level, std::ofstream &output_file);
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
