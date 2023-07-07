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
#include "base_variable.h"
#include "my_memory_pool.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <mutex>
using namespace std::placeholders;

namespace SPH
{
/** Iterator on a collection of mesh data packages. sequential computing. */
template <class DataPackageType, typename LocalFunction, typename... Args>
void package_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
                 const LocalFunction &local_function, Args &&...args)
{
    for (size_t i = 0; i != data_pkgs.size(); ++i)
        local_function(data_pkgs[i]);
};
/** Iterator on a collection of mesh data packages. parallel computing. */
template <class DataPackageType, typename LocalFunction, typename... Args>
void package_parallel_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
                          const LocalFunction &local_function, Args &&...args)
{
    parallel_for(
        IndexRange(0, data_pkgs.size()),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                local_function(data_pkgs[i]);
            }
        },
        ap);
};

/**
 * @class BaseDataPackage
 * @brief Abstract base class for a data package,
 * 		  by which the data in a derived class can be on- or off-grid.
 * 		  The data package can be defined in a cell of a background mesh so the pkg_index is
 * 		  the cell location on the mesh.
 */
class BaseDataPackage
{
  public:
    BaseDataPackage() : cell_index_on_mesh_(Arrayi::Zero()), state_indicator_(0){};
    virtual ~BaseDataPackage(){};
    void setInnerPackage() { state_indicator_ = 1; };
    bool isInnerPackage() { return state_indicator_ != 0; };
    void setCorePackage() { state_indicator_ = 2; };
    bool isCorePackage() { return state_indicator_ == 2; };
    void setCellIndexOnMesh(const Arrayi &cell_index) { cell_index_on_mesh_ = cell_index; }
    Arrayi CellIndexOnMesh() const { return cell_index_on_mesh_; }

  protected:
    Arrayi cell_index_on_mesh_; /**< index of this data package on the background mesh, zero if it is not on the mesh. */
    /** reserved value: 0 not occupying background mesh, 1 occupying.
     *  guide to use: larger for high priority of the data package. */
    int state_indicator_;
};

/**
 * @class GridDataPackage
 * @brief Abstract base class for a grid-based data package
 * whose data are defined on the grids of a small mesh patch.
 * Note that, pkg_addrs_size = pkg_size + 2 * pkg_addrs_buffer;
 * Also note that, while the mesh_lower_bound_ locates the first data address,
 * the DataLowerBound() locates the first data for initialization.
 * Also note that, as a inner package is contained in a background mesh cell,
 * it will be constructed based on the grid position (lower-left corner) of the cell.
 */
template <int PKG_SIZE, int ADDRS_BUFFER>
class GridDataPackage : public BaseDataPackage, public BaseMesh
{
  public:
    static constexpr int pkg_size = PKG_SIZE;
    static constexpr int pkg_addrs_size = PKG_SIZE + ADDRS_BUFFER * 2;
    static constexpr int pkg_addrs_buffer = ADDRS_BUFFER;
    static constexpr int pkg_ops_end = PKG_SIZE + pkg_addrs_buffer;
    template <typename DataType>
    using PackageData = PackageDataMatrix<DataType, PKG_SIZE>;
    template <typename DataType>
    using PackageDataAddress = PackageDataMatrix<DataType *, pkg_addrs_size>;
    /** Matrix data for temporary usage. Note that it is array with pkg_addrs_size.  */
    template <typename DataType>
    using PackageTemporaryData = PackageDataMatrix<DataType, pkg_addrs_size>;

    /** Default constructor for singular package */
    GridDataPackage() : BaseDataPackage(), BaseMesh(pkg_addrs_size * Arrayi::Ones()){};
    /** Constructor for inner package */
    GridDataPackage(const Vecd &container_lower_bound, Real data_spacing)
        : BaseDataPackage(),
          BaseMesh(container_lower_bound - data_spacing * Vecd::Ones() * ((Real)pkg_addrs_buffer - 0.5),
                   data_spacing, pkg_addrs_size * Arrayi::Ones()){};
    virtual ~GridDataPackage(){};
    Vecd DataPositionFromIndex(const Vecd &data_index) { return DataLowerBound() + data_index * grid_spacing_; };
    /** void (non_value_returning) function iterate on all data points by value,
     *  for function only involving the data itself. */
    template <typename FunctionOnData>
    void for_each_data(const FunctionOnData &function);
    /** void (non_value_returning) function iterate on all data points by address,
     *  for function involving operations on data neighbors. */
    template <typename FunctionOnAddress>
    void for_each_addrs(const FunctionOnAddress &function);
    /** access specific package data with mesh variable */
    template <typename DataType>
    PackageData<DataType> &getPackageData(const MeshVariable<DataType> &mesh_variable)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        return std::get<type_index>(all_pkg_data_)[mesh_variable.IndexInContainer()];
    };
    /** access specific package data address with mesh variable */
    template <typename DataType>
    PackageDataAddress<DataType> &getPackageDataAddress(const MeshVariable<DataType> &mesh_variable)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        return std::get<type_index>(all_pkg_data_addrs_)[mesh_variable.IndexInContainer()];
    };
    /** probe by applying bi and tri-linear interpolation within the package. */
    template <typename DataType>
    DataType probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position);
    /** assign value to data package according to the position of data */
    template <typename DataType, typename FunctionByPosition>
    void assignByPosition(const MeshVariable<DataType> &mesh_variable,
                          const FunctionByPosition &function_by_position);
    /** compute gradient transform within data package */
    template <typename InDataType, typename OutDataType>
    void computeGradient(const MeshVariable<InDataType> &in_variable,
                         const MeshVariable<OutDataType> &out_variable);
    /** obtain averaged value at a corner of a data cell */
    template <typename DataType>
    DataType CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs,
                           Arrayi addrs_index, Arrayi corner_direction);

  protected:
    DataContainerAssemble<PackageData> all_pkg_data_;
    DataContainerAssemble<PackageDataAddress> all_pkg_data_addrs_;

    /** lower bound coordinate for the data as reference */
    Vecd DataLowerBound() { return mesh_lower_bound_ + grid_spacing_ * Vecd::Ones() * (Real)pkg_addrs_buffer; };
    /** allocate memory for all package data */
    template <typename DataType>
    struct AllVariablesAllocation
    {
        void operator()(DataContainerAssemble<PackageData> &all_pkg_data,
                        DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs,
                        const MeshVariableAssemble &all_mesh_variables_)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            size_t total_variables = std::get<type_index>(all_mesh_variables_).size();
            std::get<type_index>(all_pkg_data).resize(total_variables);
            std::get<type_index>(all_pkg_data_addrs).resize(total_variables);
        };
    };
    DataAssembleOperation<AllVariablesAllocation> allocate_all_variables_;
    /** set the initial package data address for singular data package */
    template <typename DataType>
    struct AssignSingularPackageDataAddress
    {
        void operator()(DataContainerAssemble<PackageData> &all_pkg_data,
                        DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs);
    };
    DataAssembleOperation<AssignSingularPackageDataAddress> assign_singular_pkg_data_addrs_;
    /** assign address for all package data when the package is an inner one */
    template <typename DataType>
    struct AssignPackageDataAddress
    {
        void operator()(DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs,
                        const Arrayi &addrs_index,
                        DataContainerAssemble<PackageData> &all_pkg_data,
                        const Arrayi &data_index);
    };
    DataAssembleOperation<AssignPackageDataAddress> assign_pkg_data_addrs_;

  public:
    void allocateAllVariables(const MeshVariableAssemble &all_mesh_variables_)
    {
        allocate_all_variables_(all_pkg_data_, all_pkg_data_addrs_, all_mesh_variables_);
    };

    void assignSingularPackageDataAddress()
    {
        assign_singular_pkg_data_addrs_(all_pkg_data_, all_pkg_data_addrs_);
    };

    void assignPackageDataAddress(const Arrayi &addrs_index, GridDataPackage *src_pkg, const Arrayi &data_index)
    {
        assign_pkg_data_addrs_(all_pkg_data_addrs_, addrs_index, src_pkg->all_pkg_data_, data_index);
    };
};

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
template <class GridDataPackageType>
class MeshWithGridDataPackages : public Mesh
{
  private:
    DataContainerUniquePtrAssemble<MeshVariable> mesh_variable_ptrs_;

  public:
    template <typename... Args>
    explicit MeshWithGridDataPackages(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size)
        : Mesh(tentative_bounds, GridDataPackageType::pkg_size * data_spacing, buffer_size),
          data_spacing_(data_spacing),
          global_mesh_(this->mesh_lower_bound_ + 0.5 * data_spacing * Vecd::Ones(), data_spacing, this->all_cells_ * pkg_size)
    {
        allocateMeshDataMatrix();
    };
    virtual ~MeshWithGridDataPackages() { deleteMeshDataMatrix(); };
    /** spacing between the data, which is 1/ pkg_size of this grid spacing */
    virtual Real DataSpacing() override { return data_spacing_; };

  protected:
    MeshVariableAssemble all_mesh_variables_;              /**< all mesh variables on this mesh. */
    MyMemoryPool<GridDataPackageType> data_pkg_pool_;      /**< memory pool for all packages in the mesh. */
    MeshDataMatrix<GridDataPackageType *> data_pkg_addrs_; /**< Address of data packages. */
    ConcurrentVec<GridDataPackageType *> inner_data_pkgs_; /**< Inner data packages which is able to carry out spatial operations. */
    /** Singular data packages. provided for far field condition with usually only two values.
     * For example, when level set is considered. The first value for inner far-field and second for outer far-field */
    StdVec<GridDataPackageType *> singular_data_pkgs_addrs_;
    static constexpr int pkg_size = GridDataPackageType::pkg_size;                 /**< the size of the data package matrix*/
    static constexpr int pkg_addrs_buffer = GridDataPackageType::pkg_addrs_buffer; /**< the size of address buffer, a value less than the package size. */
    static constexpr int pkg_ops_end = GridDataPackageType::pkg_ops_end;           /**< the size of operation loops. */
    static constexpr int pkg_addrs_size = GridDataPackageType::pkg_addrs_size;     /**< the size of address matrix in the data packages. */
    const Real data_spacing_;                                                      /**< spacing of data in the data packages*/
    std::mutex mutex_my_pool;                                                      /**< mutex exclusion for memory pool */
    BaseMesh global_mesh_;                                                         /**< the mesh for the locations of all possible data points. */

    void allocateMeshDataMatrix(); /**< allocate memories for addresses of data packages. */
    void deleteMeshDataMatrix();   /**< delete memories for addresses of data packages. */

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

    template <typename InitializeSingularData>
    void initializeASingularDataPackage(
        const DataContainerAddressAssemble<MeshVariable> &all_mesh_variables_,
        const InitializeSingularData &initialize_singular_data)
    {
        GridDataPackageType *new_data_pkg = data_pkg_pool_.malloc();
        new_data_pkg->allocateAllVariables(all_mesh_variables_);
        initialize_singular_data(new_data_pkg);
        new_data_pkg->assignSingularPackageDataAddress();
        singular_data_pkgs_addrs_.push_back(new_data_pkg);
    };

    template <typename InitializePackageData>
    GridDataPackageType *createDataPackage(
        const DataContainerAddressAssemble<MeshVariable> &all_mesh_variables_,
        const Arrayi &cell_index,
        const InitializePackageData &initialize_package_data)
    {
        mutex_my_pool.lock();
        Vecd cell_position = CellPositionFromIndex(cell_index);
        Vecd grid_position = GridPositionFromCellPosition(cell_position);
        GridDataPackageType *new_data_pkg = data_pkg_pool_.malloc(grid_position, data_spacing_);
        mutex_my_pool.unlock();
        new_data_pkg->allocateAllVariables(all_mesh_variables_);
        initialize_package_data(new_data_pkg);
        new_data_pkg->setCellIndexOnMesh(cell_index);
        assignDataPackageAddress(cell_index, new_data_pkg);
        return new_data_pkg;
    };

    void assignDataPackageAddress(const Arrayi &cell_index, GridDataPackageType *data_pkg);
    /** Return data package with given cell index. */
    GridDataPackageType *DataPackageFromCellIndex(const Arrayi &cell_index);
    void initializePackageAddressesInACell(const Arrayi &cell_index);
    /** Find related cell index and data index for a data package address matrix */
    std::pair<int, int> CellShiftAndDataIndex(int data_addrs_index_component)
    {
        std::pair<int, int> shift_and_index;
        int signed_date_index = data_addrs_index_component - pkg_addrs_buffer;
        shift_and_index.first = (signed_date_index + pkg_size) / pkg_size - 1;
        shift_and_index.second = signed_date_index - shift_and_index.first * pkg_size;
        return shift_and_index;
    }
    /** This function probe a mesh value */
    template <class DataType>
    DataType probeMesh(const MeshVariable<DataType> &mesh_variable, const Vecd &position);
    /** This function find the value of data from its index from global mesh. */
    template <typename DataType>
    DataType DataValueFromGlobalIndex(const MeshVariable<DataType> &mesh_variable,
                                      const Arrayi &global_grid_index);
};
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_H
