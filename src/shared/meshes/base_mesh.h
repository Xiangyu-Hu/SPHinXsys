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

#include "base_data_type_package.h"
#include "sphinxsys_variable.h"

#include <fstream>

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
    Mesh(BoundingBoxd tentative_bounds, Real grid_spacing,
         UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset = 0);
    Mesh(Vecd mesh_lower_bound, Real grid_spacing, Arrayi all_grid_points);
    ~Mesh() {};

    Vecd MeshLowerBound() const { return mesh_lower_bound_; };
    Real GridSpacing() const { return grid_spacing_; };
    UnsignedInt BufferWidth() { return buffer_width_; };
    Arrayi AllGridPoints() const { return all_grid_points_; };
    Arrayi AllCells() const { return all_cells_; };
    UnsignedInt NumberOfGridPoints() const { return all_grid_points_.prod(); };
    UnsignedInt NumberOfCells() const { return all_cells_.prod(); };

    Arrayi CellIndexFromPosition(const Vecd &position) const;
    UnsignedInt LinearCellIndexFromPosition(const Vecd &position) const;
    UnsignedInt LinearCellIndex(const Arrayi &cell_index) const;
    Arrayi DimensionalCellIndex(UnsignedInt linear_index) const;
    Vecd CellPositionFromIndex(const Arrayi &cell_index) const;
    Vecd GridPositionFromIndex(const Arrayi &grid_index) const;
    Vecd CellLowerCornerPosition(const Arrayi &cell_index) const;
    Arrayi boundCellIndex(const Arrayi &input) const;
    //----------------------------------------------------------------------
    // Transferring between 1D mesh indexes.
    // Here, mesh size can be either AllGridPoints or AllCells.
    //----------------------------------------------------------------------
    static Array2i transfer1DtoMeshIndex(const Array2i &mesh_size, UnsignedInt i);
    static Array3i transfer1DtoMeshIndex(const Array3i &mesh_size, UnsignedInt i);
    static UnsignedInt transferMeshIndexTo1D(const Array2i &mesh_size, const Array2i &mesh_index);
    static UnsignedInt transferMeshIndexTo1D(const Array3i &mesh_size, const Array3i &mesh_index);

    /** converts mesh index into a Morton order.
     * Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros
     * https://stackoverflow.com/questions/18529057/
     * produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
     */
    static UnsignedInt transferMeshIndexToMortonOrder(const Array2i &mesh_index);
    static UnsignedInt transferMeshIndexToMortonOrder(const Array3i &mesh_index);

  protected:
    Vecd mesh_lower_bound_;                /**< mesh lower bound as reference coordinate */
    Real grid_spacing_;                    /**< grid_spacing */
    UnsignedInt buffer_width_;             /**< buffer width to avoid bound check.*/
    Arrayi all_grid_points_;               /**< number of grid points by dimension */
    Arrayi all_cells_;                     /**< number of cells by dimension */
    UnsignedInt linear_cell_index_offset_; /**< offset for linear cell index, used for sub-mesh */

    static UnsignedInt MortonCode(const UnsignedInt &i);
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
    explicit BaseMeshField(const std::string &name) : name_(name) {};
    virtual ~BaseMeshField() {};
    /** Return the mesh field name. */
    std::string Name() { return name_; };
    void setName(const std::string &new_name) { name_ = new_name; };
    /** output mesh data for Tecplot visualization */
    virtual void writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence = 0) = 0;
    /** output mesh data for Tecplot visualization */
    virtual void writeBKGMeshToPlt(const std::string &partial_file_name) {};
};

class MultiLevelMeshField : public BaseMeshField
{
    typedef DataContainerAddressAssemble<DiscreteVariable> CellVariableAssemble;
    DataContainerUniquePtrAssemble<DiscreteVariable> cell_variable_ptrs_;

  public:
    MultiLevelMeshField(
        const std::string &name, BoundingBoxd tentative_bounds,
        Real Reference_grid_spacing, UnsignedInt buffer_width, size_t total_levels = 1);
    virtual ~MultiLevelMeshField() {};

    UniquePtrsKeeper<Mesh> mesh_ptrs_keeper_;
    StdVec<Mesh *> &getMeshes() { return meshes_; };
    UnsignedInt TotalNumberOfCells() { return total_number_of_cells_; };

    template <typename DataType, typename... Args>
    DiscreteVariable<DataType> *registerCellVariable(const std::string &variable_name, Args &&...args);
    template <typename DataType>
    DiscreteVariable<DataType> *getCellVariable(const std::string &variable_name);

    template <typename DataType>
    void addCellVariableToWrite(const std::string &variable_name);
    void writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence = 0) override;

    template <class ExecutionPolicy>
    void syncCellVariablesToWrite(ExecutionPolicy &ex_policy);

  protected:
    size_t total_levels_; /**< level 0 is the coarsest */
    StdVec<Mesh *> meshes_;
    UnsignedInt total_number_of_cells_;
    CellVariableAssemble all_cell_variables_;
    CellVariableAssemble cell_variables_to_write_;

    void writeCellVariableToPltByMesh(const Mesh &mesh, std::ofstream &output_file);

  protected:
    template <template <typename> class MeshVariableType>
    struct SyncMeshVariableData
    {
        template <typename DataType, class ExecutionPolicy>
        void operator()(DataContainerAddressKeeper<MeshVariableType<DataType>> &mesh_variables,
                        ExecutionPolicy &ex_policy);
    };

    OperationOnDataAssemble<CellVariableAssemble, SyncMeshVariableData<DiscreteVariable>> sync_cell_variable_data_{};
};
} // namespace SPH
#endif // BASE_MESH_H