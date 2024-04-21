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
 * @file 	cell_linked_list.h
 * @brief 	Here gives the classes for managing cell linked lists. This is the basic class
 * 			for building the particle configurations.
 * @details The cell linked list saves for each body a list of particles
 * 			located within the cell.
 * @author	Chi Zhang, Yongchuan and Xiangyu Hu
 */

#ifndef MESH_CELL_LINKED_LIST_H
#define MESH_CELL_LINKED_LIST_H

#include "base_mesh.h"
#include "neighborhood.h"
#include "mesh_iterators.hpp"

namespace SPH
{

class BaseParticles;
class Kernel;
class SPHAdaptation;
class CellLinkedList;

/**
 * @class BaseCellLinkedList
 * @brief The Abstract class for mesh cell linked list derived from BaseMeshField.
 */
class BaseCellLinkedList : public BaseMeshField
{
  protected:
    Kernel &kernel_;

    /** clear split cell lists in this mesh*/
    virtual void clearSplitCellLists(SplitCellLists &split_cell_lists);
    /** update split particle list in this mesh */
    virtual void updateSplitCellLists(SplitCellLists &split_cell_lists) = 0;

  public:
    BaseCellLinkedList(SPHAdaptation &sph_adaptation);
    virtual ~BaseCellLinkedList(){};

    /** access concrete cell linked list levels*/
    virtual StdVec<CellLinkedList *> CellLinkedListLevels() = 0;
    virtual execution::ExecutionEvent UpdateCellLists(BaseParticles &base_particles) = 0;
    virtual SplitCellLists *getSplitCellLists();
    virtual void setUseSplitCellLists();
    /** Insert a cell-linked_list entry to the concurrent index list. */
    virtual void insertParticleIndex(size_t particle_index, const Vecd &particle_position) = 0;
    /** Insert a cell-linked_list entry of the index and particle position pair. */
    virtual void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) = 0;
    /** find the nearest list data entry */
    virtual ListData findNearestListDataEntry(const Vecd &position) = 0;
    /** computing the sequence which indicate the order of sorted particle data */
    virtual size_t* computingSequence(BaseParticles &base_particles) = 0;
    /** Tag body part by cell, call by body part */
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) = 0;
    /** Tag domain bounding cells in an axis direction, called by domain bounding classes */
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) = 0;
};


class CellLinkedListKernel {
  public:
    CellLinkedListKernel(const DeviceVecd &meshLowerBound, DeviceReal gridSpacing,
                         const DeviceArrayi &allGridPoints, const DeviceArrayi &allCells);

    ~CellLinkedListKernel();

    execution::ExecutionEvent clearCellLists();
    execution::ExecutionEvent UpdateCellLists(const BaseParticles &base_particles);

    template <typename FunctionOnEach>
    void forEachNeighbor(size_t index_i, const DeviceVecd *self_position,
                         const FunctionOnEach &function) const
    {
        const DeviceVecd pos_i = self_position[index_i];
        const size_t search_depth = 1;
        const auto target_cell_index = CellIndexFromPosition(pos_i, *mesh_lower_bound_,
                                                             *grid_spacing_, *all_grid_points_);
        mesh_for_each_array(
            VecdMax(DeviceArrayi{0}, DeviceArrayi{target_cell_index - search_depth}),
            VecdMin(*all_cells_, DeviceArrayi{target_cell_index + search_depth + 1}),
            [&](DeviceArrayi&& cell_index) {
                const auto linear_cell_index = transferCellIndexTo1D(cell_index, *all_cells_);
                size_t index_j = index_head_list_[linear_cell_index];
                // Cell list ends when index_j == 0, if index_j is already zero then cell is empty.
                while(index_j--) {  // abbreviates while(index_j != 0) { index_j -= 1; ... }
                    function(pos_i, index_j, list_data_pos_[index_j]);
                    index_j = index_list_[index_j];
                }
            });
    }

    template <typename FunctionOnEach>
    void forEachInnerNeighbor(size_t index_i, const FunctionOnEach &function) const
    {
        forEachNeighbor(index_i, list_data_pos_, function);
    }


    size_t* computingSequence(BaseParticles &baseParticles);

    template<class Type, int Dim>
    static inline DeviceArrayi CellIndexFromPosition(const sycl::vec<Type,Dim> &position, const sycl::vec<Type,Dim>& mesh_lower_bound,
                                                     const Type &grid_spacing, const sycl::vec<int,Dim> &all_grid_points)
    {
        return sycl::min(
            sycl::max(
                sycl::floor((position - mesh_lower_bound) / grid_spacing)
                    .template convert<int, sycl::rounding_mode::rtz>(), DeviceArrayi{0}),
            all_grid_points - DeviceArrayi{2});
    }

    template<class Type, int Dim>
    static inline DeviceArrayi CellIndexFromPosition(const Eigen::Matrix<Type,Dim,1> &position,
                                                     const Eigen::Matrix<Type,Dim,1>& mesh_lower_bound,
                                                     const Type &grid_spacing, const Eigen::Array<int,Dim,1> &all_grid_points)
    {
        DeviceArrayi pos_floor = floor((position - mesh_lower_bound).array() / grid_spacing).template cast<int>();
        return VecdMin(VecdMax(pos_floor, DeviceArrayi{0}), DeviceArrayi{all_grid_points - DeviceArrayi{2}});
    }

    static inline size_t transferCellIndexTo1D(const DeviceArray2i &cell_index, const DeviceArray2i &all_cells)
    {
        return cell_index[0] * all_cells[1] + cell_index[1];
    }

    static inline size_t transferCellIndexTo1D(const DeviceArray3i &cell_index, const DeviceArray3i &all_cells)
    {
        return cell_index[0] * all_cells[1] * all_cells[2] +
               cell_index[1] * all_cells[2] +
               cell_index[2];
    }

  private:
    size_t total_real_particles_;
    DeviceVecd* list_data_pos_;

    DeviceVecd *mesh_lower_bound_;
    DeviceReal *grid_spacing_;
    DeviceArrayi *all_grid_points_, *all_cells_;

    size_t* index_list_;
    size_t* index_head_list_;
    size_t index_head_list_size_;
};


/**
 * @class CellLinkedList
 * @brief Defining a mesh cell linked list for a body.
 * 		  The meshes for all bodies share the same global coordinates.
 */
class CellLinkedList : public BaseCellLinkedList, public Mesh
{
    StdVec<CellLinkedList *> single_cell_linked_list_level_;
    /**
     * @brief particle by cells lists is for parallel splitting algorithm.
     * All particles in each cell are collected together.
     * If two particles each belongs two different cell entries,
     * they have no interaction because they are too far.
     */
    SplitCellLists split_cell_lists_;
    bool use_split_cell_lists_;

  protected:
    /** using concurrent vectors due to writing conflicts when building the list */
    MeshDataMatrix<ConcurrentIndexVector> cell_index_lists_;
    /** non-concurrent list data rewritten for building neighbor list */
    MeshDataMatrix<ListDataVector> cell_data_lists_;

    void allocateMeshDataMatrix(); /**< allocate memories for addresses of data packages. */
    void deleteMeshDataMatrix();   /**< delete memories for addresses of data packages. */
    virtual void updateSplitCellLists(SplitCellLists &split_cell_lists) override;

  public:
    CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing, SPHAdaptation &sph_adaptation);
    virtual ~CellLinkedList() { deleteMeshDataMatrix(); };

    void clearCellLists();
    virtual SplitCellLists *getSplitCellLists() override { return &split_cell_lists_; };
    virtual void setUseSplitCellLists() override { use_split_cell_lists_ = true; };
    void UpdateCellListData(BaseParticles &base_particles);
    virtual execution::ExecutionEvent UpdateCellLists(BaseParticles &base_particles) override;
    void insertParticleIndex(size_t particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) override;
    virtual ListData findNearestListDataEntry(const Vecd &position) override;
    virtual size_t* computingSequence(BaseParticles &base_particles) override;
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) override;
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) override;
    virtual void writeMeshFieldToPlt(std::ofstream &output_file) override;
    virtual StdVec<CellLinkedList *> CellLinkedListLevels() override { return single_cell_linked_list_level_; };

    /** generalized particle search algorithm */
    template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
    void searchNeighborsByParticles(DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
                                    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation);

    execution::DeviceImplementation<CellLinkedListKernel> device_kernel;
};

/**
 * @class MultilevelCellLinkedList
 * @brief Defining a multilevel mesh cell linked list for a body
 * 		  for multi-resolution particle configuration.
 */
class MultilevelCellLinkedList : public MultilevelMesh<BaseCellLinkedList, CellLinkedList, RefinedMesh<CellLinkedList>>
{
  protected:
    StdLargeVec<Real> &h_ratio_; /**< Smoothing length for each level. */
    /** Update split cell list. */
    virtual void updateSplitCellLists(SplitCellLists &split_cell_lists) override{};
    /** determine mesh level from particle cutoff radius */
    inline size_t getMeshLevel(Real particle_cutoff_radius);

  public:
    MultilevelCellLinkedList(BoundingBox tentative_bounds, Real reference_grid_spacing,
                             size_t total_levels, SPHAdaptation &sph_adaptation);
    virtual ~MultilevelCellLinkedList(){};

    virtual execution::ExecutionEvent UpdateCellLists(BaseParticles &base_particles) override;
    void insertParticleIndex(size_t particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) override;
    virtual ListData findNearestListDataEntry(const Vecd &position) override { return ListData(0, Vecd::Zero()); }; // mocking, not implemented
    virtual size_t* computingSequence(BaseParticles &base_particles) override;
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) override;
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) override{};
    virtual StdVec<CellLinkedList *> CellLinkedListLevels() override { return getMeshLevels(); };
};
} // namespace SPH
#endif // MESH_CELL_LINKED_LIST_H