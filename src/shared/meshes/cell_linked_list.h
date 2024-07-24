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
    virtual void UpdateCellLists(BaseParticles &base_particles) = 0;
    virtual SplitCellLists *getSplitCellLists();
    virtual void setUseSplitCellLists();
    /** Insert a cell-linked_list entry to the concurrent index list. */
    virtual void insertParticleIndex(size_t particle_index, const Vecd &particle_position) = 0;
    /** Insert a cell-linked_list entry of the index and particle position pair. */
    virtual void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) = 0;
    /** find the nearest list data entry */
    virtual ListData findNearestListDataEntry(const Vecd &position) = 0;
    /** computing the sequence which indicate the order of sorted particle data */
    virtual StdLargeVec<size_t> &computingSequence(BaseParticles &base_particles) = 0;
    /** Tag body part by cell, call by body part */
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) = 0;
    /** Tag domain bounding cells in an axis direction, called by domain bounding classes */
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) = 0;
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
    ConcurrentIndexVector *cell_index_lists_;
    /** non-concurrent list data rewritten for building neighbor list */
    ListDataVector *cell_data_lists_;

    void allocateMeshDataMatrix(); /**< allocate memories for addresses of data packages. */
    void deleteMeshDataMatrix();   /**< delete memories for addresses of data packages. */
    virtual void updateSplitCellLists(SplitCellLists &split_cell_lists) override;
    template <typename DataListsType>
    DataListsType &getCellDataList(DataListsType *data_lists, const Arrayi &cell_index)
    {
        return data_lists[transferMeshIndexTo1D(all_cells_, cell_index)];
    };

  public:
    CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing, SPHAdaptation &sph_adaptation);
    ~CellLinkedList() { deleteMeshDataMatrix(); };

    void clearCellLists();
    virtual SplitCellLists *getSplitCellLists() override { return &split_cell_lists_; };
    virtual void setUseSplitCellLists() override { use_split_cell_lists_ = true; };
    void UpdateCellListData(BaseParticles &base_particles);
    virtual void UpdateCellLists(BaseParticles &base_particles) override;
    void insertParticleIndex(size_t particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) override;
    virtual ListData findNearestListDataEntry(const Vecd &position) override;
    virtual StdLargeVec<size_t> &computingSequence(BaseParticles &base_particles) override;
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) override;
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) override;
    virtual void writeMeshFieldToPlt(std::ofstream &output_file) override;
    virtual StdVec<CellLinkedList *> CellLinkedListLevels() override { return single_cell_linked_list_level_; };

    /** generalized particle search algorithm */
    template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
    void searchNeighborsByParticles(DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
                                    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation);
};

template <>
class RefinedMesh<CellLinkedList> : public CellLinkedList
{
  public:
    RefinedMesh(BoundingBox tentative_bounds, CellLinkedList &coarse_mesh, SPHAdaptation &sph_adaptation)
        : CellLinkedList(tentative_bounds, 0.5 * coarse_mesh.GridSpacing(), sph_adaptation){};
};

/**
 * @class MultilevelCellLinkedList
 * @brief Defining a multilevel mesh cell linked list for a body
 * 		  for multi-resolution particle configuration.
 */
class MultilevelCellLinkedList : public MultilevelMesh<BaseCellLinkedList, CellLinkedList>
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

    virtual void UpdateCellLists(BaseParticles &base_particles) override;
    void insertParticleIndex(size_t particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(size_t particle_index, const Vecd &particle_position) override;
    virtual ListData findNearestListDataEntry(const Vecd &position) override { return ListData(0, Vecd::Zero()); }; // mocking, not implemented
    virtual StdLargeVec<size_t> &computingSequence(BaseParticles &base_particles) override;
    virtual void tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included) override;
    virtual void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis) override{};
    virtual StdVec<CellLinkedList *> CellLinkedListLevels() override { return getMeshLevels(); };
};
} // namespace SPH
#endif // MESH_CELL_LINKED_LIST_H