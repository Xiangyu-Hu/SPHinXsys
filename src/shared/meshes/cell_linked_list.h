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
 * @file 	cell_linked_list.h
 * @brief 	Here gives the classes for managing cell linked lists. This is the basic class
 * 			for building the particle configurations.
 * @details The cell linked list saves for each body a list of particles
 * 			located within the cell.
 * @author	Chi Zhang, Yongchuan and Xiangyu Hu
 */

#ifndef MESH_CELL_LINKED_LIST_H
#define MESH_CELL_LINKED_LIST_H

#include "base_mesh.hpp"
#include "execution_policy.h"
#include "neighborhood.h"

namespace SPH
{

class BaseParticles;
class Kernel;
class SPHAdaptation;
class AdaptiveSmoothingLength;
class AnisotropicAdaptation;

/**
 * @class BaseCellLinkedList
 * @brief The Abstract class for mesh cell linked list derived from BaseMeshField.
 */
class BaseCellLinkedList : public MultiResolutionMeshField<Mesh>
{
  protected:
    BaseParticles &base_particles_;

  public:
    BaseCellLinkedList(BaseParticles &base_particles, SPHAdaptation &sph_adaptation,
                       BoundingBoxd tentative_bounds, Real Reference_grid_spacing, size_t total_levels);
    virtual ~BaseCellLinkedList() {};
    BaseParticles &getBaseParticles() { return base_particles_; };
    void UpdateCellLists(BaseParticles &base_particles);
    /** Insert a cell-linked_list entry to the concurrent index list. */
    virtual void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) = 0;
    /** Insert a cell-linked_list entry of the index and particle position pair. */
    virtual void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) = 0;
    /** find the nearest list data entry */
    ListData findNearestListDataEntry(const Vecd &position);
    /** computing the sequence which indicate the order of sorted particle data */
    UnsignedInt computingSequence(Vecd &position, UnsignedInt index_i);
    /** Tag body part by cell, call by body part */
    void tagBodyPartByCell(ConcurrentCellLists &cell_lists,
                           ConcurrentIndexVector &cell_indexes,
                           std::function<bool(Vecd, Real)> &check_included);
    /** Tag domain bounding cells in an axis direction, called by domain bounding classes */
    void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBoxd &bounding_bounds, int axis);

    /** split algorithm */;
    template <class ExecutionPolicy, class LocalDynamicsFunction>
    void particle_for_split(const ExecutionPolicy &ex_policy, const LocalDynamicsFunction &local_dynamics_function);

    /** generalized particle search algorithm */
    template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
    void searchNeighborsByMesh(Mesh &mesh, DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
                               GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation);
    DiscreteVariable<UnsignedInt> *dvParticleIndex() { return dv_particle_index_; };
    DiscreteVariable<UnsignedInt> *dvCellOffset() { return dv_cell_offset_; };
    Mesh getSortSequenceMesh();

  protected:
    Kernel &kernel_;
    DiscreteVariable<UnsignedInt> *dv_particle_index_;
    DiscreteVariable<UnsignedInt> *dv_cell_offset_;
    /** using concurrent vectors due to writing conflicts when building the list */
    StdVec<ConcurrentIndexVector> cell_index_lists_;
    /** non-concurrent list data rewritten for building neighbor list */
    StdVec<ListDataVector> cell_data_lists_;

    void clearCellLists();
    void UpdateCellListData(BaseParticles &base_particles);
    void tagBodyPartByCellByMesh(Mesh &mesh, ConcurrentCellLists &cell_lists,
                                 std::function<bool(Vecd, Real)> &check_included);
    virtual void tagBodyPartByCellCK(ConcurrentIndexVector &cell_indexes,
                                     std::function<bool(Vecd, Real)> &check_included) = 0;
    void tagBoundingCellsByMesh(Mesh &mesh, StdVec<CellLists> &cell_data_lists,
                                const BoundingBoxd &bounding_bounds, int axis);
    void findNearestListDataEntryByMesh(Mesh &mesh, Real &min_distance_sqr, ListData &nearest_entry,
                                        const Vecd &position);
    /** split algorithm */;
    template <class ExecutionPolicy, class LocalDynamicsFunction>
    void particle_for_split_by_mesh(const ExecutionPolicy &ex_policy, Mesh &mesh,
                                    const LocalDynamicsFunction &local_dynamics_function);
};

template <typename...>
class CellLinkedList;

template <>
class CellLinkedList<SPHAdaptation> : public BaseCellLinkedList
{
  public:
    typedef Mesh CellLinkedListMeshType;

    CellLinkedList(BoundingBoxd tentative_bounds, Real grid_spacing,
                   BaseParticles &base_particles, SPHAdaptation &sph_adaptation);
    ~CellLinkedList() {};
    Mesh &getMesh() { return *mesh_; };
    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override;
    virtual void tagBodyPartByCellCK(ConcurrentIndexVector &cell_indexes,
                                     std::function<bool(Vecd, Real)> &check_included) override;

    class NeighborSearch : public Mesh
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        NeighborSearch(const ExecutionPolicy &ex_policy, Encloser &encloser);

        template <typename FunctionOnEach>
        void forInnerSearch(const Vecd &source_pos, const FunctionOnEach &function, const Vecd &src_cut_off) const;
        template <typename FunctionOnEach>
        void forContactSearch(const Vecd &source_pos, const FunctionOnEach &function, const Vecd &src_cut_off) const;

      protected:
        UnsignedInt *particle_index_;
        UnsignedInt *cell_offset_;

        inline BoundingBoxi InnerSearchBox(const Vecd &src_cut_off) const;
        inline BoundingBoxi ContactSearchBox(const Vecd &src_cut_off) const;
        template <typename FunctionOnEach>
        void searchInRange(const FunctionOnEach &function, const BoundingBoxi &rang_box) const;
    };

    Mesh &getCellLinkedListMesh() { return cell_linked_list_mesh_; };

  protected:
    Mesh *mesh_;
    Mesh cell_linked_list_mesh_;
};

template <>
class CellLinkedList<AdaptiveSmoothingLength> : public BaseCellLinkedList
{
  protected:
    Real *h_ratio_; /**< Smoothing length for each level. */
    int *h_level_;  /**< Smoothing length level for each particle. */

    /** determine mesh level from particle cutoff radius */
    inline UnsignedInt getMeshLevel(Real particle_cutoff_radius);

  public:
    CellLinkedList(BoundingBoxd tentative_bounds,
                   Real reference_grid_spacing, UnsignedInt total_levels,
                   BaseParticles &base_particles, SPHAdaptation &sph_adaptation);
    virtual ~CellLinkedList() {};
    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override;
    virtual void tagBodyPartByCellCK(ConcurrentIndexVector &cell_indexes,
                                     std::function<bool(Vecd, Real)> &check_included) override;

    class CellLinkedListMesh : public Mesh
    {
      public:
        CellLinkedListMesh(BaseCellLinkedList &base_cell_linked_list);
        Real CoarsestGridSpacing() const { return coarsest_grid_spacing_; };

      protected:
        Real coarsest_grid_spacing_;
    };

    typedef CellLinkedListMesh CellLinkedListMeshType;

    class NeighborSearch : public CellLinkedListMesh
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        NeighborSearch(const ExecutionPolicy &ex_policy, Encloser &encloser);

        template <typename FunctionOnEach>
        void forInnerSearch(const Vecd &source_pos, const FunctionOnEach &function, const Vecd &src_cut_off) const;
        template <typename FunctionOnEach>
        void forContactSearch(const Vecd &source_pos, const FunctionOnEach &function, const Vecd &src_cut_off) const;

      protected:
        UnsignedInt *particle_index_;
        UnsignedInt *cell_offset_;

        inline BoundingBoxi InnerSearchBox(const Vecd &src_cut_off) const;
        inline BoundingBoxi ContactSearchBox(const Vecd &src_cut_off) const;
        template <typename FunctionOnEach>
        void searchInRange(const FunctionOnEach &function, const BoundingBoxi &rang_box) const;
    };

    CellLinkedListMesh &getCellLinkedListMesh() { return cell_linked_list_mesh_; };

  protected:
    CellLinkedListMesh cell_linked_list_mesh_;
};

template <>
class CellLinkedList<AnisotropicAdaptation> : public BaseCellLinkedList
{

  public:
    CellLinkedList(BoundingBoxd tentative_bounds, Real grid_spacing,
                   BaseParticles &base_particles, SPHAdaptation &sph_adaptation);

  protected:
    Mesh *mesh_;
    Mesh cell_linked_list_mesh_;
};
} // namespace SPH
#endif // MESH_CELL_LINKED_LIST_H
