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

#include "base_mesh.h"
#include "execution_policy.h"
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
    BaseParticles &base_particles_;
    DataContainerUniquePtrAssemble<DiscreteVariable> all_discrete_variable_ptrs_;
    UniquePtrsKeeper<Entity> unique_variable_ptrs_;
    UniquePtrsKeeper<SingularVariable<Mesh>> mesh_ptrs_keeper_;
    StdVec<Mesh *> meshes_;
    Mesh *coarsest_mesh_;
    Mesh *finest_mesh_;

  public:
    BaseCellLinkedList(BaseParticles &base_particles, SPHAdaptation &sph_adaptation);
    virtual ~BaseCellLinkedList();
    StdVec<Mesh *> &getMeshes() { return meshes_; };
    BaseParticles &getBaseParticles() { return base_particles_; };
    void writeMeshFieldToPlt(const std::string &partial_file_name);
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
    void tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis);

    /** split algorithm */;
    template <class ExecutionPolicy, class LocalDynamicsFunction>
    void particle_for_split(const ExecutionPolicy &ex_policy, const LocalDynamicsFunction &local_dynamics_function);

    /** generalized particle search algorithm */
    template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
    void searchNeighborsByMesh(Mesh &mesh, DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
                               GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation);
    DiscreteVariable<UnsignedInt> *dvParticleIndex() { return dv_particle_index_; };
    DiscreteVariable<UnsignedInt> *dvCellOffset() { return dv_cell_offset_; };

    UnsignedInt TotalNumberOfCells() { return total_number_of_cells_; };
    template <typename DataType>
    DataType *initializeVariable(DiscreteVariable<DataType> *variable, DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType, typename... Args>
    DiscreteVariable<DataType> *registerDiscreteVariable(const std::string &name, size_t data_size, Args &&...args);

  protected:
    Kernel &kernel_;
    UnsignedInt total_number_of_cells_;
    UnsignedInt cell_offset_list_size_;
    UnsignedInt index_list_size_; // at least number_of_cells_pluse_one_
    DiscreteVariable<UnsignedInt> *dv_particle_index_;
    DiscreteVariable<UnsignedInt> *dv_cell_offset_;
    /** using concurrent vectors due to writing conflicts when building the list */
    ConcurrentIndexVector *cell_index_lists_;
    /** non-concurrent list data rewritten for building neighbor list */
    ListDataVector *cell_data_lists_;
    ParticleVariables all_discrete_variables_;

    void initialize(BaseParticles &base_particles);
    void clearCellLists();
    void UpdateCellListData(BaseParticles &base_particles);
    void tagBodyPartByCellByMesh(Mesh &mesh, ConcurrentCellLists &cell_lists,
                                 ConcurrentIndexVector &cell_indexes,
                                 std::function<bool(Vecd, Real)> &check_included);
    void writeMeshFieldToPltByMesh(Mesh &mesh, std::ofstream &output_file);
    void tagBoundingCellsByMesh(Mesh &mesh, StdVec<CellLists> &cell_data_lists,
                                const BoundingBox &bounding_bounds, int axis);
    void findNearestListDataEntryByMesh(Mesh &mesh, Real &min_distance_sqr, ListData &nearest_entry,
                                        const Vecd &position);
    /** split algorithm */;
    template <class ExecutionPolicy, class LocalDynamicsFunction>
    void particle_for_split_by_mesh(const ExecutionPolicy &ex_policy, Mesh &mesh,
                                    const LocalDynamicsFunction &local_dynamics_function);
};

class NeighborSearch : public Mesh
{
  public:
    template <class ExecutionPolicy>
    NeighborSearch(const ExecutionPolicy &ex_policy, CellLinkedList &cell_linked_list);

    template <typename FunctionOnEach>
    void forEachSearch(UnsignedInt source_index, const Vecd *source_pos,
                       const FunctionOnEach &function, int depth = 1) const;

  protected:
    UnsignedInt *particle_index_;
    UnsignedInt *cell_offset_;
};

/**
 * @class CellLinkedList
 * @brief Defining a mesh cell linked list for a body.
 * 		  The meshes for all bodies share the same global coordinates.
 */
class CellLinkedList : public BaseCellLinkedList
{
  protected:
    Mesh *mesh_;
    SingularVariable<Mesh> *sv_mesh_;

  public:
    CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
                   BaseParticles &base_particles, SPHAdaptation &sph_adaptation);
    ~CellLinkedList() {};
    Mesh &getMesh() { return *mesh_; };
    SingularVariable<Mesh> *svMesh() { return sv_mesh_; };
    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override;

    template <class ExecutionPolicy>
    NeighborSearch createNeighborSearch(const ExecutionPolicy &ex_policy);
    UnsignedInt getCellOffsetListSize() { return cell_offset_list_size_; };
};

/**
 * @class MultilevelCellLinkedList
 * @brief Defining a multilevel mesh cell linked list for a body
 * 		  for multi-resolution particle configuration.
 */
class MultilevelCellLinkedList : public BaseCellLinkedList
{
  protected:
    Real *h_ratio_; /**< Smoothing length for each level. */
    int *level_;    /**< Mesh level for each particle. */

    /** determine mesh level from particle cutoff radius */
    inline UnsignedInt getMeshLevel(Real particle_cutoff_radius);

  public:
    MultilevelCellLinkedList(BoundingBox tentative_bounds,
                             Real reference_grid_spacing, UnsignedInt total_levels,
                             BaseParticles &base_particles, SPHAdaptation &sph_adaptation);
    virtual ~MultilevelCellLinkedList() {};
    void insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position) override;
    void InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position) override;
};
} // namespace SPH
#endif // MESH_CELL_LINKED_LIST_H
