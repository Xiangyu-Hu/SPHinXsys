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
 * @file    level_set_initialization.h
 * @brief   This class encapsulates functions related to mesh dynamics,
 *          including initialization and manipulation of level set data.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef LEVEL_SET_INITIALIZATION_H
#define LEVEL_SET_INITIALIZATION_H

#include "base_geometry.h"
#include "base_local_mesh_dynamics.h"
#include "mesh_dynamics_algorithm.h"

namespace SPH
{
class InitialCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit InitialCellTagging(SparseMeshField<4> &data_mesh, Shape &shape);
    virtual ~InitialCellTagging() {};
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const Arrayi &cell_index);

      protected:
        ConcurrentVec<std::pair<UnsignedInt, int>> *occupied_data_pkgs_;
        UnsignedInt *cell_pkg_index_;
        IndexHandler index_handler_;
        Real grid_spacing_;
        Shape *shape_;
        int *cell_contain_id_;
    };

  private:
    Shape &shape_;
    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs_;
    CellVariable<UnsignedInt> &mcv_cell_pkg_index_;
    CellVariable<int> &mcv_cell_contain_id_;
};

class InitialCellTaggingFromCoarse : public BaseMeshLocalDynamics
{
  using ProbeCoarsePhi = SparseMeshField<4>::ProbeMesh<Real>;

  public:
    InitialCellTaggingFromCoarse(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
        SparseMeshField<4> &coarse_mesh, UnsignedInt coarse_resolution_level,
        Shape &shape);
    virtual ~InitialCellTaggingFromCoarse() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const Arrayi &cell_index);

      protected:
        Shape *shape_;
        ConcurrentVec<std::pair<UnsignedInt, int>> *occupied_data_pkgs_;
        UnsignedInt boundary_pkg_index_offset_;
        UnsignedInt *cell_pkg_index_;
        IndexHandler index_handler_;
        IndexHandler coarse_index_handler_;
        Real grid_spacing_;
        Real far_field_distance_;
        ProbeCoarsePhi probe_coarse_phi_;
        int *cell_contain_id_;
        UnsignedInt *cell_pkg_index_coarse_;
        int *pkg_type_coarse_;
    };

  private:
    SparseMeshField<4> &coarse_mesh_;
    UnsignedInt coarse_resolution_level_;
    Shape &shape_;
    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs_;
    UnsignedInt boundary_pkg_index_offset_;
    CellVariable<UnsignedInt> &mcv_cell_pkg_index_;
    CellVariable<int> &mcv_cell_contain_id_;
    CellVariable<UnsignedInt> &mcv_cell_pkg_index_coarse_;
    MetaVariable<int> &dv_pkg_type_coarse_;
};

/**
 * @class InnerCellTagging
 * @brief Distinguish and categorize the inner data packages within the level set mesh.
 */
class InnerCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit InnerCellTagging(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
    virtual ~InnerCellTagging() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const Arrayi &cell_index);

      protected:
        ConcurrentVec<std::pair<UnsignedInt, int>> *occupied_data_pkgs_;
        IndexHandler index_handler_;
        UnsignedInt *cell_pkg_index_;

        bool isNearInitiallyTagged(const Arrayi &cell_index);
        bool isInitiallyTagged(const Arrayi &cell_index);
    };

    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs_;
    CellVariable<UnsignedInt> &mcv_cell_pkg_index_;
};

/**
 * @class InitializeCellNeighborhood
 * @brief Store the indices of neighboring cells in a 1-D array for each occupied cell.
 */
class InitializeCellNeighborhood : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeCellNeighborhood(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
    virtual ~InitializeCellNeighborhood() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &package_index);

      protected:
        IndexHandler index_handler_;
        UnsignedInt *pkg_1d_cell_index_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *cell_pkg_index_;
        UnsignedInt num_boundary_pkgs_;
    };

  protected:
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    CellVariable<UnsignedInt> &mcv_cell_pkg_index_;
};

/**
 * @class InitializeBasicPackageData
 * @brief Initialize the `phi` and `near_interface_id` mesh data for each occupied cell.
 */
class InitializeBasicPackageData : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeBasicPackageData(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level, Shape &shape);
    virtual ~InitializeBasicPackageData() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        IndexHandler index_handler_;
        UnsignedInt *pkg_1d_cell_index_;
        Shape *shape_;
        PackageVariableData<Real> *phi_;
        PackageVariableData<int> *near_interface_id_;
    };

  private:
    Shape &shape_;
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    PackageVariable<Real> &mv_phi_;
    PackageVariable<Vecd> &mv_phi_gradient_;
    PackageVariable<int> &mv_near_interface_id_;
};

class NearInterfaceCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit NearInterfaceCellTagging(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
    virtual ~NearInterfaceCellTagging() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        UnsignedInt *pkg_1d_cell_index_;
        int *cell_contain_id_;
        PackageVariableData<Real> *phi_;
    };

  protected:
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    CellVariable<int> &mcv_cell_contain_id_;
    PackageVariable<Real> &mv_phi_;
};

class CellContainDiffusion : public BaseMeshLocalDynamics
{
  public:
    explicit CellContainDiffusion(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
        SingularVariable<UnsignedInt> &sv_count_modified);
    virtual ~CellContainDiffusion() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const Arrayi &cell_index);

      protected:
        IndexHandler index_handler_;
        UnsignedInt boundary_pkg_index_offset_;
        int *cell_contain_id_;
        UnsignedInt *cell_pkg_index_;
        UnsignedInt *count_modified_;
    };

  protected:
    UnsignedInt boundary_pkg_index_offset_;
    CellVariable<int> &mcv_cell_contain_id_;
    CellVariable<UnsignedInt> &mcv_cell_package_index_;
    SingularVariable<UnsignedInt> &sv_count_modified_;
};

class FinishPackageDatas : public BaseDynamics<void>
{
  public:
    FinishPackageDatas(SparseMeshField<4> &mesh_data,
                       UnsignedInt resolution_level, Shape &shape);
    virtual ~FinishPackageDatas() {};
    void exec(Real dt = 0.0) override;

  private:
    SparseMeshField<4> &mesh_data_;
    UnsignedInt resolution_level_;
    Shape &shape_;
    SingularVariable<UnsignedInt> sv_count_modified_{"CountModifiedCell", 1};

    MeshInnerDynamics<execution::ParallelPolicy, InitializeCellNeighborhood> initialize_cell_neighborhood;
    MeshInnerDynamics<execution::ParallelPolicy, InitializeBasicPackageData> initialize_basic_data_for_a_package;
    MeshInnerDynamics<execution::ParallelPolicy, NearInterfaceCellTagging> near_interface_cell_tagging;
    MeshAllDynamics<execution::ParallelPolicy, CellContainDiffusion> cell_contain_diffusion;
};
} // namespace SPH
#endif // LEVEL_SET_INITIALIZATION_H
