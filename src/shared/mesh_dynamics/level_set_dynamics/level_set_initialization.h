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
#include "level_set_probe.h"
#include "mesh_dynamics_algorithm.h"

namespace SPH
{
class InitialCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit InitialCellTagging(MeshWithGridDataPackagesType &data_mesh, Shape &shape);
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
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
    BKGMeshVariable<int> &bmv_cell_contain_id_;
};

class InitialCellTaggingFromCoarse : public BaseMeshLocalDynamics
{
  public:
    explicit InitialCellTaggingFromCoarse(
        MeshWithGridDataPackagesType &data_mesh,
        MeshWithGridDataPackagesType &coarse_mesh, Shape &shape);
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
        UnsignedInt *cell_pkg_index_;
        IndexHandler index_handler_;
        IndexHandler coarse_index_handler_;
        Real grid_spacing_;
        Real far_field_distance_;
        ProbeSignedDistance probe_coarse_phi_;
        int *cell_contain_id_;
        UnsignedInt *cell_pkg_index_coarse_;
        int *pkg_type_coarse_;
    };

  private:
    MeshWithGridDataPackagesType &coarse_mesh_;
    Shape &shape_;
    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
    BKGMeshVariable<int> &bmv_cell_contain_id_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_coarse_;
    MetaVariable<int> &dv_pkg_type_coarse_;
};

/**
 * @class InnerCellTagging
 * @brief Distinguish and categorize the inner data packages within the level set mesh.
 */
class InnerCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit InnerCellTagging(MeshWithGridDataPackagesType &data_mesh);
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
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
};

/**
 * @class InitializeCellNeighborhood
 * @brief Store the indices of neighboring cells in a 1-D array for each occupied cell.
 */
class InitializeCellNeighborhood : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeCellNeighborhood(MeshWithGridDataPackagesType &data_mesh);
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
        UnsignedInt num_singular_pkgs_;
    };

  protected:
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
};

/**
 * @class InitializeBasicPackageData
 * @brief Initialize the `phi` and `near_interface_id` mesh data for each occupied cell.
 */
class InitializeBasicPackageData : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeBasicPackageData(MeshWithGridDataPackagesType &data_mesh, Shape &shape);
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
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
    };

  private:
    Shape &shape_;
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    MeshVariable<int> &mv_near_interface_id_;
    void initializeSingularPackages(UnsignedInt package_index, Real far_field_level_set);
};

class NearInterfaceCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit NearInterfaceCellTagging(MeshWithGridDataPackagesType &data_mesh);
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
        MeshVariableData<Real> *phi_;
    };

  protected:
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    BKGMeshVariable<int> &bmv_cell_contain_id_;
    MeshVariable<Real> &mv_phi_;
};

class CellContainDiffusion : public BaseMeshLocalDynamics
{
  public:
    explicit CellContainDiffusion(MeshWithGridDataPackagesType &data_mesh,
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
        int *cell_contain_id_;
        UnsignedInt *cell_package_index_;
        UnsignedInt *count_modified_;
    };

  protected:
    BKGMeshVariable<int> &bmv_cell_contain_id_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_package_index_;
    SingularVariable<UnsignedInt> &sv_count_modified_;
};

class FinishDataPackages
{
  public:
    explicit FinishDataPackages(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : mesh_data_(mesh_data), shape_(shape) {};
    virtual ~FinishDataPackages() {};

    void exec()
    {
        initialize_basic_data_for_a_package.exec();

        near_interface_cell_tagging.exec();
        while (sv_count_modified_.getValue() > 0)
        {
            sv_count_modified_.setValue(0);
            cell_contain_diffusion.exec();
        }

        initialize_cell_neighborhood.exec();
    };

  private:
    MeshWithGridDataPackagesType &mesh_data_;
    Shape &shape_;
    SingularVariable<UnsignedInt> sv_count_modified_{"CountModifiedCell", 1};

    MeshInnerDynamics<execution::ParallelPolicy, InitializeCellNeighborhood> initialize_cell_neighborhood{mesh_data_};
    MeshInnerDynamics<execution::ParallelPolicy, InitializeBasicPackageData> initialize_basic_data_for_a_package{mesh_data_, shape_};
    MeshInnerDynamics<execution::ParallelPolicy, NearInterfaceCellTagging> near_interface_cell_tagging{mesh_data_};
    MeshAllDynamics<execution::ParallelPolicy, CellContainDiffusion> cell_contain_diffusion{mesh_data_, sv_count_modified_};
};
} // namespace SPH
#endif // LEVEL_SET_INITIALIZATION_H
