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
 * @file    mesh_local_dynamics.h
 * @brief   This class encapsulates functions related to mesh dynamics,
 *          including initialization and manipulation of level set data.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_LOCAL_DYNAMICS_H
#define MESH_LOCAL_DYNAMICS_H

#include "base_dynamics.h"
#include "base_geometry.h"
#include "base_implementation.h"
#include "execution_policy.h"
#include "kernel_tabulated_ck.h"
#include "mesh_iterators.hpp"
#include "mesh_with_data_packages.hpp"
#include "neighbor_method.h"
#include "data_package_function.hpp"

namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

template <typename DataType>
using MeshVariableData = MeshWithGridDataPackagesType::MeshVariableData<DataType>;

template <typename DataType>
using MeshVariable = MeshWithGridDataPackagesType::MeshVariable<DataType>;

template <typename DataType>
using BKGMeshVariable = MeshWithGridDataPackagesType::BKGMeshVariable<DataType>;

template <typename DataType>
using MetaVariable = MeshWithGridDataPackagesType::MetaVariable<DataType>;

using MeshVariableAssemble = MeshWithGridDataPackagesType::MeshVariableAssemble;
using BKGMeshVariableAssemble = MeshWithGridDataPackagesType::BKGMeshVariableAssemble;
using MetaVariableAssemble = MeshWithGridDataPackagesType::MetaVariableAssemble;
using IndexHandler = MeshWithGridDataPackagesType::IndexHandler;

/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackagesType &data_mesh);
    virtual ~BaseMeshLocalDynamics() {};

    MeshWithGridDataPackagesType &data_mesh_;
    IndexHandler &index_handler_;
    static constexpr int pkg_size = MeshWithGridDataPackagesType::DataPackageSize();
    static constexpr int pkg_size_minus1 = pkg_size - 1;
};

class ProbeSignedDistance : public ProbeMesh<Real, 4>
{
  public:
    template <class ExecutionPolicy>
    ProbeSignedDistance(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "LevelSet"){};
    ~ProbeSignedDistance() {};
};

class ProbeLevelSetGradient : public ProbeMesh<Vecd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeLevelSetGradient(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "LevelSetGradient"){};
    ~ProbeLevelSetGradient() {};
};

class ProbeNormalDirection : public ProbeMesh<Vecd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeNormalDirection(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh<Vecd, 4>(ex_policy, mesh_data, "LevelSetGradient"){};
    virtual ~ProbeNormalDirection() {};

    Vecd operator()(const Vecd &position)
    {
        return ProbeMesh<Vecd, 4>::operator()(position).normalized();
    }
};

class ProbeKernelIntegral : public ProbeMesh<Real, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "KernelWeight"){};
    ~ProbeKernelIntegral() {};
};

class ProbeKernelGradientIntegral : public ProbeMesh<Vecd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelGradientIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "KernelGradient"){};
    ~ProbeKernelGradientIntegral() {};
};

class ProbeKernelSecondGradientIntegral : public ProbeMesh<Matd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelSecondGradientIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "KernelSecondGradient"){};
    ~ProbeKernelSecondGradientIntegral() {};
};

/**
 * @class InitialCellTagging
 * @brief Distinguish and categorize the core data packages within the level set mesh.
 */
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
        bool isInitiallyTagged(const Arrayi &cell_index)
        {
            UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
            return cell_pkg_index_[index_1d] == 2;
        };
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

/**
 * @class UpdateLevelSetGradient
 * @brief Compute `phi_gradient_` mesh data base on `phi_` state for each occupied cell.
 */
class UpdateLevelSetGradient : public BaseMeshLocalDynamics
{
  public:
    explicit UpdateLevelSetGradient(MeshWithGridDataPackagesType &data_mesh);
    virtual ~UpdateLevelSetGradient() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class UpdateKernelIntegrals : public BaseMeshLocalDynamics
{
    using SmoothingKernel =
        typename NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel;

  public:
    explicit UpdateKernelIntegrals(
        MeshWithGridDataPackagesType &data_mesh,
        NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method);
    virtual ~UpdateKernelIntegrals() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &package_index);

      protected:
        Real global_h_ratio_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        MeshVariableData<Real> *kernel_weight_;
        MeshVariableData<Vecd> *kernel_gradient_;
        MeshVariableData<Matd> *kernel_second_gradient_;
        UnsignedInt *pkg_1d_cell_index_;

        SmoothingKernel kernel_;
        IndexHandler index_handler_;
        Real data_spacing_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *cell_pkg_index_;
        ProbeSignedDistance probe_signed_distance_;

        Real cutoff_radius_, depth_;
        BoundingBoxi bounding_box_;
        Real computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &data_index);
        Vecd computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index);
        Matd computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index);

        /** a cut cell is a cut by the level set. */
        /** "Multi-scale modeling of compressible multi-fluid flows with conservative interface method."
         * Hu, X. Y., et al., Proceedings of the Summer Program. Vol. 301. Stanford, CA, USA:
         * Center for Turbulence Research, Stanford University, 2010.*/
        Real CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing)
        {
            Real squared_norm_inv = 1.0 / (phi_gradient.squaredNorm() + TinyReal);
            Real volume_fraction(0);
            for (UnsignedInt i = 0; i != Dimensions; ++i)
            {
                volume_fraction += phi_gradient[i] * phi_gradient[i] * squared_norm_inv *
                                   Heaviside(phi / (ABS(phi_gradient[i]) + TinyReal), 0.5 * data_spacing);
            }
            return volume_fraction;
        }

        template <typename DataType, typename FunctionByDataIndex>
        DataType computeIntegral(Real phi, const UnsignedInt &package_index, const Arrayi &grid_index,
                                 const DataType &initial_value, const FunctionByDataIndex &function_by_index);
    };

  private:
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method_;
    Real global_h_ratio_;
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    MetaVariable<UnsignedInt> &dv_pkg_1d_cell_index_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
    MeshVariable<Real> &mv_kernel_weight_;
    MeshVariable<Vecd> &mv_kernel_gradient_;
    MeshVariable<Matd> &mv_kernel_second_gradient_;

    void initializeSingularPackages(UnsignedInt package_index, Real far_field_level_set);
};

class ReinitializeLevelSet : public BaseMeshLocalDynamics
{
  public:
    explicit ReinitializeLevelSet(MeshWithGridDataPackagesType &data_mesh);
    virtual ~ReinitializeLevelSet() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;

        Real upwindDifference(Real sign, Real df_p, Real df_n);
    };

  protected:
    MeshVariable<Real> &mv_phi_;
    MeshVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class MarkCutInterfaces : public BaseMeshLocalDynamics
{
  public:
    explicit MarkCutInterfaces(MeshWithGridDataPackagesType &data_mesh, Real perturbation_ratio);
    virtual ~MarkCutInterfaces() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index, Real dt = 0.0);

      protected:
        IndexHandler index_handler_;
        Real threshold_, perturbation_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    Real perturbation_ratio_;
    MeshVariable<Real> &mv_phi_;
    MeshVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class MarkNearInterface : public BaseMeshLocalDynamics
{
  public:
    explicit MarkNearInterface(MeshWithGridDataPackagesType &data_mesh);
    virtual ~MarkNearInterface() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index, Real dt = 0.0);

      protected:
        IndexHandler index_handler_;
        Real threshold_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    MeshVariable<Real> &mv_phi_;
    MeshVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class RedistanceInterface : public BaseMeshLocalDynamics
{
  public:
    explicit RedistanceInterface(MeshWithGridDataPackagesType &data_mesh);
    virtual ~RedistanceInterface() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    MeshVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class DiffuseLevelSetSign : public BaseMeshLocalDynamics
{
  public:
    explicit DiffuseLevelSetSign(MeshWithGridDataPackagesType &data_mesh,
                                 SingularVariable<UnsignedInt> &sv_count_modified);
    virtual ~DiffuseLevelSetSign() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *count_modified_;
    };

  protected:
    MeshVariable<Real> &mv_phi_;
    MeshVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    SingularVariable<UnsignedInt> &sv_count_modified_;
};
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_H
