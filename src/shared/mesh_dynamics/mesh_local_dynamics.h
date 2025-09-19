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

#include "base_geometry.h"
#include "base_implementation.h"
#include "execution_policy.h"
#include "kernel_tabulated_ck.h"
#include "mesh_iterators.hpp"
#include "mesh_with_data_packages.hpp"
#include "neighbor_method.h"

namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

template <typename DataType>
using MeshVariableData = MeshWithGridDataPackagesType::MeshVariableData<DataType>;

template <typename DataType>
using MeshVariable = MeshWithGridDataPackagesType::MeshVariable<DataType>;

template <typename DataType>
using BKGMeshVariable = MeshWithGridDataPackagesType::BKGMeshVariable<DataType>;
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
    static constexpr int pkg_size = 4;
    static constexpr int pkg_size_minus1 = pkg_size - 1;
    Arrayi all_cells_;
    Real grid_spacing_;
    Real data_spacing_;
    UnsignedInt num_singular_pkgs_;

    UnsignedInt SortIndexFromCellIndex(const Arrayi &cell_index);
    Arrayi CellIndexFromSortIndex(const UnsignedInt &sort_index);
    static void registerComputingKernel(execution::Implementation<Base> *implementation) {};
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_), grid_spacing_(encloser.grid_spacing_),
              shape_(&encloser.shape_), base_dynamics(&encloser),
              cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)){};
        void update(const Arrayi &cell_index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        Real grid_spacing_;
        Shape *shape_;
        BaseMeshLocalDynamics *base_dynamics;
        int *cell_contain_id_;
    };

  private:
    Shape &shape_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : shape_(&encloser.shape_), grid_spacing_(encloser.grid_spacing_),
              far_field_distance_(encloser.far_field_distance_),
              data_mesh_(&encloser.data_mesh_), coarse_mesh_(&encloser.coarse_mesh_),
              base_dynamics_(&encloser), probe_coarse_phi_(ex_policy, coarse_mesh_),
              cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
              cell_pkg_index_coarse_(encloser.bmv_cell_pkg_index__coarse_.DelegatedData(ex_policy)),
              pkg_cell_info_coarse_(encloser.dv_pkg_cell_info_coarse_.DelegatedData(ex_policy)){};
        void update(const Arrayi &cell_index);

      protected:
        Shape *shape_;
        Real grid_spacing_;
        Real far_field_distance_;
        MeshWithGridDataPackagesType *data_mesh_;
        MeshWithGridDataPackagesType *coarse_mesh_;
        BaseMeshLocalDynamics *base_dynamics_;
        ProbeSignedDistance probe_coarse_phi_;
        int *cell_contain_id_;
        UnsignedInt *cell_pkg_index_coarse_;
        std::pair<Arrayi, int> *pkg_cell_info_coarse_;
    };

  private:
    MeshWithGridDataPackagesType &coarse_mesh_;
    Shape &shape_;
    Real far_field_distance_;
    BKGMeshVariable<int> &bmv_cell_contain_id_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index__coarse_;
    DiscreteVariable<std::pair<Arrayi, int>> &dv_pkg_cell_info_coarse_;
};

/**
 * @class InnerCellTagging
 * @brief Distinguish and categorize the inner data packages within the level set mesh.
 */
class InnerCellTagging : public BaseMeshLocalDynamics
{
  public:
    explicit InnerCellTagging(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~InnerCellTagging() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : all_cells_(encloser.all_cells_),
              data_mesh_(&encloser.data_mesh_),
              base_dynamics(&encloser){};
        void update(const Arrayi &cell_index);

      protected:
        Arrayi all_cells_;
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;

        bool isInnerPackage(const Arrayi &cell_index);
    };
};

/**
 * @class InitializeCellPackageInfo
 * @brief Store the 1-D array package index for each occupied cell on the mesh.
 */
class InitializeCellPackageInfo : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeCellPackageInfo(MeshWithGridDataPackagesType &data_mesh);
    virtual ~InitializeCellPackageInfo() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_), base_dynamics(&encloser),
              num_singular_pkgs_(encloser.num_singular_pkgs_),
              pkg_cell_info_(encloser.dv_pkg_cell_info_.DelegatedData(ex_policy)),
              cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &package_index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
        UnsignedInt num_singular_pkgs_;
        std::pair<Arrayi, int> *pkg_cell_info_;
        UnsignedInt *cell_pkg_index_;
    };

  protected:
    DiscreteVariable<std::pair<Arrayi, int>> &dv_pkg_cell_info_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pkg_cell_info_(encloser.dv_pkg_cell_info_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
              cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
              data_mesh_(&encloser.data_mesh_), base_dynamics(&encloser),
              num_singular_pkgs_(encloser.num_singular_pkgs_){};
        void update(const UnsignedInt &package_index);

      protected:
        std::pair<Arrayi, int> *pkg_cell_info_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *cell_pkg_index_;
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
        UnsignedInt num_singular_pkgs_;
    };

  protected:
    DiscreteVariable<std::pair<Arrayi, int>> &dv_pkg_cell_info_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : index_handler_(encloser.data_mesh_.index_handler_.DelegatedData(ex_policy)),
              pkg_cell_info_(encloser.dv_pkg_cell_info_.DelegatedData(ex_policy)),
              shape_(&encloser.shape_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index);

      protected:
        MeshWithGridDataPackagesType::IndexHandler *index_handler_;
        std::pair<Arrayi, int> *pkg_cell_info_;
        Shape *shape_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
    };

  private:
    Shape &shape_;
    Real far_field_distance;
    DiscreteVariable<std::pair<Arrayi, int>> &dv_pkg_cell_info_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_), base_dynamics(&encloser),
              num_singular_pkgs_(encloser.num_singular_pkgs_),
              cell_contain_id_(
                  encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
        UnsignedInt num_singular_pkgs_;
        int *cell_contain_id_;
        MeshVariableData<Real> *phi_;
    };

  protected:
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_), base_dynamics(&encloser),
              cell_contain_id_(encloser.bmv_cell_contain_id_.DelegatedData(ex_policy)),
              cell_package_index_(encloser.bmv_cell_package_index_.DelegatedData(ex_policy)),
              count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)){};
        void update(const Arrayi &cell_index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)){};
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
    using SmoothingKernel = typename NeighborMethod<SingleValued>::SmoothingKernel;

  public:
    explicit UpdateKernelIntegrals(
        MeshWithGridDataPackagesType &data_mesh, NeighborMethod<SingleValued> &neighbor_method);
    virtual ~UpdateKernelIntegrals() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
              kernel_weight_(encloser.mv_kernel_weight_.DelegatedData(ex_policy)),
              kernel_gradient_(encloser.mv_kernel_gradient_.DelegatedData(ex_policy)),
              kernel_second_gradient_(encloser.mv_kernel_second_gradient_.DelegatedData(ex_policy)),
              pkg_cell_info_(encloser.dv_pkg_cell_info_.DelegatedData(ex_policy)),
              kernel_(ex_policy, encloser.neighbor_method_),
              index_handler_(encloser.data_mesh_.index_handler_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
              cell_pkg_index_(encloser.bmv_cell_pkg_index_.DelegatedData(ex_policy)),
              probe_signed_distance_(ex_policy, &encloser.data_mesh_),
              cutoff_radius_(encloser.neighbor_method_.CutOffRadius()),
              depth_(static_cast<int>(std::ceil((cutoff_radius_ - Eps) / data_spacing_))){};
        void update(const UnsignedInt &package_index)
        {
            Arrayi cell_index = pkg_cell_info_[package_index].first;
            assignByGrid(
                kernel_weight_, cell_index, [&](const Arrayi &grid_index) -> Real
                { return computeKernelIntegral(package_index, grid_index); });
            assignByGrid(
                kernel_gradient_, cell_index, [&](const Arrayi &grid_index) -> Vecd
                { return computeKernelGradientIntegral(package_index, grid_index); });
            assignByGrid(
                kernel_second_gradient_, cell_index, [&](const Arrayi &grid_index) -> Matd
                { return computeKernelSecondGradientIntegral(package_index, grid_index); });
        }

      protected:
        Real data_spacing_;
        Real global_h_ratio_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        MeshVariableData<Real> *kernel_weight_;
        MeshVariableData<Vecd> *kernel_gradient_;
        MeshVariableData<Matd> *kernel_second_gradient_;
        std::pair<Arrayi, int> *pkg_cell_info_;

        SmoothingKernel kernel_;
        MeshWithGridDataPackagesType::IndexHandler *index_handler_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *cell_pkg_index_;
        ProbeSignedDistance probe_signed_distance_;

        Real cutoff_radius_, depth_;
        Real computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &grid_index);
        Vecd computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &grid_index);
        Matd computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &grid_index);
        template <typename DataType, typename FunctionByGrid>
        void assignByGrid(MeshVariableData<DataType> *mesh_variable, const Arrayi &cell_index,
                          const FunctionByGrid &function_by_grid);

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
    };

  private:
    NeighborMethod<SingleValued> &neighbor_method_;
    Real global_h_ratio_;
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    DiscreteVariable<std::pair<Arrayi, int>> &dv_pkg_cell_info_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    BKGMeshVariable<UnsignedInt> &bmv_cell_pkg_index_;
    MeshVariable<Real> &mv_kernel_weight_;
    MeshVariable<Vecd> &mv_kernel_gradient_;
    MeshVariable<Matd> &mv_kernel_second_gradient_;
    Real far_field_distance;

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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;

        Real upwindDifference(Real sign, Real df_p, Real df_n)
        {
            if (sign * df_p >= 0.0 && sign * df_n >= 0.0)
                return df_n;
            if (sign * df_p <= 0.0 && sign * df_n <= 0.0)
                return df_p;
            if (sign * df_p > 0.0 && sign * df_n < 0.0)
                return 0.0;

            Real df = df_p;
            if (sign * df_p < 0.0 && sign * df_n > 0.0)
            {
                Real ss = sign * (fabs(df_p) - fabs(df_n)) / (df_p - df_n);
                if (ss > 0.0)
                    df = df_n;
            }

            return df;
        }
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : threshold_(encloser.threshold_),
              perturbation_(threshold_ * encloser.perturbation_ratio_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index, Real dt = 0.0);

      protected:
        Real threshold_, perturbation_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    Real threshold_, perturbation_ratio_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : threshold_(encloser.threshold_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index, Real dt = 0.0);

      protected:
        Real threshold_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    Real threshold_;
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;

        Real upwindDifference(Real sign, Real df_p, Real df_n)
        {
            if (sign * df_p >= 0.0 && sign * df_n >= 0.0)
                return df_n;
            if (sign * df_p <= 0.0 && sign * df_n <= 0.0)
                return df_p;
            if (sign * df_p > 0.0 && sign * df_n < 0.0)
                return 0.0;

            Real df = df_p;
            if (sign * df_p < 0.0 && sign * df_n > 0.0)
            {
                Real ss = sign * (fabs(df_p) - fabs(df_n)) / (df_p - df_n);
                if (ss > 0.0)
                    df = df_n;
            }

            return df;
        }
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
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.mv_near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
              count_modified_(encloser.sv_count_modified_.DelegatedData(ex_policy)){};
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

class ProbeIsWithinMeshBound : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeIsWithinMeshBound(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~ProbeIsWithinMeshBound() {};

    bool update(const Vecd &position)
    {
        bool is_bounded = true;
        Arrayi cell_pos = data_mesh_.CellIndexFromPosition(position);
        for (int i = 0; i != position.size(); ++i)
        {
            if (cell_pos[i] < 2)
                is_bounded = false;
            if (cell_pos[i] > (all_cells_[i] - 2))
                is_bounded = false;
        }
        return is_bounded;
    }
};
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_H
