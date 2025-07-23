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
 * @file    mesh_local_dynamics.h
 * @brief   This class encapsulates functions related to mesh dynamics,
 *          including initialization and manipulation of level set data.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_LOCAL_DYNAMICS_H
#define MESH_LOCAL_DYNAMICS_H

#include "base_geometry.h"
#include "base_implementation.h"
#include "base_kernel.h"
#include "data_type.h"
#include "execution_policy.h"
#include "mesh_iterators.hpp"
#include "mesh_with_data_packages.h"
#include "sphinxsys_variable.h"
namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

template <typename DataType>
using MeshVariableData = MeshWithGridDataPackagesType::MeshVariableData<DataType>;

template <typename DataType>
using MeshVariable = MeshWithGridDataPackagesType::MeshVariable<DataType>;
/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackagesType &data_mesh)
        : data_mesh_(data_mesh),
          all_cells_(data_mesh.AllCells()),
          grid_spacing_(data_mesh.GridSpacing()),
          data_spacing_(data_mesh.DataSpacing()),
          meta_data_cell_(data_mesh.meta_data_cell_),
          cell_neighborhood_(data_mesh.cell_neighborhood_),
          cell_package_index_(data_mesh.cell_package_index_),
          phi_(*data_mesh.getMeshVariable<Real>("Levelset")),
          phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelsetGradient")),
          near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
          kernel_weight_(*data_mesh.getMeshVariable<Real>("KernelWeight")),
          kernel_gradient_(*data_mesh.getMeshVariable<Vecd>("KernelGradient")),
          kernel_second_gradient_(*data_mesh.getMeshVariable<Matd>("KernelSecondGradient")) {};
    virtual ~BaseMeshLocalDynamics() {};

    MeshWithGridDataPackagesType &data_mesh_;
    static constexpr int pkg_size = 4;
    Arrayi all_cells_;
    Real grid_spacing_;
    Real data_spacing_;
    DiscreteVariable<std::pair<Arrayi, int>> &meta_data_cell_;
    DiscreteVariable<CellNeighborhood> &cell_neighborhood_;
    DiscreteVariable<size_t> &cell_package_index_;

    MeshVariable<Real> &phi_;
    MeshVariable<Vecd> &phi_gradient_;
    MeshVariable<int> &near_interface_id_;
    MeshVariable<Real> &kernel_weight_;
    MeshVariable<Vecd> &kernel_gradient_;
    MeshVariable<Matd> &kernel_second_gradient_;

    size_t SortIndexFromCellIndex(const Arrayi &cell_index);
    Arrayi CellIndexFromSortIndex(const size_t &sort_index);

    static void registerComputingKernel(execution::Implementation<Base> *implementation) {};
};

class ProbeSignedDistance : public ProbeMesh<Real, 4>
{
  public:
    template <class ExecutionPolicy>
    ProbeSignedDistance(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "Levelset"){};
    ~ProbeSignedDistance() {};
};

class ProbeLevelSetGradient : public ProbeMesh<Vecd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeLevelSetGradient(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *data_mesh)
        : ProbeMesh(ex_policy, data_mesh, "LevelsetGradient"){};
    ~ProbeLevelSetGradient() {};
};

class ProbeNormalDirection : public ProbeMesh<Vecd, 4>
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeNormalDirection(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh<Vecd, 4>(ex_policy, mesh_data, "LevelsetGradient"){};
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
 * @class InitializeDataForSingularPackage
 * @brief Update function for singular data initialization.
 */
class InitializeDataForSingularPackage : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataForSingularPackage(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~InitializeDataForSingularPackage() {};

    void update(const size_t package_index, Real far_field_level_set);
};

/**
 * @class InitializeDataInACell
 * @brief Distinguish and categorize the core data packages within the level set mesh.
 */
class InitializeDataInACell : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataInACell(MeshWithGridDataPackagesType &data_mesh, Shape &shape)
        : BaseMeshLocalDynamics(data_mesh),
          shape_(shape) {};
    virtual ~InitializeDataInACell() {};
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_),
              grid_spacing_(encloser.grid_spacing_),
              shape_(&encloser.shape_),
              base_dynamics(&encloser){};
        void update(const Arrayi &cell_index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        Real grid_spacing_;
        Shape *shape_;
        BaseMeshLocalDynamics *base_dynamics;
    };

  private:
    Shape &shape_;
};

/**
 * @class TagACellIsInnerPackage
 * @brief Distinguish and categorize the inner data packages within the level set mesh.
 */
class TagACellIsInnerPackage : public BaseMeshLocalDynamics
{
  public:
    explicit TagACellIsInnerPackage(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~TagACellIsInnerPackage() {};

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
 * @class InitializeIndexMesh
 * @brief Store the 1-D array package index for each occupied cell on the mesh.
 */
class InitializeIndexMesh : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeIndexMesh(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~InitializeIndexMesh() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_mesh_(&encloser.data_mesh_),
              base_dynamics(&encloser){};
        void update(const size_t &package_index);

      protected:
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
    };
};

/**
 * @class InitializeCellNeighborhood
 * @brief Store the indices of neighboring cells in a 1-D array for each occupied cell.
 */
class InitializeCellNeighborhood : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeCellNeighborhood(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~InitializeCellNeighborhood() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : meta_data_cell_(encloser.meta_data_cell_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)),
              cell_package_index_(encloser.cell_package_index_.DelegatedData(ex_policy)),
              data_mesh_(&encloser.data_mesh_),
              base_dynamics(&encloser){};
        void update(const size_t &package_index);

      protected:
        std::pair<Arrayi, int> *meta_data_cell_;
        CellNeighborhood *cell_neighborhood_;
        size_t *cell_package_index_;
        MeshWithGridDataPackagesType *data_mesh_;
        BaseMeshLocalDynamics *base_dynamics;
    };
};

/**
 * @class InitializeBasicDataForAPackage
 * @brief Initialize the `phi` and `near_interface_id` mesh data for each occupied cell.
 */
class InitializeBasicDataForAPackage : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeBasicDataForAPackage(MeshWithGridDataPackagesType &data_mesh, Shape &shape)
        : BaseMeshLocalDynamics(data_mesh),
          shape_(shape) {};
    virtual ~InitializeBasicDataForAPackage() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : index_handler_(encloser.data_mesh_.index_handler_.DelegatedData(ex_policy)),
              meta_data_cell_(encloser.meta_data_cell_.DelegatedData(ex_policy)),
              shape_(&encloser.shape_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.near_interface_id_.DelegatedData(ex_policy)){};
        void update(const size_t &index);

      protected:
        MeshWithGridDataPackagesType::IndexHandler *index_handler_;
        std::pair<Arrayi, int> *meta_data_cell_;
        Shape *shape_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
    };

  private:
    Shape &shape_;
};

/**
 * @class UpdateLevelSetGradient
 * @brief Compute `phi_gradient_` mesh data base on `phi_` state for each occupied cell.
 */
class UpdateLevelSetGradient : public BaseMeshLocalDynamics
{
  public:
    explicit UpdateLevelSetGradient(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~UpdateLevelSetGradient() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.phi_gradient_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const size_t &index);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        CellNeighborhood *cell_neighborhood_;
    };
};

template <class KernelType>
class UpdateKernelIntegrals : public BaseMeshLocalDynamics
{
  public:
    explicit UpdateKernelIntegrals(MeshWithGridDataPackagesType &data_mesh, KernelType *kernel, Real global_h_ratio)
        : BaseMeshLocalDynamics(data_mesh),
          kernel_(kernel),
          global_h_ratio_(global_h_ratio) {};
    virtual ~UpdateKernelIntegrals() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              global_h_ratio_(encloser.global_h_ratio_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.phi_gradient_.DelegatedData(ex_policy)),
              kernel_weight_(encloser.kernel_weight_.DelegatedData(ex_policy)),
              kernel_gradient_(encloser.kernel_gradient_.DelegatedData(ex_policy)),
              kernel_second_gradient_(encloser.kernel_second_gradient_.DelegatedData(ex_policy)),
              meta_data_cell_(encloser.meta_data_cell_.DelegatedData(ex_policy)),
              kernel_(encloser.kernel_),
              index_handler_(encloser.data_mesh_.index_handler_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)),
              cell_package_index_(encloser.cell_package_index_.DelegatedData(ex_policy)),
              probe_signed_distance_(ex_policy, &encloser.data_mesh_){};
        void update(const size_t &package_index)
        {
            Arrayi cell_index = meta_data_cell_[package_index].first;
            assignByPosition(
                kernel_weight_, cell_index, [&](const Vecd &position, const Arrayi &grid_index) -> Real
                { return computeKernelIntegral(position, package_index, grid_index); });
            assignByPosition(
                kernel_gradient_, cell_index, [&](const Vecd &position, const Arrayi &grid_index) -> Vecd
                { return computeKernelGradientIntegral(position, package_index, grid_index); });
            assignByPosition(
                kernel_second_gradient_, cell_index, [&](const Vecd &position, const Arrayi &grid_index) -> Matd
                { return computeKernelSecondGradientIntegral(position, package_index, grid_index); });
        }

      protected:
        Real data_spacing_;
        Real global_h_ratio_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<Vecd> *phi_gradient_;
        MeshVariableData<Real> *kernel_weight_;
        MeshVariableData<Vecd> *kernel_gradient_;
        MeshVariableData<Matd> *kernel_second_gradient_;
        std::pair<Arrayi, int> *meta_data_cell_;

        KernelType *kernel_;
        MeshWithGridDataPackagesType::IndexHandler *index_handler_;
        CellNeighborhood *cell_neighborhood_;
        size_t *cell_package_index_;
        ProbeSignedDistance probe_signed_distance_;

        Real cutoff_radius_;
        Real threshold;
        Real computeKernelIntegral(const Vecd &position, const size_t &package_index, const Arrayi &grid_index);
        Vecd computeKernelGradientIntegral(const Vecd &position, const size_t &package_index, const Arrayi &grid_index);
        Matd computeKernelSecondGradientIntegral(const Vecd &position, const size_t &package_index, const Arrayi &grid_index);
        template <typename DataType, typename FunctionByPosition>
        void assignByPosition(MeshVariableData<DataType> *mesh_variable, const Arrayi &cell_index,
                              const FunctionByPosition &function_by_position);

        /** a cut cell is a cut by the level set. */
        /** "Multi-scale modeling of compressible multi-fluid flows with conservative interface method."
         * Hu, X. Y., et al., Proceedings of the Summer Program. Vol. 301. Stanford, CA, USA:
         * Center for Turbulence Research, Stanford University, 2010.*/
        Real CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing)
        {
            Real squared_norm_inv = 1.0 / (phi_gradient.squaredNorm() + TinyReal);
            Real volume_fraction(0);
            for (size_t i = 0; i != Dimensions; ++i)
            {
                volume_fraction += phi_gradient[i] * phi_gradient[i] * squared_norm_inv *
                                   Heaviside(phi / (ABS(phi_gradient[i]) + TinyReal), 0.5 * data_spacing);
            }
            return volume_fraction;
        }
    };

  private:
    KernelType *kernel_;
    Real global_h_ratio_;
};

class ReinitializeLevelSet : public BaseMeshLocalDynamics
{
  public:
    explicit ReinitializeLevelSet(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~ReinitializeLevelSet() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const size_t &index);

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
};

class MarkNearInterface : public BaseMeshLocalDynamics
{
  public:
    explicit MarkNearInterface(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~MarkNearInterface() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const size_t &index, Real small_shift_factor);

      protected:
        Real data_spacing_;
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };
};

class RedistanceInterface : public BaseMeshLocalDynamics
{
  public:
    explicit RedistanceInterface(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~RedistanceInterface() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : data_spacing_(encloser.data_spacing_),
              phi_(encloser.phi_.DelegatedData(ex_policy)),
              phi_gradient_(encloser.phi_gradient_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const size_t &index);

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
};

class DiffuseLevelSetSign : public BaseMeshLocalDynamics
{
  public:
    explicit DiffuseLevelSetSign(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~DiffuseLevelSetSign() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : phi_(encloser.phi_.DelegatedData(ex_policy)),
              near_interface_id_(encloser.near_interface_id_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)){};
        void update(const size_t &index);

      protected:
        MeshVariableData<Real> *phi_;
        MeshVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };
};

class InitializeDataInACellFromCoarse : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataInACellFromCoarse(MeshWithGridDataPackagesType &data_mesh, MeshWithGridDataPackagesType &coarse_mesh, Shape &shape)
        : BaseMeshLocalDynamics(data_mesh),
          coarse_mesh_(coarse_mesh),
          shape_(shape) {};
    virtual ~InitializeDataInACellFromCoarse() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : shape_(&encloser.shape_),
              grid_spacing_(encloser.grid_spacing_),
              data_mesh_(&encloser.data_mesh_),
              coarse_mesh_(&encloser.coarse_mesh_),
              base_dynamics_(&encloser),
              probe_coarse_phi_(ex_policy, coarse_mesh_){};
        void update(const Arrayi &cell_index);

      protected:
        Shape *shape_;
        Real grid_spacing_;
        MeshWithGridDataPackagesType *data_mesh_;
        MeshWithGridDataPackagesType *coarse_mesh_;
        BaseMeshLocalDynamics *base_dynamics_;
        ProbeSignedDistance probe_coarse_phi_;
    };

  private:
    MeshWithGridDataPackagesType &coarse_mesh_;
    Shape &shape_;
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

class WriteMeshFieldToPlt : public BaseMeshLocalDynamics
{
  public:
    explicit WriteMeshFieldToPlt(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh) {};
    virtual ~WriteMeshFieldToPlt() {};

    void update(std::ofstream &output_file);
};
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_H
