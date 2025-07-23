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
/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackagesType &mesh_data)
        : mesh_data_(mesh_data),
          all_cells_(mesh_data.AllCells()),
          grid_spacing_(mesh_data.GridSpacing()),
          data_spacing_(mesh_data.DataSpacing()),
          meta_data_cell_(mesh_data.meta_data_cell_),
          cell_neighborhood_(mesh_data.cell_neighborhood_),
          cell_package_index_(mesh_data.cell_package_index_),
          phi_(*mesh_data.getMeshVariable<Real>("Levelset")),
          phi_gradient_(*mesh_data.getMeshVariable<Vecd>("LevelsetGradient")),
          near_interface_id_(*mesh_data.getMeshVariable<int>("NearInterfaceID")),
          kernel_weight_(*mesh_data.getMeshVariable<Real>("KernelWeight")),
          kernel_gradient_(*mesh_data.getMeshVariable<Vecd>("KernelGradient")),
          kernel_second_gradient_(*mesh_data.getMeshVariable<Matd>("KernelSecondGradient")) {};
    virtual ~BaseMeshLocalDynamics() {};

    MeshWithGridDataPackagesType &mesh_data_;
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

class ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeMesh(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : index_handler_(mesh_data->index_handler_.DelegatedData(ex_policy)),
          cell_package_index_(mesh_data->cell_package_index_.DelegatedData(ex_policy)),
          cell_neighborhood_(mesh_data->cell_neighborhood_.DelegatedData(ex_policy)){};

    /** This function probe a mesh value */
    template <class DataType>
    DataType probeMesh(MeshVariableData<DataType> *mesh_variable_data, const Vecd &position);

  protected:
    static constexpr int pkg_size = 4;
    MeshWithGridDataPackagesType::IndexHandler *index_handler_;
    size_t *cell_package_index_;
    CellNeighborhood *cell_neighborhood_;

    /** probe by applying bi and tri-linear interpolation within the package. */
    template <class DataType>
    DataType probeDataPackage(MeshVariableData<DataType> *mesh_variable_data,
                              size_t package_index,
                              const Arrayi &cell_index, const Vecd &position);
};

class ProbeSignedDistance : public ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeSignedDistance(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh(ex_policy, mesh_data),
          phi_(mesh_data->getMeshVariable<Real>("Levelset")->DelegatedData(ex_policy)){};
    ~ProbeSignedDistance() {};

    Real update(const Vecd &position) { return probeMesh(phi_, position); };

  private:
    MeshVariableData<Real> *phi_;
};

class ProbeLevelSetGradient : public ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeLevelSetGradient(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh(ex_policy, mesh_data),
          phi_gradient_(mesh_data->getMeshVariable<Vecd>("LevelsetGradient")->DelegatedData(ex_policy)){};
    ~ProbeLevelSetGradient() {};

    Vecd update(const Vecd &position) { return probeMesh(phi_gradient_, position); };

  private:
    MeshVariableData<Vecd> *phi_gradient_;
};

class ProbeKernelIntegral : public ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh(ex_policy, mesh_data),
          kernel_weight_(mesh_data->getMeshVariable<Real>("KernelWeight")->DelegatedData(ex_policy)){};
    ~ProbeKernelIntegral() {};

    Real update(const Vecd &position) { return probeMesh(kernel_weight_, position); };

  private:
    MeshVariableData<Real> *kernel_weight_;
};

class ProbeKernelGradientIntegral : public ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelGradientIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh(ex_policy, mesh_data),
          kernel_gradient_(mesh_data->getMeshVariable<Vecd>("KernelGradient")->DelegatedData(ex_policy)){};
    ~ProbeKernelGradientIntegral() {};

    Vecd update(const Vecd &position) { return probeMesh(kernel_gradient_, position); };

  private:
    MeshVariableData<Vecd> *kernel_gradient_;
};

class ProbeKernelSecondGradientIntegral : public ProbeMesh
{
  public:
    template <class ExecutionPolicy>
    explicit ProbeKernelSecondGradientIntegral(const ExecutionPolicy &ex_policy, MeshWithGridDataPackagesType *mesh_data)
        : ProbeMesh(ex_policy, mesh_data),
          kernel_second_gradient_(mesh_data->getMeshVariable<Matd>("KernelSecondGradient")->DelegatedData(ex_policy)){};
    ~ProbeKernelSecondGradientIntegral() {};

    Matd update(const Vecd &position) { return probeMesh(kernel_second_gradient_, position); };

  private:
    MeshVariableData<Matd> *kernel_second_gradient_;
};

/**
 * @class InitializeDataForSingularPackage
 * @brief Update function for singular data initialization.
 */
class InitializeDataForSingularPackage : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataForSingularPackage(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit InitializeDataInACell(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          shape_(shape) {};
    virtual ~InitializeDataInACell() {};
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : mesh_data_(&encloser.mesh_data_),
              grid_spacing_(encloser.grid_spacing_),
              shape_(&encloser.shape_),
              base_dynamics(&encloser){};
        void update(const Arrayi &cell_index);

      protected:
        MeshWithGridDataPackagesType *mesh_data_;
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
    explicit TagACellIsInnerPackage(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~TagACellIsInnerPackage() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : all_cells_(encloser.all_cells_),
              mesh_data_(&encloser.mesh_data_),
              base_dynamics(&encloser){};
        void update(const Arrayi &cell_index);

      protected:
        Arrayi all_cells_;
        MeshWithGridDataPackagesType *mesh_data_;
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
    explicit InitializeIndexMesh(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~InitializeIndexMesh() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : mesh_data_(&encloser.mesh_data_),
              base_dynamics(&encloser){};
        void update(const size_t &package_index);

      protected:
        MeshWithGridDataPackagesType *mesh_data_;
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
    explicit InitializeCellNeighborhood(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~InitializeCellNeighborhood() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : meta_data_cell_(encloser.meta_data_cell_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)),
              cell_package_index_(encloser.cell_package_index_.DelegatedData(ex_policy)),
              mesh_data_(&encloser.mesh_data_),
              base_dynamics(&encloser){};
        void update(const size_t &package_index);

      protected:
        std::pair<Arrayi, int> *meta_data_cell_;
        CellNeighborhood *cell_neighborhood_;
        size_t *cell_package_index_;
        MeshWithGridDataPackagesType *mesh_data_;
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
    explicit InitializeBasicDataForAPackage(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          shape_(shape) {};
    virtual ~InitializeBasicDataForAPackage() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : index_handler_(encloser.mesh_data_.index_handler_.DelegatedData(ex_policy)),
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
    explicit UpdateLevelSetGradient(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit UpdateKernelIntegrals(MeshWithGridDataPackagesType &mesh_data, KernelType *kernel, Real global_h_ratio)
        : BaseMeshLocalDynamics(mesh_data),
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
              index_handler_(encloser.mesh_data_.index_handler_.DelegatedData(ex_policy)),
              cell_neighborhood_(encloser.cell_neighborhood_.DelegatedData(ex_policy)),
              cell_package_index_(encloser.cell_package_index_.DelegatedData(ex_policy)),
              probe_signed_distance_(ex_policy, &encloser.mesh_data_){};
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
    explicit ReinitializeLevelSet(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit MarkNearInterface(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit RedistanceInterface(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit DiffuseLevelSetSign(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
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
    explicit InitializeDataInACellFromCoarse(MeshWithGridDataPackagesType &mesh_data, MeshWithGridDataPackagesType &coarse_mesh, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
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
              mesh_data_(&encloser.mesh_data_),
              coarse_mesh_(&encloser.coarse_mesh_),
              base_dynamics_(&encloser),
              probe_coarse_phi_(ex_policy, coarse_mesh_){};
        void update(const Arrayi &cell_index);

      protected:
        Shape *shape_;
        Real grid_spacing_;
        MeshWithGridDataPackagesType *mesh_data_;
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
    explicit ProbeIsWithinMeshBound(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeIsWithinMeshBound() {};

    bool update(const Vecd &position)
    {
        bool is_bounded = true;
        Arrayi cell_pos = mesh_data_.CellIndexFromPosition(position);
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
    explicit WriteMeshFieldToPlt(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~WriteMeshFieldToPlt() {};

    void update(std::ofstream &output_file);
};
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_H
