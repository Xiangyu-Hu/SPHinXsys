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
 * @file 	base_relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_LOCAL_DYNAMICS_H
#define MESH_LOCAL_DYNAMICS_H

#include "base_geometry.h"
#include "base_kernel.h"
#include "data_type.h"
#include "mesh_with_data_packages.hpp"
#include "sphinxsys_variable.h"

namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackages<4> &mesh_data)
        : mesh_data_(mesh_data),
          all_cells_(mesh_data.AllCells()),
          grid_spacing_(mesh_data.GridSpacing()),
          data_spacing_(mesh_data.DataSpacing()),
          phi_(*mesh_data.getMeshVariable<Real>("Levelset")),
          phi_gradient_(*mesh_data.getMeshVariable<Vecd>("LevelsetGradient")),
          near_interface_id_(*mesh_data.getMeshVariable<int>("NearInterfaceID")),
          kernel_weight_(*mesh_data.getMeshVariable<Real>("KernelWeight")),
          kernel_gradient_(*mesh_data.getMeshVariable<Vecd>("KernelGradient")),
          kernel_second_gradient_(*mesh_data.getMeshVariable<Matd>("KernelSecondGradient")) {};
    virtual ~BaseMeshLocalDynamics() {};

  protected:
    MeshWithGridDataPackages<4> &mesh_data_;
    Arrayi all_cells_;
    Real grid_spacing_;
    Real data_spacing_;

    MeshVariable<Real> &phi_;
    MeshVariable<Vecd> &phi_gradient_;
    MeshVariable<int> &near_interface_id_;
    MeshVariable<Real> &kernel_weight_;
    MeshVariable<Vecd> &kernel_gradient_;
    MeshVariable<Matd> &kernel_second_gradient_;

    size_t SortIndexFromCellIndex(const Arrayi &cell_index);
    Arrayi CellIndexFromSortIndex(const size_t &sort_index);
};

class InitializeDataForSingularPackage : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataForSingularPackage(MeshWithGridDataPackages<4> &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~InitializeDataForSingularPackage() {};

    void update(const size_t package_index, Real far_field_level_set);
};

class InitializeDataInACell : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataInACell(MeshWithGridDataPackages<4> &mesh_data, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          shape_(shape) {};
    virtual ~InitializeDataInACell() {};

    void update(const Arrayi &index);

  private:
    Shape &shape_;
};

class TagACellIsInnerPackage : public BaseMeshLocalDynamics
{
  public:
    explicit TagACellIsInnerPackage(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~TagACellIsInnerPackage() {};

    void update(const Arrayi &index);

  private:
    bool isInnerPackage(const Arrayi &cell_index);
};

class InitializeIndexMesh : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeIndexMesh(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~InitializeIndexMesh() {};

    void update(const size_t &index);
};

class InitializeCellNeighborhood : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeCellNeighborhood(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~InitializeCellNeighborhood() {};

    void update(const size_t &index);
};

class InitializeBasicDataForAPackage : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeBasicDataForAPackage(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          shape_(shape) {};
    virtual ~InitializeBasicDataForAPackage() {};

    void update(const size_t &index);

  private:
    Shape &shape_;
};

class UpdateLevelSetGradient : public BaseMeshLocalDynamics
{
  public:
    explicit UpdateLevelSetGradient(MeshWithGridDataPackages<4> &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~UpdateLevelSetGradient() {};

    void update(const size_t &index);
};

class UpdateKernelIntegrals : public BaseMeshLocalDynamics
{
  public:
    explicit UpdateKernelIntegrals(MeshWithGridDataPackagesType &mesh_data, Kernel &kernel, Real global_h_ratio)
        : BaseMeshLocalDynamics(mesh_data),
          kernel_(kernel),
          global_h_ratio_(global_h_ratio) {};
    virtual ~UpdateKernelIntegrals() {};

    void update(const size_t &package_index);

  private:
    Kernel &kernel_;
    Real global_h_ratio_;

    Real probeSignedDistance(const Vecd &position) { return mesh_data_.probeMesh(phi_, position); };
    Real computeKernelIntegral(const Vecd &position);
    Vecd computeKernelGradientIntegral(const Vecd &position);
    Matd computeKernelSecondGradientIntegral(const Vecd &position);

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

class ReinitializeLevelSet : public BaseMeshLocalDynamics
{
  public:
    explicit ReinitializeLevelSet(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ReinitializeLevelSet() {};

    void update(const size_t &package_index);

  private:
    // upwind algorithm choosing candidate difference by the sign
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

class MarkNearInterface : public BaseMeshLocalDynamics
{
  public:
    explicit MarkNearInterface(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~MarkNearInterface() {};

    void update(const size_t &package_index);
    void setSmallShiftFactor(Real small_shift_factor) { small_shift = data_spacing_ * small_shift_factor; };

  private:
    Real small_shift;
};

class RedistanceInterface : public BaseMeshLocalDynamics
{
  public:
    explicit RedistanceInterface(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~RedistanceInterface() {};

    void update(const size_t &package_index);
};

class DiffuseLevelSetSign : public BaseMeshLocalDynamics
{
  public:
    explicit DiffuseLevelSetSign(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~DiffuseLevelSetSign() {};

    void update(const size_t &package_index);
};

class InitializeDataInACellFromCoarse : public BaseMeshLocalDynamics
{
  public:
    explicit InitializeDataInACellFromCoarse(MeshWithGridDataPackagesType &mesh_data, MeshWithGridDataPackagesType &coarse_mesh, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          coarse_mesh_(coarse_mesh),
          shape_(shape) {};
    virtual ~InitializeDataInACellFromCoarse() {};

    void update(const Arrayi &cell_index);

  private:
    MeshWithGridDataPackagesType &coarse_mesh_;
    Shape &shape_;
};

class ProbeSignedDistance : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeSignedDistance(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeSignedDistance() {};

    Real update(const Vecd &position) { return mesh_data_.probeMesh(phi_, position); };
};

class ProbeLevelSetGradient : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeLevelSetGradient(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeLevelSetGradient() {};

    Vecd update(const Vecd &position) { return mesh_data_.probeMesh(phi_gradient_, position); };
};

class ProbeKernelIntegral : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeKernelIntegral(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeKernelIntegral() {};

    Real update(const Vecd &position) { return mesh_data_.probeMesh(kernel_weight_, position); };
};

class ProbeKernelGradientIntegral : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeKernelGradientIntegral(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeKernelGradientIntegral() {};

    Vecd update(const Vecd &position) { return mesh_data_.probeMesh(kernel_gradient_, position); };
};

class ProbeKernelSecondGradientIntegral : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeKernelSecondGradientIntegral(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~ProbeKernelSecondGradientIntegral() {};

    Matd update(const Vecd &position) { return mesh_data_.probeMesh(kernel_second_gradient_, position); };
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

class IsWithinCorePackage : public BaseMeshLocalDynamics
{
  public:
    explicit IsWithinCorePackage(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data) {};
    virtual ~IsWithinCorePackage() {};

    bool update(Vecd position)
    {
        Arrayi cell_index = mesh_data_.CellIndexFromPosition(position);
        return mesh_data_.isCoreDataPackage(cell_index);
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
