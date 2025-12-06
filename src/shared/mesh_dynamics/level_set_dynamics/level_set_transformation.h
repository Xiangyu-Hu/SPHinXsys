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
 * @file    level_set_transformation.h
 * @brief   This class encapsulates functions related to mesh dynamics,
 *          including initialization and manipulation of level set data.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef LEVEL_SET_TRANSFORMATION_H
#define LEVEL_SET_TRANSFORMATION_H

#include "base_local_mesh_dynamics.h"
#include "neighbor_method.h"

namespace SPH
{
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
        SmoothingKernel kernel_;
        Real data_spacing_, data_cell_volume_;
        CellNeighborhood *cell_neighborhood_;
        Real cutoff_radius_;
        BoundingBoxi bounding_box_;

        Real computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &data_index);
        Vecd computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index);
        Matd computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index);
        /** a cut cell is a cut by the level set. */
        /** "Multi-scale modeling of compressible multi-fluid flows with conservative interface method."
         * Hu, X. Y., et al., Proceedings of the Summer Program. Vol. 301. Stanford, CA, USA:
         * Center for Turbulence Research, Stanford University, 2010.*/
        Real CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing);

        template <typename DataType, typename FunctionByDataIndex>
        DataType computeIntegral(Real phi, const UnsignedInt &package_index, const Arrayi &grid_index,
                                 const DataType &initial_value, const FunctionByDataIndex &function_by_index);
    };

  private:
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method_;
    Real global_h_ratio_;
    MeshVariable<Real> &mv_phi_;
    MeshVariable<Vecd> &mv_phi_gradient_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    MeshVariable<Real> &mv_kernel_weight_;
    MeshVariable<Vecd> &mv_kernel_gradient_;
    MeshVariable<Matd> &mv_kernel_second_gradient_;

    void initializeSingularPackages(UnsignedInt package_index, Real far_field_level_set);
};
} // namespace SPH
#endif // LEVEL_SET_TRANSFORMATION_H
