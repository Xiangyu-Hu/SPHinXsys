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
 * @file    level_set_correction.h
 * @brief   This class encapsulates functions related to mesh dynamics,
 *          including initialization and manipulation of level set data.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef LEVEL_SET_CORRECTION_H
#define LEVEL_SET_CORRECTION_H

#include "base_local_mesh_dynamics.h"
#include "level_set_transformation.h"
#include "mesh_dynamics_algorithm.h"

namespace SPH
{
class ReinitializeLevelSet : public BaseMeshLocalDynamics
{
  public:
    explicit ReinitializeLevelSet(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
    virtual ~ReinitializeLevelSet() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        PackageVariableData<Real> *phi_;
        PackageVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;

        Real upwindDifference(Real sign, Real df_p, Real df_n);
    };

  protected:
    PackageVariable<Real> &mv_phi_;
    PackageVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class MarkCutInterfaces : public BaseMeshLocalDynamics
{
  public:
    explicit MarkCutInterfaces(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
        Real perturbation_ratio);
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
        PackageVariableData<Real> *phi_;
        PackageVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    Real perturbation_ratio_;
    PackageVariable<Real> &mv_phi_;
    PackageVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class MarkNearInterface : public BaseMeshLocalDynamics
{
  public:
    explicit MarkNearInterface(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
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
        PackageVariableData<Real> *phi_;
        PackageVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    PackageVariable<Real> &mv_phi_;
    PackageVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class RedistanceInterface : public BaseMeshLocalDynamics
{
  public:
    RedistanceInterface(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level);
    virtual ~RedistanceInterface() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        Real data_spacing_;
        PackageVariableData<Real> *phi_;
        PackageVariableData<Vecd> *phi_gradient_;
        PackageVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
    };

  protected:
    PackageVariable<Real> &mv_phi_;
    PackageVariable<Vecd> &mv_phi_gradient_;
    PackageVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
};

class DiffuseLevelSetSign : public BaseMeshLocalDynamics
{
  public:
    explicit DiffuseLevelSetSign(
        SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
        SingularVariable<UnsignedInt> &sv_count_modified);
    virtual ~DiffuseLevelSetSign() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(const UnsignedInt &index);

      protected:
        PackageVariableData<Real> *phi_;
        PackageVariableData<int> *near_interface_id_;
        CellNeighborhood *cell_neighborhood_;
        UnsignedInt *count_modified_;
    };

  protected:
    PackageVariable<Real> &mv_phi_;
    PackageVariable<int> &mv_near_interface_id_;
    DiscreteVariable<CellNeighborhood> &dv_cell_neighborhood_;
    SingularVariable<UnsignedInt> &sv_count_modified_;
};

class RepeatTimes
{
  public:
    explicit RepeatTimes() {};
    virtual ~RepeatTimes() {};
    void operator()(UnsignedInt repeat_times) { repeat_times_ = repeat_times; }

  protected:
    UnsignedInt repeat_times_ = 1;
};

template <class ExecutionPolicy>
class CleanInterface : public RepeatTimes, public BaseDynamics<void>
{
  public:
    CleanInterface(SparseMeshField<4> &mesh_data, UnsignedInt resolution_level,
                   NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method,
                   Real refinement_ratio);
    virtual ~CleanInterface() {};

    void exec(Real dt = 0.0) override
    {
        for (UnsignedInt k = 0; k != 2 * repeat_times_; ++k)
        {
            for (UnsignedInt l = 0; l != 2; ++l)
            {
                mark_cut_interfaces.exec();
                redistance_interface.exec();
            }

            for (UnsignedInt m = 0; m != 10; ++m)
            {
                reinitialize_level_set.exec();
            }
        }
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method_;
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient;
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals;
    MeshInnerDynamics<ExecutionPolicy, MarkCutInterfaces> mark_cut_interfaces;
    MeshCoreDynamics<ExecutionPolicy, RedistanceInterface> redistance_interface;
    MeshInnerDynamics<ExecutionPolicy, ReinitializeLevelSet> reinitialize_level_set;
};

template <class ExecutionPolicy>
class CorrectTopology : public BaseDynamics<void>
{
  public:
    CorrectTopology(
        SparseMeshField<4> &mesh_data, UnsignedInt resolution_level,
        NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method);
    virtual ~CorrectTopology() {};

    void exec(Real dt = 0.0) override
    {
        mark_near_interface.exec();
        while (sv_count_modified_.getValue() > 0)
        {
            sv_count_modified_.setValue(0);
            diffuse_level_set_sign.exec();
        }
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method_;
    SingularVariable<UnsignedInt> sv_count_modified_{"CountModifiedData", 1};
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient;
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals;
    MeshInnerDynamics<ExecutionPolicy, MarkNearInterface> mark_near_interface;
    MeshInnerDynamics<ExecutionPolicy, DiffuseLevelSetSign> diffuse_level_set_sign;
};
} // namespace SPH
#endif // LEVEL_SET_CORRECTION_H
