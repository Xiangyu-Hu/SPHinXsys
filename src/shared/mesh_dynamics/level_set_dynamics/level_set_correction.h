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
class CleanInterface : public RepeatTimes, public BaseMeshDynamics, public BaseDynamics<void>
{
  public:
    explicit CleanInterface(MeshWithGridDataPackagesType &mesh_data,
                            NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method,
                            Real refinement_ratio)
        : RepeatTimes(), BaseMeshDynamics(mesh_data), BaseDynamics<void>(),
          neighbor_method_(neighbor_method), refinement_ratio_(refinement_ratio) {};
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
    Real refinement_ratio_;
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, neighbor_method_};
    MeshInnerDynamics<ExecutionPolicy, MarkCutInterfaces> mark_cut_interfaces{mesh_data_, 0.5 * refinement_ratio_};
    MeshCoreDynamics<ExecutionPolicy, RedistanceInterface> redistance_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, ReinitializeLevelSet> reinitialize_level_set{mesh_data_};
};

template <class ExecutionPolicy>
class CorrectTopology : public BaseMeshDynamics, public BaseDynamics<void>
{
  public:
    explicit CorrectTopology(MeshWithGridDataPackagesType &mesh_data, NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method)
        : BaseMeshDynamics(mesh_data),
          BaseDynamics<void>(),
          neighbor_method_(neighbor_method) {};
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
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, neighbor_method_};
    MeshInnerDynamics<ExecutionPolicy, MarkNearInterface> mark_near_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, DiffuseLevelSetSign> diffuse_level_set_sign{mesh_data_, sv_count_modified_};
};
} // namespace SPH
#endif // LEVEL_SET_CORRECTION_H
