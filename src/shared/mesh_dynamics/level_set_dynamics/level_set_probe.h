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
 * @file level_set_probe.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef LEVEL_SET_PROBE_H
#define LEVEL_SET_PROBE_H

#include "base_local_mesh_dynamics.h"

namespace SPH
{
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
} // namespace SPH
#endif // LEVEL_SET_PROBE_H
