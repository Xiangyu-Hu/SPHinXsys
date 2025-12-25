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
 * @file 	level_set.h
 * @brief 	Level set is a function which is defined as signed distance to a surface or interface.
 * @author	Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */

#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include "adaptation.h"
#include "base_geometry.h"
#include "kernel_tabulated_ck.h"
#include "level_set_correction.hpp"
#include "level_set_initialization.hpp"
#include "level_set_transformation.hpp"
#include "mesh_data_package_sort.h"
#include "mesh_dynamics_algorithm.h"
#include "sphinxsys_variable.h"
namespace SPH
{
enum class UsageType
{
    Volumetric,
    Surface,
};
/**
 * @class LevelSet
 * @brief Defining a multilevel level set for a complex region.
 */
class LevelSet : public SparseMeshField<4>
{
  public:
    LevelSet(BoundingBoxd tentative_bounds, Real reference_data_spacing, size_t total_levels,
             Shape &shape, const SPHAdaptation &sph_adaptation, Real refinement = 1.0);
    LevelSet(BoundingBoxd tentative_bounds, SparseMeshField<4> *coarse_data,
             Shape &shape, const SPHAdaptation &sph_adaptation, Real refinement = 1.0);
    ~LevelSet() {};

    template <class ExecutionPolicy>
    void finishInitialization(const ExecutionPolicy &ex_policy, UsageType usage_type);
    void cleanInterface(UnsignedInt repeat_times);
    void correctTopology();
    Real probeSignedDistance(const Vecd &position);
    Vecd probeNormalDirection(const Vecd &position);
    Vecd probeLevelSetGradient(const Vecd &position);
    Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0);
    Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0);
    Matd probeKernelSecondGradientIntegral(const Vecd &position, Real h_ratio = 1.0);
    void writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence = 0) override;

    template <typename DataType>
    class ProbeLevelSet : public SparseMeshField<4>::ProbeMesh<DataType>
    {
        using BaseProbe = SparseMeshField<4>::ProbeMesh<DataType>;
        Real *global_h_ratio_;

      public:
        template <class ExecutionPolicy>
        ProbeLevelSet(const ExecutionPolicy &ex_policy, LevelSet &encloser, const std::string &variable_name);
        DataType operator()(const Vecd &position);
        DataType operator()(const Vecd &position, Real h_ratio);
    };

  protected:
    void initializeLevel(UnsignedInt level, SparseMeshField<4> *coarse_data = nullptr, UnsignedInt coarse_level = 0);
    template <class ExecutionPolicy>
    void initializePackageVariables(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void initializeKernelIntegralVariables(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void registerProbes(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void registerKernelIntegralProbes(const ExecutionPolicy &ex_policy);

    Shape &shape_; /**< the geometry is described by the level set. */
    Real refinement_;
    ConstantArray<Real> *ca_global_h_ratio_; /**< the ratio of the reference spacing to the data spacing */
    StdVec<NeighborMethod<SPHAdaptation, SPHAdaptation> *> neighbor_method_set_;
    ProbeLevelSet<Real> *probe_signed_distance_;
    ProbeLevelSet<Vecd> *probe_level_set_gradient_;
    ProbeLevelSet<Real> *probe_kernel_integral_;
    ProbeLevelSet<Vecd> *probe_kernel_gradient_integral_;
    ProbeLevelSet<Matd> *probe_kernel_second_gradient_integral_;
    UniquePtrKeeper<ProbeLevelSet<Real>> probe_signed_distance_keeper_;
    UniquePtrKeeper<ProbeLevelSet<Vecd>> probe_level_set_gradient_keeper_;
    UniquePtrKeeper<ProbeLevelSet<Real>> probe_kernel_integral_keeper_;
    UniquePtrKeeper<ProbeLevelSet<Vecd>> probe_kernel_gradient_integral_keeper_;
    UniquePtrKeeper<ProbeLevelSet<Matd>> probe_kernel_second_gradient_integral_keeper_;

    UniquePtr<BaseDynamics<void>> correct_topology_keeper_;
    UniquePtr<BaseDynamics<void>> clean_interface_keeper_;
    UniquePtrsKeeper<NeighborMethod<SPHAdaptation, SPHAdaptation>> neighbor_method_keeper_;
    std::function<void()> sync_mesh_variables_to_write_, sync_mesh_variables_to_probe_;

    template <class ExecutionPolicy>
    void configLevelSetPostProcesses(const ExecutionPolicy &ex_policy);
};
} // namespace SPH
#endif // LEVEL_SET_H
