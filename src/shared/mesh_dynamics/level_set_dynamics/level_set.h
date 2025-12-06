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
class LevelSet : public BaseMeshField
{
  public:
    LevelSet(BoundingBoxd tentative_bounds, Real reference_data_spacing, size_t total_levels,
             Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio = 1.0);
    LevelSet(BoundingBoxd tentative_bounds, MeshWithGridDataPackagesType *coarse_data,
             Shape &shape, SPHAdaptation &sph_adaptation, Real refinement_ratio = 1.0);
    ~LevelSet() {};

    template <class ExecutionPolicy>
    ProbeSignedDistance getProbeSignedDistance(const ExecutionPolicy &ex_policy);

    template <class ExecutionPolicy>
    ProbeNormalDirection getProbeNormalDirection(const ExecutionPolicy &ex_policy);

    template <class ExecutionPolicy>
    ProbeKernelGradientIntegral getProbeKernelGradientIntegral(const ExecutionPolicy &ex_policy);

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
    StdVec<MeshWithGridDataPackagesType *> getMeshLevels() { return mesh_data_set_; };

    template <typename DataType>
    void addMeshVariableToWrite(const std::string &variable_name);
    void writeMeshFieldToPlt(const std::string &partial_file_name, size_t sequence = 0) override;
    template <typename DataType>
    void addBKGMeshVariableToWrite(const std::string &variable_name);
    void writeBKGMeshToPlt(const std::string &partial_file_name) override;
    template <class ExecutionPolicy>
    void syncMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncBKGMeshVariablesToWrite(ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void syncMeshVariablesToProbe(ExecutionPolicy &ex_policy);

  protected:
    inline size_t getProbeLevel(const Vecd &position);
    inline size_t getCoarseLevel(Real h_ratio);

    void initializeLevel(Real reference_data_spacing, BoundingBoxd tentative_bounds,
                         MeshWithGridDataPackagesType *coarse_data = nullptr);
    template <class ExecutionPolicy>
    void initializeMeshVariables(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void initializeKernelIntegralVariables(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void registerProbes(const ExecutionPolicy &ex_policy);
    template <class ExecutionPolicy>
    void registerKernelIntegralProbes(const ExecutionPolicy &ex_policy);

    size_t total_levels_; /**< level 0 is the coarsest */
    Shape &shape_;        /**< the geometry is described by the level set. */
    Real refinement_ratio_;
    StdVec<UnsignedInt *> cell_pkg_index_set_;
    StdVec<int *> pkg_type_set_;
    StdVec<Real> global_h_ratio_vec_; /**< the ratio of the reference spacing to the data spacing */
    StdVec<NeighborMethod<SPHAdaptation, SPHAdaptation> *> neighbor_method_set_;
    StdVec<MeshWithGridDataPackagesType *> mesh_data_set_;
    StdVec<MeshWithGridDataPackagesType::IndexHandler *> mesh_index_handler_set_;
    StdVec<ProbeSignedDistance *> probe_signed_distance_set_;
    StdVec<ProbeNormalDirection *> probe_normal_direction_set_;
    StdVec<ProbeLevelSetGradient *> probe_level_set_gradient_set_;
    StdVec<ProbeKernelIntegral *> probe_kernel_integral_set_;
    StdVec<ProbeKernelGradientIntegral *> probe_kernel_gradient_integral_set_;
    StdVec<ProbeKernelSecondGradientIntegral *> probe_kernel_second_gradient_integral_set_;
    UniquePtrsKeeper<MeshWithGridDataPackagesType> mesh_data_ptr_vector_keeper_;
    UniquePtrsKeeper<ProbeSignedDistance> probe_signed_distance_vector_keeper_;
    UniquePtrsKeeper<ProbeNormalDirection> probe_normal_direction_vector_keeper_;
    UniquePtrsKeeper<ProbeLevelSetGradient> probe_level_set_gradient_vector_keeper_;
    UniquePtrsKeeper<ProbeKernelIntegral> probe_kernel_integral_vector_keeper_;
    UniquePtrsKeeper<ProbeKernelGradientIntegral> probe_kernel_gradient_integral_vector_keeper_;
    UniquePtrsKeeper<ProbeKernelSecondGradientIntegral> probe_kernel_second_gradient_integral_vector_keeper_;

    UniquePtr<BaseDynamics<void>> correct_topology_keeper_;
    UniquePtr<BaseDynamics<void>> clean_interface_keeper_;
    UniquePtrsKeeper<NeighborMethod<SPHAdaptation, SPHAdaptation>> neighbor_method_keeper_;
    std::function<void()> sync_mesh_variables_to_write_, sync_bkg_mesh_variables_to_write_, sync_mesh_variables_to_probe_;

    template <class ExecutionPolicy>
    void configLevelSetPostProcesses(const ExecutionPolicy &ex_policy);
};
} // namespace SPH
#endif // LEVEL_SET_H
