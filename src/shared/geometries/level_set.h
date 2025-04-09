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
 * @file 	level_set.h
 * @brief 	Level set is a function which is defined as signed distance to a surface or interface.
 * @author	Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */

#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include "adaptation.h"
#include "all_mesh_dynamics.h"
#include "base_geometry.h"
#include "mesh_dynamics.h"
#include "mesh_local_dynamics.h"
#include "mesh_with_data_packages.hpp"

namespace SPH
{
/**
 * @class MultilevelLevelSet
 * @brief Defining a multilevel level set for a complex region.
 */
class MultilevelLevelSet : public BaseMeshField
{
  public:
    MultilevelLevelSet(BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation);
    MultilevelLevelSet(BoundingBox tentative_bounds, MeshWithGridDataPackagesType *coarse_data, Shape &shape, SPHAdaptation &sph_adaptation);
    ~MultilevelLevelSet() {};

    void cleanInterface(Real small_shift_factor);
    void correctTopology(Real small_shift_factor);
    bool probeIsWithinMeshBound(const Vecd &position);
    Real probeSignedDistance(const Vecd &position);
    Vecd probeNormalDirection(const Vecd &position);
    Vecd probeLevelSetGradient(const Vecd &position);
    Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0);
    Real probeKernelIntegral(const Vecd &position);
    Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0);
    Vecd probeKernelGradientIntegral(const Vecd &position);
    Matd probeKernelSecondGradientIntegral(const Vecd &position, Real h_ratio = 1.0);
    Matd probeKernelSecondGradientIntegral(const Vecd &position);
    StdVec<MeshWithGridDataPackagesType *> getMeshLevels() { return mesh_data_set_; };

    void writeMeshFieldToPlt(const std::string &partial_file_name) override
    {
        for (size_t l = 0; l != total_levels_; ++l)
        {
            std::string full_file_name = partial_file_name + "_" + std::to_string(l) + ".dat";
            std::ofstream out_file(full_file_name.c_str(), std::ios::app);
            WriteMeshFieldToPlt(*mesh_data_set_[l]).update(out_file);
            out_file.close();
        }
    }

  protected:
    inline size_t getProbeLevel(const Vecd &position);
    inline size_t getCoarseLevel(Real h_ratio);

    void initializeLevel(size_t level, Real reference_data_spacing, Real global_h_ratio, BoundingBox tentative_bounds, MeshWithGridDataPackagesType *coarse_data = nullptr);
    void registerProbes(size_t level);

    Kernel &kernel_;
    Shape &shape_;        /**< the geometry is described by the level set. */
    size_t total_levels_; /**< level 0 is the coarsest */
    StdVec<Real> global_h_ratio_vec_;
    StdVec<MeshWithGridDataPackagesType *> mesh_data_set_;
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

    UniquePtr<CleanInterface> clean_interface;
    UniquePtr<CorrectTopology> correct_topology;
};
} // namespace SPH
#endif // LEVEL_SET_H
