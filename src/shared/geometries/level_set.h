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
#include "base_geometry.h"
#include "mesh_with_data_packages.hpp"
#include "mesh_dynamics.h"
#include "mesh_local_dynamics.h"
#include "all_mesh_dynamics.h"

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
    MultilevelLevelSet(BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data, Shape &shape, SPHAdaptation &sph_adaptation);
    ~MultilevelLevelSet(){};

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
    StdVec<MeshWithGridDataPackagesType *> getMeshLevels() { return mesh_data_set_; };

    void writeMeshFieldToPlt(std::ofstream &output_file) override
    {
        for (size_t l = 0; l != total_levels_; ++l)
        {
            MeshCalculateDynamics<void, WriteMeshFieldToPlt> write_mesh_field_to_plt(*mesh_data_set_[l]);
            write_mesh_field_to_plt.exec(output_file);
        }
    }

  protected:
    inline size_t getProbeLevel(const Vecd &position);
    inline size_t getCoarseLevel(Real h_ratio);

    Kernel &kernel_;
    Shape &shape_;                           /**< the geometry is described by the level set. */
    size_t total_levels_;                    /**< level 0 is the coarsest */
    StdVec<Real> global_h_ratio_vec_;
    StdVec<MeshWithGridDataPackagesType *> mesh_data_set_;
    StdVec<MeshCalculateDynamics<Real, ProbeSignedDistance> *> probe_signed_distance_set_;
    StdVec<MeshCalculateDynamics<Vecd, ProbeNormalDirection> *> probe_normal_direction_set_;
    StdVec<MeshCalculateDynamics<Vecd, ProbeLevelSetGradient> *> probe_level_set_gradient_set_;
    UniquePtrsKeeper<MeshWithGridDataPackagesType> mesh_data_ptr_vector_keeper_;
    UniquePtrsKeeper<MeshCalculateDynamics<Real, ProbeSignedDistance>> probe_signed_distance_vector_keeper_;
    UniquePtrsKeeper<MeshCalculateDynamics<Vecd, ProbeNormalDirection>> probe_normal_direction_vector_keeper_;
    UniquePtrsKeeper<MeshCalculateDynamics<Vecd, ProbeLevelSetGradient>> probe_level_set_gradient_vector_keeper_;

    CleanInterface *clean_interface;
    CorrectTopology *correct_topology;
};
} // namespace SPH
#endif // LEVEL_SET_H
