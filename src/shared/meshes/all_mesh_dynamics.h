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
 * @file 	all_mesh_dynamics.h
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_MESH_DYNAMICS_H
#define ALL_MESH_DYNAMICS_H

#include "mesh_dynamics.h"
#include "mesh_local_dynamics.h"

namespace SPH
{
class FinishDataPackages : public BaseMeshDynamics
{
  public:
    explicit FinishDataPackages(MeshWithGridDataPackages<4> &mesh_data, Shape &shape, Kernel &kernel, Real global_h_ratio)
        : BaseMeshDynamics(mesh_data),
          shape_(shape),
          kernel_(kernel),
          global_h_ratio_(global_h_ratio),
          grid_spacing_(mesh_data.GridSpacing()),
          buffer_width_(mesh_data.BufferWidth()){};
    virtual ~FinishDataPackages(){};

    void exec(){
        tag_a_cell_is_inner_package.exec();

        mesh_data_.organizeOccupiedPackages();
        initialize_index_mesh.exec();
        initialize_cell_neighborhood.exec();
        mesh_data_.resizeMeshVariableData();

        Real far_field_distance = grid_spacing_ * (Real)buffer_width_;
        initialize_data_for_singular_package.exec(0, -far_field_distance);
        initialize_data_for_singular_package.exec(1, far_field_distance);

        initialize_basic_data_for_a_package.exec();
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    };

  private:
    Shape &shape_;
    Kernel &kernel_;
    Real global_h_ratio_;
    Real grid_spacing_;
    size_t buffer_width_;

    MeshSingleDynamics<InitializeDataForSingularPackage> initialize_data_for_singular_package{mesh_data_};
    MeshAllDynamics<TagACellIsInnerPackage> tag_a_cell_is_inner_package{mesh_data_};
    MeshInnerDynamics<InitializeIndexMesh> initialize_index_mesh{mesh_data_};
    MeshInnerDynamics<InitializeCellNeighborhood> initialize_cell_neighborhood{mesh_data_};
    MeshInnerDynamics<InitializeBasicDataForAPackage> initialize_basic_data_for_a_package{mesh_data_, shape_};
    MeshInnerDynamics<UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, kernel_, global_h_ratio_};
};

class ProbeNormalDirection : public BaseMeshLocalDynamics
{
  public:
    explicit ProbeNormalDirection(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data){};
    virtual ~ProbeNormalDirection(){};

    Vecd update(const Vecd &position)
    {
        Vecd probed_value = probe_level_set_gradient.exec(position);

        Real threshold = 1.0e-2 * data_spacing_;
        while (probed_value.norm() < threshold)
        {
            Vecd jittered = position; // jittering
            for (int l = 0; l != position.size(); ++l)
                jittered[l] += rand_uniform(-0.5, 0.5) * 0.5 * data_spacing_;
            probed_value = probe_level_set_gradient.exec(jittered);
        }
        return probed_value.normalized();
    }

  private:
    MeshCalculateDynamics<Vecd, ProbeLevelSetGradient> probe_level_set_gradient{mesh_data_};
};

// class CleanInterface
//     : BaseMeshDynamics
// {
//   public:
//     explicit CleanInterface(MeshWithGridDataPackages &mesh_data, Real small_shift_factor)
//         : BaseMeshDynamics(mesh_data),
//           small_shift_factor_(small_shift_factor){};
//     virtual ~CleanInterface(){};

//     virtual void exec() override;

//   private:
//     Real small_shift_factor_;

//     MeshInnerDynamics<MarkNearInterface> mark_near_interface(small_shift_factor_);
//     MeshCoreDynamics<RedistanceInterface> redistance_interface;
//     MeshInnerDynamics<ReinitializeLevelSet> reinitialize_level_set;
//     MeshInnerDynamics<UpdateLevelSetGradient> update_level_set_gradient;
//     MeshInnerDynamics<UpdateKernelIntegrals> update_kernel_integrals;
// }

// class CorrectTopology
//     : BaseMeshDynamics
// {
//   public:
//     explicit CorrectTopology(MeshWithGridDataPackages &mesh_data, Real small_shift_factor)
//         : BaseMeshDynamics(mesh_data),
//           small_shift_factor_(small_shift_factor){};
//     virtual ~CorrectTopology(){};

//     virtual void exec() override;

//   private:
//     Real small_shift_factor_;

//     MeshInnerDynamics<MarkNearInterface> mark_near_interface(small_shift_factor_);
//     MeshInnerDynamics<DiffuseLevelSetSign> diffuse_level_set_sign();
//     MeshInnerDynamics<UpdateLevelSetGradient> update_level_set_gradient;
//     MeshInnerDynamics<UpdateKernelIntegrals> update_kernel_integrals;
// }
} // namespace SPH
#endif // ALL_MESH_DYNAMICS_H