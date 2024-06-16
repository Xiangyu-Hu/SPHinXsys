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
#include "all_body_relations.h"
#include "base_body.h"
#include "base_data_package.h"
#include "neighborhood.h"
#include "sph_data_containers.h"

#include <functional>

using namespace std::placeholders;

namespace SPH
{
// class DefineShapeOnMeshWithGridDataPackages
//     : BaseMeshDynamics
// {
//   public:
//     explicit DefineShapeOnMeshWithGridDataPackages(MeshWithGridDataPackages &mesh_data, Shape &shape)
//         : BaseMeshDynamics(mesh_data),
//           shape_(shape){};
//     virtual ~DefineShapeOnMeshWithGridDataPackages(){};

//     virtual void exec() override;

//   private:
//     Shape &shape_;

//     MeshAllDynamics<InitializeDataInACell> initialize_data_in_a_cell(mesh_data_, shape_);
//     MeshAllDynamics<TagACellIsInnerPackage> tag_a_cell_is_inner_package(mesh_data_);
//     MeshInnerDynamics<InitializeIndexMesh>  initialize_index_mesh(mesh_data_);
//     MeshInnerDynamics<InitializeCellNeighborhood> initialize_cell_neighborhood(mesh_data_);

//     bool isInnerPackage();
// };

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