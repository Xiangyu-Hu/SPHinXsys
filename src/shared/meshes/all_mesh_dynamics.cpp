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
 * @file 	base_relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "all_mesh_dynamics.h"

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
// void DefineShapeOnMeshWithGridDataPackages::exec()
// {
//     initialize_data_in_a_cell.exec();
//     tag_a_cell_is_inner_package.exec();
//     initialize_index_mesh.exec();
//     initialize_cell_neighborhood.exec();
//     mesh_data_.resizeMeshVariableData();
// }
//=================================================================================================//
// void CleanInterface::exec()
// {
//     mark_near_interface.exec();
//     redistance_interface.exec();
//     reinitialize_level_set.exec();
//     update_level_set_gradient.exec();
//     update_kernel_integrals.exec();
// }
//=================================================================================================//
// void CorrectTopology::exec()
// {
//     mark_near_interface.exec();
//     for (size_t i = 0; i != 10; ++i)
//         diffuse_level_set_sign.exec();
//     update_level_set_gradient.exec();
//     update_kernel_integrals.exec();
// }
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
