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
 * @file    all_shared_physical_dynamics_ck.h
 * @brief   Head file for all shared physics dynamics for both 2- and 3D build.
 *          This is the header file that user code should include to pick up all
            particle dynamics capabilities.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_SHARED_PHYSICAL_DYNAMICS_CK_H
#define ALL_SHARED_PHYSICAL_DYNAMICS_CK_H

#include "all_fluid_structure_interactions.h"
#include "all_general_dynamics_ck.h"
#include "all_shared_fluid_dynamics_ck.h"
#include "all_solid_dynamics_ck.h"
#include "complex_algorithms_ck.h"
#include "diffusion_dynamics_ck.hpp"
#include "interaction_algorithms_ck.hpp"
#include "particle_functors_ck.h"
#include "particle_sort_ck.hpp"
#include "simple_algorithms_ck.h"
#include "all_continum_dynamics.h"
#include "sph_solver.h"

#endif // ALL_SHARED_PHYSICAL_DYNAMICS_CK_H
