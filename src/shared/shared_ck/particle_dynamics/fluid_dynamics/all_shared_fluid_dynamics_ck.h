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
 * @file    all_shared_fluid_dynamics_ck.h
 * @brief   Head file for all shared physics dynamics for both 2- and 3D build.
 *          This is the header file that user code should include to pick up all
            particle dynamics capabilities.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_SHARED_FLUID_DYNAMICS_CK_H
#define ALL_SHARED_FLUID_DYNAMICS_CK_H

#include "acoustic_step_1st_half.hpp"
#include "acoustic_step_2nd_half.hpp"
#include "density_regularization.hpp"
#include "all_fluid_boundary_condition_ck.h"
#include "fluid_time_step_ck.hpp"
#include "transport_velocity_correction_ck.hpp"
#include "viscous_force.hpp"

#endif // ALL_SHARED_FLUID_DYNAMICS_CK_H
