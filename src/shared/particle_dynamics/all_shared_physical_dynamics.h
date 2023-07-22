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
 * @file    all_shared_physical_dynamics.h
 * @brief   Head file for all shared physics dynamics for both 2- and 3D build.
 *          This is the header file that user code should include to pick up all
            particle dynamics capabilities.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_SHARED_PHYSICAL_DYNAMICS_H
#define ALL_SHARED_PHYSICAL_DYNAMICS_H

#include "active_muscle_dynamics.h"
#include "all_diffusion_reaction_dynamics.h"
#include "all_fluid_dynamics.h"
#include "all_general_dynamics.h"
#include "all_solid_dynamics.h"
#include "electro_physiology.h"
#include "external_force.h"
#include "particle_dynamics_dissipation.h"
#include "particle_dynamics_dissipation.hpp"
#include "relax_dynamics.h"

#endif // ALL_SHARED_PHYSICAL_DYNAMICS_H
