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
 * @file    all_fluid_dynamics.h
 * @brief   This is the header file that user code should include to pick up all
 *          fluid dynamics used in SPHinXsys.
 * @details The fluid dynamics algorithms begin for fluid bulk without boundary condition,
 *          then algorithm interacting with wall is defined, further algorithms
 *          for multiphase flow interaction built upon these basic algorithms.
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "all_eulerian_fluid_dynamics.h"
#include "all_fluid_boundaries.h"
#include "density_summation.hpp"
#include "fluid_integration.hpp"
#include "fluid_time_step.h"
#include "non_newtonian_dynamics.h"
#include "shape_confinement.h"
#include "surface_tension.hpp"
#include "transport_velocity_correction.hpp"
#include "velocity_gradient.hpp"
#include "viscous_dynamics.hpp"
