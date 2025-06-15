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
 * @file    all_general_dynamics_ck.h
 * @brief   This is the header file that user code should include to pick up all
 *          general dynamics used in SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "adapt_criterion.h"
#include "adapt_indication.h"
#include "all_surface_indication_ck.h"
#include "force_prior_ck.hpp"
#include "general_constraint_ck.h"
#include "general_gradient.hpp"
#include "general_assignment.h"
#include "general_reduce_ck.hpp"
#include "geometric_dynamics.hpp"
#include "hessian_correction_ck.hpp"
#include "interpolation_dynamics.hpp"
#include "kernel_correction_ck.hpp"