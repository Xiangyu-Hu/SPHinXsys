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
 * @file 	sphinxsys_sycl.h
 * @brief 	All SPHinXsys capabilities and SYCL.
 * @author	Xiangyu Hu
 */
#ifndef SPHINXSYS_SYCL_H
#define SPHINXSYS_SYCL_H

#include "base_configuration_dynamics_sycl.h"
#include "particle_iterators_sycl.h"
#include "particle_sort_sycl.h"
#include "particle_sort_sycl.hpp"
#include "device_copyable_variable.h"
#include "all_mesh_dynamics_sycl.h"
#include "sphinxsys_ck.h"
#include "sphinxsys_constant_sycl.hpp"
#include "sphinxsys_variable_array_sycl.hpp"
#include "sphinxsys_variable_sycl.hpp"

#endif // SPHINXSYS_SYCL_H
