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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	sphinxsys.h
 * @brief 	This is the header file that user code should include to pick up all SPHinXsys capabilities.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef SPHINXSYS_H
#define SPHINXSYS_H

#if SPHINXSYS_USE_SYCL
#include "base_configuration_dynamics_sycl.h"
#include "device_copyable_variable.h"
#include "mesh_iterators_sycl.hpp"
#include "particle_iterators_sycl.h"
#include "sphinxsys_buffer_array_sycl.hpp"
#include "sphinxsys_constant_sycl.hpp"
#include "sphinxsys_variable_array_sycl.hpp"
#include "sphinxsys_variable_sycl.hpp"
#endif // SPHINXSYS_USE_SYCL

#include "all_bodies.h"
#include "all_body_relations.h"
#include "all_closures.h"
#include "all_geometries.h"
#include "all_io.h"
#include "all_kernels.h"
#include "all_particle_generators.h"
#include "all_particles.h"
#include "all_physical_dynamics.h"
#include "all_regression_test_methods.h"
#include "all_simbody.h"
#include "parameterization.h"
#include "particle_method_container.h"
#include "sph_solver.h"
#include "sph_system.hpp"

#endif // SPHINXSYS_H
