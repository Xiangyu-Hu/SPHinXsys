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
 * @file 	sphinxsys.h
 * @brief 	This is the header file that user code should include to pick up all SPHinXsys capabilities.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef SPHINXSYS_H
#define SPHINXSYS_H

#include "all_bodies.h"
#include "all_body_relations.h"
#include "all_geometries.h"
#include "all_kernels.h"
#include "all_materials.h"
#include "all_particle_generators.h"
#include "all_particles.h"
#include "all_physical_dynamics.h"
#include "all_simbody.h"
#include "io_all.h"
#include "parameterization.h"
#include "all_regression_test_methods.h"
#include "sph_system.h"

#endif // SPHINXSYS_H
