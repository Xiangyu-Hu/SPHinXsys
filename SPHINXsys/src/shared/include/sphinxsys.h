/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/

#ifndef SPHINXSYS_H
#define SPHINXSYS_H

/** @file
This is the header file that user code should include to pick up all SPHinXsys
capabilities. **/

#include "all_kernels.h"
#include "all_particles.h"
#include "all_particle_generators.h"
#include "all_geometries.h"
#include "all_bodies.h"
#include "generative_structures.h"
#include "sph_system.h"
#include "all_materials.h"
#include "all_physical_dynamics.h"
#include "all_simbody.h"
#include "in_output.h"
#include "parameterization.h"
#include "regression_test.h"
#include "all_boundaries.h"

#endif //SPHINXSYS_H
