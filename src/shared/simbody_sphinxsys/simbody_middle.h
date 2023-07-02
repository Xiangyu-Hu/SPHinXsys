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
 * @file 	simbody_middle.h
 * @brief 	file to include Simbody headers and suppress their warnings
 * @author	Wen-Yang Chu
 */

#ifndef SIMBODY_MIDDLE_H
#define SIMBODY_MIDDLE_H

#ifdef __linux__
#pragma GCC system_header // for GCC/CLANG
#elif _WIN32
#pragma warning(push, 0) // for MSVC
#endif

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

#ifdef _WIN32
#pragma warning(pop) // for MSVC
#endif
#endif // SIMBODY_MIDDLE_H