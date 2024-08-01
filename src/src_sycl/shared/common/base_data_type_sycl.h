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
 * @file 	base_data_type_sycl.h
 * @brief 	This is the base data type definition for SPHinXsysSYCL.
 * @details For simplicity, the same data types are used in both CPU and device.
 * The only constraint is that eigen matrix and vector types are degraded so that
 * they can be used in device.
 * @author	Xiangyu Hu
 */

#ifndef BASE_DATA_TYPE_SYCL_H
#define BASE_DATA_TYPE_SYCL_H

#include <CL/sycl.hpp>

#ifndef SYCL_DEVICE_ONLY
#define SYCL_DEVICE_ONLY // Eigen matrix and vector types degraded for SYCL
#endif

#include "base_data_type.h"
namespace SPH
{
} // namespace SPH

#endif // BASE_DATA_TYPE_SYCL_H
