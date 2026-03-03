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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	device_copyable_variable.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */
#ifndef DEVICE_COPYABLE_VARIABLE_H
#define DEVICE_COPYABLE_VARIABLE_H

#include "base_data_type_package.h"
#include "simtk_wrapper.h"

namespace sycl
{
template <>
struct is_device_copyable<SPH::SimTKVec3> : std::true_type
{
};

template <>
struct is_device_copyable<SimTK::Vec<2, SPH::SimTKVec3>> : std::true_type
{
};

template <int N, int M>
struct is_device_copyable<Eigen::Matrix<SPH::Real, N, M>> : std::true_type
{
};

template <>
struct is_device_copyable<SPH::Mat2d> : std::true_type
{
};
} // namespace sycl
#endif // DEVICE_COPYABLE_VARIABLE_H
