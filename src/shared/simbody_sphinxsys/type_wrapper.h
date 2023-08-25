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
 * @file 	type_wrapper.h
 * @brief 	type wrapper between eigen and SimTK.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef TYPE_WRAPPER_H
#define TYPE_WRAPPER_H

#include "base_data_type.h"
#include "simbody_middle.h"

namespace SPH
{
using SimTKVec2 = SimTK::Vec2;
using SimTKVec3 = SimTK::Vec3;
using SimTKMat22 = SimTK::Mat22;
using SimTKMat33 = SimTK::Mat33;

SimTKVec2 EigenToSimTK(const Vec2d &eigen_vector);
SimTKVec3 EigenToSimTK(const Vec3d &eigen_vector);
Vec2d SimTKToEigen(const SimTKVec2 &simTK_vector);
Vec3d SimTKToEigen(const SimTKVec3 &simTK_vector);

SimTKMat22 EigenToSimTK(const Mat2d &eigen_matrix);
SimTKMat33 EigenToSimTK(const Mat3d &eigen_matrix);
Mat2d SimTKToEigen(const SimTKMat22 &simTK_matrix);
Mat3d SimTKToEigen(const SimTKMat33 &simTK_matrix);

} // namespace SPH
#endif // TYPE_WRAPPER_H
