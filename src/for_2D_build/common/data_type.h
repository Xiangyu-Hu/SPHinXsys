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
 * @file 	data_type.h
 * @brief 	This is the date type definition in 2D for SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef DATA_TYPE_2D_H
#define DATA_TYPE_2D_H

#include "base_data_type.h"
#include "geometric_primitive.h"
#include "scalar_functions.h"

namespace SPH
{
using Arrayi = Array2i;
using Vecd = Vec2d;
using Matd = Mat2d;
using VecMatd = Vec3d;           // vectorized symmetric 2x2 matrix
using MatTend = Mat3d;           // matricized symmetric 2x2x2x2 tensor
using VecMatGrad = VecMatGrad2d; // gradient of vectorized symmetric 2x2 matrix
using AngularVecd = Real;
using Rotation = Rotation2d;
using BoundingBoxd = BoundingBox<VecdBound, 2>;
using BoundingBoxi = BoundingBox<ArrayiBound, 2>;
using Transform = BaseTransform<Rotation2d, Vec2d>;

/** only works for smoothing length ratio less or equal than 1.3*/
constexpr int MaximumNeighborhoodSize = int(M_PI * 9);
constexpr int Dimensions = 2;
/** correction matrix, only works for thin structure dynamics. */
const Matd reduced_unit_matrix{
    {1.0, 0.0}, // First row
    {0.0, 0.0}, // Second row
};

/** initial local normal, only works for thin structure dynamics. */
const Vecd local_pseudo_n_0 = Vecd(0.0, 1.0);
const Vecd ZeroVecd = Vec2d::Zero();

inline Vecd degradeToVecd(const Vec3d &input) { return Vecd(input[0], input[1]); };
inline Matd degradeToMatd(const Mat3d &input) { return input.block<2, 2>(0, 0); };

} // namespace SPH

#endif // DATA_TYPE_2D_H