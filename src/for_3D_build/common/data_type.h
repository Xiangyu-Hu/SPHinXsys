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
 * @file 	data_type.h
 * @brief 	This is the date type definition in 3D for SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef DATA_TYPE_3D_H
#define DATA_TYPE_3D_H

#include "base_data_type.h"
#include "scalar_functions.h"
#include "vector_functions.h"

namespace SPH
{
using Arrayi = Array3i;
using Vecd = Vec3d;
using Matd = Mat3d;
using AlignedBox = AlignedBox3d;
using AngularVecd = Vec3d;
using Rotation = Rotation3d;
using BoundingBox = BaseBoundingBox<Vec3d>;

template <class DataType, int array_size>
using PackageDataMatrix = std::array<std::array<std::array<DataType, array_size>, array_size>, array_size>;

template <class DataType>
using MeshDataMatrix = DataType ***;

/** only works for smoothing length ratio less or equal than 1.3*/
constexpr int MaximumNeighborhoodSize = int(1.33 * M_PI * 27);
constexpr int Dimensions = 3;
/** correction matrix, only works for thin structure dynamics. */
const Matd reduced_unit_matrix{
    {1, 0, 0}, // 0 row
    {0, 1, 0}, // 1 row
    {0, 0, 0}, // 2 row
};
/** initial local normal, only works for thin structure dynamics. */
const Vecd local_pseudo_n_0 = Vecd(0.0, 0.0, 1.0);
const Vecd local_pseudo_b_n_0 = Vecd(0.0, 1.0, 0.0);

const Vecd ZeroVecd = Vec3d::Zero();
} // namespace SPH
#endif // DATA_TYPE_3D_H