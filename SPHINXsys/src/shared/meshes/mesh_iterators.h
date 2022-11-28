/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4	.                                                *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file mesh_iterators.h
 * @brief This is for the base functions for mesh iterator.
 * @author  Xiangyu Hu
 */

#ifndef MESH_ITERATORS_H
#define MESH_ITERATORS_H

#include "base_data_package.h"

namespace SPH
{
    using MeshRange = std::pair<Vecu, Vecu>;

    /** Iterator on the mesh by looping index. sequential computing. */
    template <typename LocalFunction, typename... Args>
    void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);
    /** Iterator on the mesh by looping index. parallel computing. */
    template <typename LocalFunction, typename... Args>
    void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);
}
#endif // MESH_ITERATORS_H