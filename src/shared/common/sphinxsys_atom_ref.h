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
 * @file 	sphinxsys_atom_ref.h
 * @brief 	This is the date type definition for SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPHINXSYS_ATOM_REF_H
#define SPHINXSYS_ATOM_REF_H

#if !SPHINXSYS_USE_SYCL
#include <boost/atomic/atomic_ref.hpp>
#endif // !SPHINXSYS_USE_SYCL

namespace SPH
{
#if SPHINXSYS_USE_SYCL
template <typename T>
using AtomicRef = sycl::atomic_ref<
    T, sycl::memory_order_relaxed, sycl::memory_scope_device,
    sycl::access::address_space::global_space>;
#else
template <typename T>
using AtomicRef = boost::atomic_ref<T>;
#endif // SPHINXSYS_USE_SYCL
} // namespace SPH
#endif // SPHINXSYS_ATOM_REF_H
