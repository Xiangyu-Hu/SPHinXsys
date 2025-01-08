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
 * @file 	base_implementation.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef BASE_IMPLEMENTATION_H
#define BASE_IMPLEMENTATION_H

#include "sphinxsys_containers.h"

namespace SPH
{
namespace execution
{
template <typename... T>
class Implementation;

template <>
class Implementation<Base>
{
  public:
    explicit Implementation() {}
    ~Implementation() {}

    bool isUpdated() { return is_updated_; };
    void resetUpdated() { is_updated_ = false; };

  protected:
    bool is_updated_ = false;
    void setUpdated() { is_updated_ = true; };
};

} // namespace execution
} // namespace SPH
#endif // BASE_IMPLEMENTATION_H
