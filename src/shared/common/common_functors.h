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
 * @file    common_functors.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef COMMON_FUNCTORS_H
#define COMMON_FUNCTORS_H

#include "base_data_type.h"
#include "scalar_functions.h"

namespace SPH
{

template <typename... Args>
struct ConstructArgs : public std::tuple<Args...>
{
    ConstructArgs(Args... args) : std::tuple<Args...>(args...) {}
};

template <class DataType>
struct ReturnFunction
{
    typedef DataType ReturnType;
};

struct AssignIndex
{
    UnsignedInt operator()(UnsignedInt i) const { return i; }
};

/**
 * @class Limiter
 * Base class introduce the concept of limiter,
 * which limits the magnitude of a value with a fraction.
 * The derived class should implement the operator Real()
 * to indicate limiting fraction.
 * Generally, the object of the derived class
 * should be named as "limiter" or "limiter_" (class member)
 * so that the code can be more readable.
 */
class Limiter
{
};

class NoLimiter : public Limiter
{
  public:
    template <typename... Args>
    NoLimiter(Args &&...args) : Limiter(){};

    template <typename... Args>
    Real operator()(Args &&...args)
    {
        return 1.0;
    };
};

class TruncatedLinear : public Limiter
{
    Real ref_, slope_;

  public:
    TruncatedLinear(Real ref, Real slope = 100.0)
        : Limiter(), ref_(ref), slope_(slope) {};
    Real operator()(Real measure)
    {
        Real measure_scale = measure * ref_;
        return SMIN(slope_ * measure_scale, Real(1));
    };
};

} // namespace SPH
#endif // COMMON_FUNCTORS_H
