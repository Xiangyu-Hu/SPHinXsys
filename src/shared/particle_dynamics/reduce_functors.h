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
 * @file    reduce_functors.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef REDUCE_FUNCTORS_H
#define REDUCE_FUNCTORS_H

#include "base_data_package.h"

namespace SPH
{
template <typename...>
struct ReduceReference;

template <class DataType>
struct ReduceSum : ReturnFunction<DataType>
{
    DataType operator()(const DataType &x, const DataType &y) const { return x + y; };
};

template <typename DataType>
struct ReduceReference<ReduceSum<DataType>>
{
    static inline DataType value = ZeroData<DataType>::value;
};

struct ReduceMax : ReturnFunction<Real>
{
    Real operator()(Real x, Real y) const { return SMAX(x, y); };
};

template <>
struct ReduceReference<ReduceMax>
{
    static inline Real value = MinReal;
};

struct ReduceMin : ReturnFunction<Real>
{
    Real operator()(Real x, Real y) const { return SMIN(x, y); };
};

template <>
struct ReduceReference<ReduceMin>
{
    static inline Real value = MaxReal;
};
struct ReduceOR : ReturnFunction<bool>
{
    bool operator()(bool x, bool y) const { return x || y; };
};

template <>
struct ReduceReference<ReduceOR>
{
    static inline bool value = false;
};

struct ReduceAND : ReturnFunction<bool>
{
    bool operator()(bool x, bool y) const { return x && y; };
};

template <>
struct ReduceReference<ReduceAND>
{
    static inline bool value = true;
};
struct ReduceLowerBound : ReturnFunction<Vecd>
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd lower_bound;
        for (int i = 0; i < lower_bound.size(); ++i)
            lower_bound[i] = SMIN(x[i], y[i]);
        return lower_bound;
    };
};

template <>
struct ReduceReference<ReduceLowerBound>
{
    static inline Vecd value = MaxReal * Vecd::Ones();
};

struct ReduceUpperBound : ReturnFunction<Vecd>
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd upper_bound;
        for (int i = 0; i < upper_bound.size(); ++i)
            upper_bound[i] = SMAX(x[i], y[i]);
        return upper_bound;
    };
};

template <>
struct ReduceReference<ReduceUpperBound>
{
    static inline Vecd value = MinReal * Vecd::Ones();
};

} // namespace SPH
#endif // REDUCE_FUNCTORS_H
