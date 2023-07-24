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
 * @file 	scalar_functions.h
 * @brief 	Base scalar function definition.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef SCALAR_FUNCTIONS_H
#define SCALAR_FUNCTIONS_H

#include "base_data_type.h"

#include <chrono>
#include <random>

namespace SPH
{

template <class T>
inline T SMIN(T a1, T a2)
{
    return (a1 <= a2 ? a1 : a2);
}

template <class T>
inline T SMIN(T a1, T a2, T a3)
{
    return SMIN(a1, SMIN(a2, a3));
}

template <class T>
inline T SMIN(T a1, T a2, T a3, T a4)
{
    return SMIN(SMIN(a1, a2), SMIN(a3, a4));
}

template <class T>
inline T SMIN(T a1, T a2, T a3, T a4, T a5)
{
    return SMIN(SMIN(a1, a2), SMIN(a3, a4), a5);
}

template <class T>
inline T SMIN(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return SMIN(SMIN(a1, a2), SMIN(a3, a4), SMIN(a5, a6));
}

template <class T>
inline T SMAX(T a1, T a2)
{
    return (a1 >= a2 ? a1 : a2);
}

template <class T>
inline T SMAX(T a1, T a2, T a3)
{
    return SMAX(a1, SMAX(a2, a3));
}

template <class T>
inline T SMAX(T a1, T a2, T a3, T a4)
{
    return SMAX(SMAX(a1, a2), SMAX(a3, a4));
}

template <class T>
inline T SMAX(T a1, T a2, T a3, T a4, T a5)
{
    return SMAX(SMAX(a1, a2), SMAX(a3, a4), a5);
}

template <class T>
inline T SMAX(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return SMAX(SMAX(a1, a2), SMAX(a3, a4), SMAX(a5, a6));
}

template <class T>
inline void update_minmax(T a1, T &amin, T &amax)
{
    amin = SMIN(a1, amin);
    amax = SMAX(a1, amax);
}

template <class T>
inline void update_minmax(T a1, T a2, T &amin, T &amax)
{
    if (a1 > a2)
    {
        amin = a2;
        amax = a1;
    }
    else
    {
        amin = a1;
        amax = a2;
    }
}

template <class T>
inline T ABS(const T &x)
{
    return SMAX(x, -x);
}

template <class T>
inline T SGN(const T &x)
{
    return (x < 0) ? -1 : ((x > 0) ? 1 : 0);
}

template <class T>
inline T clamp(T a, T lower, T upper)
{
    if (a < lower)
        return lower;
    else if (a > upper)
        return upper;
    else
        return a;
}

template <class T>
inline bool Not_a_number(T a)
{
    return (std::isnan(a) || !(std::isfinite(a))) ? true : false;
}

inline Real rand_norm(Real u, Real std)
{
    unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<Real> distribution(u, std);
    return distribution(generator);
}
/** rotating axis once according to right hand rule.
 * The first_axis must be 0, 1 for 2d and 0, 1, 2 for 3d
 */
int SecondAxis(int first_axis);
/** rotating axis twice  according to right hand rule.
 * The first_axis must be 0, 1 for 2d and 0, 1, 2 for 3d
 */
int ThirdAxis(int first_axis);

Real getLeftStateInWeno(Real v1, Real v2, Real v3, Real v4);
Real getRightStateInWeno(Real v1, Real v2, Real v3, Real v4);

/** linear heaviside function.*/
Real Heaviside(Real phi, Real half_width);
} // namespace SPH
#endif // SCALAR_FUNCTIONS_H
