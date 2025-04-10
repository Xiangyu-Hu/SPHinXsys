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
 * @file 	scalar_numerics.h
 * @brief 	Use as scalar for high-level computation using eigen library.
 * @author	Xiangyu Hu
 */
#ifndef SCALAR_NUMERICS_H
#define SCALAR_NUMERICS_H

#include "base_data_type.h"

namespace SPH
{
template <typename T>
struct Scalar
{
    T value;

    // Constructors
    Scalar() = default;
    Scalar(const T &v) : value(v) {}
    Scalar(const Real &s) : value(s * UnitData<T>::value) {}
    // Implicit conversion
    operator const T &() const { return value; }
    operator T &() { return value; }

    // Access
    T &ref() { return value; }
    const T &get() const { return value; }

    // Arithmetic (only if T supports it)
    Scalar &operator+=(const Scalar &rhs)
    {
        value += rhs.value;
        return *this;
    }
    Scalar &operator-=(const Scalar &rhs)
    {
        value -= rhs.value;
        return *this;
    }
    Scalar &operator*=(const typename T::Scalar &s)
    {
        value *= s;
        return *this;
    }
    Scalar &operator/=(const typename T::Scalar &s)
    {
        value /= s;
        return *this;
    }

    friend Scalar operator+(Scalar a, const Scalar &b) { return a += b; }
    friend Scalar operator-(Scalar a, const Scalar &b) { return a -= b; }
    friend Scalar operator*(Scalar a, const typename T::Scalar &s) { return a *= s; }
    friend Scalar operator*(const typename T::Scalar &s, Scalar a) { return a *= s; }
    friend Scalar operator/(Scalar a, const typename T::Scalar &s) { return a /= s; }

    // Comparison
    bool operator==(const Scalar &rhs) const { return value == rhs.value; }

    // Output
    friend std::ostream &operator<<(std::ostream &os, const Scalar &s)
    {

        return os << s.value;
    }
};

template <typename T>
struct ZeroData<Scalar<T>>
{
    static inline const Scalar<T> value = Scalar<T>(ZeroData<T>::value);
};

template <typename T, int N>
Scalar<T> dotProduct(const Eigen::Matrix<Real, N, 1> &real_vector, const Eigen::Matrix<Scalar<T>, N, 1> &scalar_vector)
{
    Scalar<T> scalar = Scalar<T>(0);
    for (UnsignedInt i = 0; i < N; ++i)
    {
        scalar += real_vector[i] * scalar_vector[i];
    }
    return scalar;
};

template <typename T, int N, int M>
Eigen::Matrix<Scalar<T>, N, M> scalarProduct(const Eigen::Matrix<Real, N, M> &real_matrix, const Scalar<T> &scalar)
{
    Eigen::Matrix<Scalar<T>, N, M> scalar_matrix = Eigen::Matrix<Scalar<T>, N, M>::Zero();
    for (UnsignedInt i = 0; i < N; ++i)
        for (UnsignedInt j = 0; j < M; ++j)
        {
            scalar_matrix(i, j) = real_matrix(i, j) * scalar;
        }

    return scalar_matrix;
};

} // namespace SPH

namespace Eigen
{
template <typename T>
struct NumTraits<SPH::Scalar<T>> : GenericNumTraits<SPH::Scalar<T>>
{
    using Real = SPH::Real;
    using NonInteger = SPH::Real;
    using Nested = SPH::Real;

    enum
    {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 1,
        MulCost = 1
    };
};
} // namespace Eigen
#endif // SCALAR_NUMERICS_H
