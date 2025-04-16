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
    T value_;

    // Constructors
    Scalar() = default;
    Scalar(const T &v) : value_(v) {}
    Scalar(int) : value_(initZero()) {}

    static T initZero()
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            return 0.0;
        }
        else
        {
            return T::Zero(); // Assumes T has a static Zero() method
        }
    }
    // Implicit conversion
    operator const T &() const { return value_; }
    operator T &() { return value_; }

    // Access
    T &ref() { return value_; }
    const T &get() const { return value_; }

    // Arithmetic (only if T supports it)
    Scalar &operator+=(const Scalar &rhs)
    {
        value_ += rhs.value_;
        return *this;
    }
    Scalar &operator-=(const Scalar &rhs)
    {
        value_ -= rhs.value_;
        return *this;
    }
    template <typename ArithmeticType>
    Scalar &operator*=(const ArithmeticType &s)
    {
        value_ *= s;
        return *this;
    }
    template <typename ArithmeticType>
    Scalar &operator/=(const ArithmeticType &s)
    {
        value_ /= s;
        return *this;
    }

    friend Scalar operator+(Scalar a, const Scalar &b) { return a += b; }
    friend Scalar operator-(Scalar a, const Scalar &b) { return a -= b; }
    template <typename ArithmeticType>
    friend Scalar operator*(Scalar a, const ArithmeticType &s) { return a *= s; }
    template <typename ArithmeticType>
    friend Scalar operator*(const ArithmeticType &s, Scalar a) { return a *= s; }
    template <typename ArithmeticType>
    friend Scalar operator/(Scalar a, const ArithmeticType &s) { return a /= s; }

    // Comparison
    bool operator==(const Scalar &rhs) const { return value_ == rhs.value_; }

    // Output
    friend std::ostream &operator<<(std::ostream &os, const Scalar &s)
    {

        return os << s.value_;
    }
};

template <typename T>
struct ZeroData<Scalar<T>>
{
    static inline const Scalar<T> value = Scalar<T>(ZeroData<T>::value_);
};

template <typename T, int N>
using ScalarVec = Eigen::Matrix<Scalar<T>, N, 1>; // vector of generalized scalar

template <typename T, int N>
struct ZeroData<ScalarVec<T, N>>
{
    static inline const ScalarVec<T, N> value = ScalarVec<T, N>::Zero();
};

template <int N, int M>
using ScalarVecVec = ScalarVec<Eigen::Matrix<Real, M, 1>, N>; // vector of generalized scalar of vector

template <int N, int M>
ScalarVecVec<N, M> MatrixToScalarVecVec(const Eigen::Matrix<Real, N, M> &input)
{
    ScalarVecVec<N, M> result = ScalarVecVec<N, M>::Zero();
    for (UnsignedInt i = 0; i < N; ++i)
    {
        result[i] = Scalar<Eigen::Matrix<Real, M, 1>>(input.row(i).transpose());
    }
    return result;
};

template <int N>
ScalarVec<Real, N> MatrixToScalarVecVec(const Eigen::Matrix<Real, N, 1> &input)
{
    ScalarVec<Real, N> result = ScalarVec<Real, N>::Zero();
    for (UnsignedInt i = 0; i < N; ++i)
    {
        result[i] = Scalar<Real>(input[i]);
    }
    return result;
};

template <int N, int M>
Eigen::Matrix<Real, N, M> ScalarVecVecToMatrix(const ScalarVecVec<N, M> &input)
{
    Eigen::Matrix<Real, N, M> result = Eigen::Matrix<Real, N, M>::Zero();
    for (UnsignedInt i = 0; i < N; ++i)
    {
        result.row(i) = input[i].get().transpose();
    }
    return result;
};

template <int N>
Eigen::Matrix<Real, N, 1> ScalarVecVecToMatrix(const ScalarVec<Real, N> &input)
{
    Eigen::Matrix<Real, N, 1> result = Eigen::Matrix<Real, N, 1>::Zero();
    for (UnsignedInt i = 0; i < N; ++i)
    {
        result[i] = input[i].get();
    }
    return result;
};

template <typename T, int N>
Scalar<T> dotProduct(const Eigen::Matrix<Real, N, 1> &real_vector, const ScalarVec<T, N> &scalar_vector)
{
    Scalar<T> scalar = Scalar<T>(0);
    for (UnsignedInt i = 0; i < N; ++i)
    {
        scalar += real_vector[i] * scalar_vector[i];
    }
    return scalar;
};

template <typename T, int N, int M>
using ScalarMat = Eigen::Matrix<Scalar<T>, N, M>; // matrix of generalized scalar

template <typename T, int N, int M>
ScalarMat<T, N, M> scalarProduct(const Eigen::Matrix<Real, N, M> &real_matrix, const Scalar<T> &scalar)
{
    ScalarMat<T, N, M> scalar_matrix = ScalarMat<T, N, M>::Zero();
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
