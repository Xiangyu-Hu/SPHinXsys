// MIT License
//
// Copyright (c) 2017 Martin Bisson
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#ifndef __POLAR_DECOMPOSITION_3X3_MATRIX_H__
#define __POLAR_DECOMPOSITION_3X3_MATRIX_H__

#include <cassert>
#include <cmath>
#include <limits>

namespace polar
{

namespace detail
{

// Helper methods for standard mathematical function calls.
template <typename TReal>
struct math_utils
{
    static inline TReal sqrt(const TReal x);
    static inline TReal fabs(const TReal x);
    static inline TReal clamp(const TReal x, const TReal a, const TReal b);
    static inline TReal cos(const TReal x);
    static inline TReal acos(const TReal x);
    static inline TReal ceil(const TReal x);
    static inline TReal log10(const TReal x);
};

template <>
inline float math_utils<float>::sqrt(const float x)
{
    return ::sqrtf(x);
}
template <>
inline double math_utils<double>::sqrt(const double x)
{
    return ::sqrt(x);
}

template <>
inline float math_utils<float>::fabs(const float x)
{
    return ::fabs(x);
}
template <>
inline double math_utils<double>::fabs(const double x)
{
    return ::fabs(x);
}

template <typename TReal>
inline TReal math_utils<TReal>::clamp(const TReal x, const TReal a, const TReal b)
{
    assert(a <= b);
    if (x < a)
        return a;
    else if (x > b)
        return b;
    else
        return x;
}

template <>
inline float math_utils<float>::cos(const float x)
{
    return ::cosf(x);
}
template <>
inline double math_utils<double>::cos(const double x)
{
    return ::cos(x);
}

template <>
inline float math_utils<float>::acos(const float x)
{
    return ::acosf(x);
}
template <>
inline double math_utils<double>::acos(const double x)
{
    return ::acos(x);
}

template <>
inline float math_utils<float>::ceil(const float x)
{
    return ::ceilf(x);
}
template <>
inline double math_utils<double>::ceil(const double x)
{
    return ::ceil(x);
}

template <>
inline float math_utils<float>::log10(const float x)
{
    return ::log10f(x);
}
template <>
inline double math_utils<double>::log10(const double x)
{
    return ::log10(x);
}

}; // End of namespace detail.

// Implementation of the basic matrix type used in the algorithm.
namespace detail
{

// Definition of the column-major matrix type.
template <typename TReal, int NbColumns, int NbRows>
class matrix
{
  public:
    inline TReal operator()(int columnIndex, int rowIndex) const;
    inline TReal &operator()(int columnIndex, int rowIndex);
    inline TReal operator()(int index) const;
    inline TReal &operator()(int index);

  private:
    TReal data[NbColumns * NbRows];
};

// We rely on the compiler to be efficient in inlining the calls to this
// method, and especially optimize the double-indices access and convert them
// to single access ones.
template <typename TReal, int NbRows, int NbColumns>
inline TReal matrix<TReal, NbRows, NbColumns>::operator()(int columnIndex, int rowIndex) const
{
    assert(0 <= rowIndex);
    assert(rowIndex < NbRows);
    assert(0 <= columnIndex);
    assert(columnIndex < NbColumns);
    return data[columnIndex * NbRows + rowIndex];
}

template <typename TReal, int NbRows, int NbColumns>
inline TReal &matrix<TReal, NbRows, NbColumns>::operator()(int columnIndex, int rowIndex)
{
    assert(0 <= rowIndex);
    assert(rowIndex < NbRows);
    assert(0 <= columnIndex);
    assert(columnIndex < NbColumns);
    return data[columnIndex * NbRows + rowIndex];
}

template <typename TReal, int NbRows, int NbColumns>
inline TReal matrix<TReal, NbRows, NbColumns>::operator()(int index) const
{
    assert(0 <= index);
    assert(index < NbRows * NbColumns);
    return data[index];
}

template <typename TReal, int NbRows, int NbColumns>
inline TReal &matrix<TReal, NbRows, NbColumns>::operator()(int index)
{
    assert(0 <= index);
    assert(index < NbRows * NbColumns);
    return data[index];
}

// Definition of the vector type.
template <typename TReal, int NbElements>
class vector
{
  public:
    inline TReal operator()(int index) const;
    inline TReal &operator()(int index);

  private:
    TReal data[NbElements];
};

template <typename TReal, int NbElements>
inline TReal vector<TReal, NbElements>::operator()(int index) const
{
    assert(0 <= index);
    assert(index < NbElements);
    return data[index];
}

template <typename TReal, int NbElements>
inline TReal &vector<TReal, NbElements>::operator()(int index)
{
    assert(0 <= index);
    assert(index < NbElements);
    return data[index];
}

//
// The following re-implements matrix operations that are very common in
// most linear algebra packages such as GLM or Eigen.  However, by design
// this library does not want to depend on such packages, so the simple
// matrix operations are re-implemented here.
//
// Most methods are not implemented in the most generic way; they are
// specialized for the ways they are used in this algorithm specific
// implementation.
//
// They definitely could be implemented in a "cleaner", shorter way
// using the proper libraries.
//

template <typename TReal>
inline void normalize(matrix<TReal, 3, 3> &m)
{
    TReal length = m(0) * m(0);
    length += m(1) * m(1);
    length += m(2) * m(2);
    length += m(3) * m(3);
    length += m(4) * m(4);
    length += m(5) * m(5);
    length += m(6) * m(6);
    length += m(7) * m(7);
    length += m(8) * m(8);

    // Add small numerical value to make sure we are ok with 0-length.
    const TReal factor = 1 / (math_utils<TReal>::sqrt(length) + (std::numeric_limits<TReal>::min)());

    m(0) *= factor;
    m(1) *= factor;
    m(2) *= factor;
    m(3) *= factor;
    m(4) *= factor;
    m(5) *= factor;
    m(6) *= factor;
    m(7) *= factor;
    m(8) *= factor;
}

template <typename TReal>
inline void normalize(vector<TReal, 4> &m)
{
    TReal length = m(0) * m(0);
    length += m(1) * m(1);
    length += m(2) * m(2);
    length += m(3) * m(3);

    // Add small numerical value to make sure we are ok with 0-length.
    const TReal factor = 1 / (math_utils<TReal>::sqrt(length) + (std::numeric_limits<TReal>::min)());

    m(0) *= factor;
    m(1) *= factor;
    m(2) *= factor;
    m(3) *= factor;
}

template <typename TReal>
inline void transpose_multiply(
    matrix<TReal, 3, 3> &result,
    const matrix<TReal, 3, 3> &a,
    const matrix<TReal, 3, 3> &b)
{
    // Perform transposition at the same time as the multiplication.
    result(0, 0) = a(0, 0) * b(0, 0) + a(0, 1) * b(0, 1) + a(0, 2) * b(0, 2);
    result(0, 1) = a(1, 0) * b(0, 0) + a(1, 1) * b(0, 1) + a(1, 2) * b(0, 2);
    result(0, 2) = a(2, 0) * b(0, 0) + a(2, 1) * b(0, 1) + a(2, 2) * b(0, 2);
    result(1, 0) = a(0, 0) * b(1, 0) + a(0, 1) * b(1, 1) + a(0, 2) * b(1, 2);
    result(1, 1) = a(1, 0) * b(1, 0) + a(1, 1) * b(1, 1) + a(1, 2) * b(1, 2);
    result(1, 2) = a(2, 0) * b(1, 0) + a(2, 1) * b(1, 1) + a(2, 2) * b(1, 2);
    result(2, 0) = a(0, 0) * b(2, 0) + a(0, 1) * b(2, 1) + a(0, 2) * b(2, 2);
    result(2, 1) = a(1, 0) * b(2, 0) + a(1, 1) * b(2, 1) + a(1, 2) * b(2, 2);
    result(2, 2) = a(2, 0) * b(2, 0) + a(2, 1) * b(2, 1) + a(2, 2) * b(2, 2);
}

template <typename TReal>
inline void multiply(
    matrix<TReal, 4, 4> &result,
    const TReal factor)
{
    result(0) *= factor;
    result(1) *= factor;
    result(2) *= factor;
    result(3) *= factor;
    result(4) *= factor;
    result(5) *= factor;
    result(6) *= factor;
    result(7) *= factor;
    result(8) *= factor;
    result(9) *= factor;
    result(10) *= factor;
    result(11) *= factor;
    result(12) *= factor;
    result(13) *= factor;
    result(14) *= factor;
    result(15) *= factor;
}

template <typename TReal>
inline void multiply(
    matrix<TReal, 2, 2> &result,
    const TReal factor)
{
    result(0) *= factor;
    result(1) *= factor;
    result(2) *= factor;
    result(3) *= factor;
}

template <typename TReal>
inline TReal dot(const vector<TReal, 4> &a, const vector<TReal, 4> &b)
{
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2) + a(3) * b(3);
}

}; // End of namespace detail.

}; // End of namespace polar.

#endif // __POLAR_DECOMPOSITION_3X3_MATRIX_H__
