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
#ifndef __POLAR_DECOMPOSITION_3X3_IMPL_H__
#define __POLAR_DECOMPOSITION_3X3_IMPL_H__

#include "polar_decomposition_3x3_matrix.h"

#include <algorithm>
#include <limits>

namespace polar
{

// This part implements portions of the algorithms that are tailored
// to just the specifics of what is needed.
namespace detail
{

// This method is based on the Matlab implementation.
// It computes the determinant of B matrix from an LU factorization with partial pivoting,
// except that it actually computes it from the values of A, (B matrix is a function of A).
template <typename TReal>
inline TReal compute_b_determinant_from_a_matrix_lu_partial(const matrix<TReal, 3, 3> &A)
{
    TReal b = 0;
    {
        TReal temp;
        temp = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2));
        b += (temp * temp);
        temp = (A(0, 1) * A(2, 2) - A(2, 1) * A(0, 2));
        b += (temp * temp);
        temp = (A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2));
        b += (temp * temp);
        temp = (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2));
        b += (temp * temp);
        temp = (A(0, 0) * A(2, 2) - A(2, 0) * A(0, 2));
        b += (temp * temp);
        temp = (A(1, 0) * A(2, 2) - A(2, 0) * A(1, 2));
        b += (temp * temp);
        temp = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1));
        b += (temp * temp);
        temp = (A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1));
        b += (temp * temp);
        temp = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1));
        b += (temp * temp);
        b = -4 * b + 1;
    }
    return b;
}

// This method is based on the Matlab implementation.
// It computes the determinant of A matrix from an LU factorization with partial pivoting.
template <typename TReal>
inline TReal compute_determinant_lu_partial(const matrix<TReal, 3, 3> &A, TReal &d)
{
    TReal dd; // Determinant, d is the sign of the determinant.
    // ANSME: Do we really need to keep track of the sign, or can't we just check it at the end?

    matrix<TReal, 3, 3> AA;
    const TReal absA00 = math_utils<TReal>::fabs(A(0, 0));
    const TReal absA01 = math_utils<TReal>::fabs(A(0, 1));
    const TReal absA02 = math_utils<TReal>::fabs(A(0, 2));
    if (absA01 > absA02)
    {
        if (absA00 > absA01)
        {
            AA = A;
            dd = 1;
        }
        else
        {
            AA(0, 0) = A(0, 1);
            AA(1, 0) = A(1, 1);
            AA(2, 0) = A(2, 1);
            AA(0, 1) = A(0, 0);
            AA(1, 1) = A(1, 0);
            AA(2, 1) = A(2, 0);
            AA(0, 2) = A(0, 2);
            AA(1, 2) = A(1, 2);
            AA(2, 2) = A(2, 2);
            dd = -1;
        }
    }
    else
    {
        if (absA00 > absA02)
        {
            AA = A;
            dd = 1;
        }
        else
        {
            AA(0, 0) = A(0, 2);
            AA(1, 0) = A(1, 2);
            AA(2, 0) = A(2, 2);
            AA(0, 1) = A(0, 1);
            AA(1, 1) = A(1, 1);
            AA(2, 1) = A(2, 1);
            AA(0, 2) = A(0, 0);
            AA(1, 2) = A(1, 0);
            AA(2, 2) = A(2, 0);
            dd = -1;
        }
    }

    d = dd;
    vector<TReal, 3> U;
    U(0) = AA(0, 0);
    if (U(0) < 0)
        d = -d;

    const TReal m1 = AA(1, 0) / AA(0, 0);
    const TReal m2 = AA(2, 0) / AA(0, 0);
    const TReal AA00 = AA(1, 1) - AA(0, 1) * m1;
    const TReal AA10 = AA(2, 1) - AA(0, 1) * m2;
    const TReal AA01 = AA(1, 2) - AA(0, 2) * m1;
    const TReal AA11 = AA(2, 2) - AA(0, 2) * m2;

    if (math_utils<TReal>::fabs(AA00) < math_utils<TReal>::fabs(AA01))
    {
        U(1) = AA01;
        U(2) = AA10 - AA00 * AA11 / AA01;
        dd = -dd;
        d = -d;
        if (U(1) < 0)
            d = -d;
        if (U(2) < 0)
            d = -d;
    }
    else if (AA00 == 0)
    {
        U(1) = 0;
        U(2) = 0;
    }
    else
    {
        U(1) = AA00;
        U(2) = AA11 - AA01 * AA10 / AA00;
        if (U(1) < 0)
            d = -d;
        if (U(2) < 0)
            d = -d;
    }

    dd = dd * U(0) * U(1) * U(2);
    if (d == 0)
        d = 1;

    assert(((d < 0) && (dd < 0)) || (dd >= 0));
    assert((d == 1) || (d == -1));

    return dd;
}

// Swap rows.
template <typename TReal>
inline void swap_rows(matrix<TReal, 3, 3> &m, const int row0, const int row1)
{
    assert(row0 != row1);
    std::swap(m(0, row0), m(0, row1));
    std::swap(m(1, row0), m(1, row1));
    std::swap(m(2, row0), m(2, row1));
}
template <typename TReal>
inline void swap_rows(matrix<TReal, 4, 4> &m, const int row0, const int row1)
{
    assert(row0 != row1);
    std::swap(m(0, row0), m(0, row1));
    std::swap(m(1, row0), m(1, row1));
    std::swap(m(2, row0), m(2, row1));
    std::swap(m(3, row0), m(3, row1));
}

// Swap columns.
template <typename TReal>
inline void swap_columns(matrix<TReal, 3, 3> &m, const int column0, const int column1)
{
    assert(column0 != column1);
    std::swap(m(column0, 0), m(column1, 0));
    std::swap(m(column0, 1), m(column1, 1));
    std::swap(m(column0, 2), m(column1, 2));
}
template <typename TReal>
inline void swap_columns(matrix<TReal, 4, 4> &m, const int column0, const int column1)
{
    assert(column0 != column1);
    std::swap(m(column0, 0), m(column1, 0));
    std::swap(m(column0, 1), m(column1, 1));
    std::swap(m(column0, 2), m(column1, 2));
    std::swap(m(column0, 3), m(column1, 3));
}

// This method is based on the Matlab implementation.
// It computes an LDL^T factorization with diagonal pivoting, P^T Bs P = L D L^T.
// This method modifies the Bs matrix.
// Note: The method does not fill the whole L matrix, just the lower left part.
// The caller should assume:
//       - L(i,i) == 1
//       - L(i,j) == 0 for i > j
template <typename TReal>
inline void compute_ldlt_factorization_diagonal(matrix<TReal, 4, 4> &L, vector<TReal, 4> &D, vector<int, 4> &p, matrix<TReal, 4, 4> &Bs)
{
    p(0) = 0;
    p(1) = 1;
    p(2) = 2;
    p(3) = 3;

    // ANSME: Should we compare absolute values to pick which row to pivot?
    //        This whole code could be refactored and improved to be more coherent
    //        from one step to another.

    // First step.
    {
        int r = 3;
        if (Bs(3, 3) < Bs(2, 2))
            r = 2;
        if (Bs(r, r) < Bs(1, 1))
            r = 1;

        if (Bs(r, r) > Bs(0, 0))
        {
            std::swap(p(0), p(r));
            swap_rows(Bs, 0, r);
            swap_columns(Bs, 0, r);
        }

        D(0) = Bs(0, 0);
        L(0, 1) = Bs(0, 1) / D(0);
        L(0, 2) = Bs(0, 2) / D(0);
        L(0, 3) = Bs(0, 3) / D(0);

        Bs(1, 1) = Bs(1, 1) - L(0, 1) * Bs(1, 0);
        Bs(1, 2) = Bs(1, 2) - L(0, 1) * Bs(2, 0);
        Bs(2, 1) = Bs(1, 2);
        Bs(1, 3) = Bs(1, 3) - L(0, 1) * Bs(3, 0);
        Bs(3, 1) = Bs(1, 3);

        Bs(2, 2) = Bs(2, 2) - L(0, 2) * Bs(2, 0);
        Bs(2, 3) = Bs(2, 3) - L(0, 2) * Bs(3, 0);
        Bs(3, 2) = Bs(2, 3);

        Bs(3, 3) = Bs(3, 3) - L(0, 3) * Bs(3, 0);
    }

    // Second step.
    {
        int r = 3;
        if (Bs(3, 3) < Bs(2, 2))
            r = 2;

        if (Bs(r, r) > Bs(1, 1))
        {
            std::swap(p(1), p(r));
            swap_rows(Bs, 1, r);
            swap_columns(Bs, 1, r);
#if 0
                    swap_rows(L, 1, r);
                    swap_columns(L, 1, r);
#else
            // Here, only the first column has been written, so we can swap just that.
            //                 swap(1,2)      swap(1,3)
            // | 1 0 0 0 |    | 1 0 0 0 |    | 1 0 0 0 |
            // | a 1 0 0 |    | b 1 0 0 |    | c 1 0 0 |
            // | b 0 1 0 |    | a 0 1 0 |    | b 0 1 0 |
            // | c 0 0 1 |    | c 0 0 1 |    | a 0 0 1 |
            std::swap(L(0, 1), L(0, r));
#endif
        }

        D(1) = Bs(1, 1);
        L(1, 2) = Bs(1, 2) / D(1);
        L(1, 3) = Bs(1, 3) / D(1);

        Bs(2, 2) = Bs(2, 2) - L(1, 2) * Bs(2, 1);
        Bs(2, 3) = Bs(2, 3) - L(1, 2) * Bs(3, 1);
        Bs(3, 2) = Bs(2, 3);

        Bs(3, 3) = Bs(3, 3) - L(1, 3) * Bs(3, 1);
    }

    // Third step.
    {
        if (Bs(2, 2) < Bs(3, 3))
        {
            D(2) = Bs(3, 3);
            std::swap(p(2), p(3));
            swap_rows(Bs, 2, 3);
            swap_columns(Bs, 2, 3);
#if 0
                    swap_rows(L, 2, 3);
                    swap_columns(L, 2, 3);
#else
            // Here, only the first two columns has been written, so we can swap just that.
            //                 swap(2,3)
            // | 1 0 0 0 |    | 1 0 0 0 |
            // | a 1 0 0 |    | a 1 0 0 |
            // | b e 1 0 |    | c f 1 0 |
            // | c f 0 1 |    | b e 0 1 |
            std::swap(L(0, 2), L(0, 3));
            std::swap(L(1, 2), L(1, 3));
#endif
        }
        else
        {
            D(2) = Bs(2, 2);
        }

        L(2, 3) = Bs(2, 3) / D(2);
    }
}

// This method is based on the Matlab implementation.
// It computes the determinant of A matrix from an LU factorization with complete pivoting.
template <typename TReal>
inline TReal compute_determinant_lu_complete(const matrix<TReal, 3, 3> &A, TReal &d, TReal &u22)
{
    TReal dd; // Determinant, d is the sign of the determinant.
    // ANSME: Do we really need to keep track of the sign, or can't we just check it at the end?

    matrix<TReal, 3, 3> AA = A;
    {
        int r = 0;
        int c = 0;
        dd = 1;
        if (math_utils<TReal>::fabs(A(0, 1)) > math_utils<TReal>::fabs(A(0, 0)))
        {
            r = 1;
        }
        if (math_utils<TReal>::fabs(A(0, 2)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 2;
        }
        if (math_utils<TReal>::fabs(A(1, 0)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 0;
            c = 1;
        }
        if (math_utils<TReal>::fabs(A(1, 1)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 1;
            c = 1;
        }
        if (math_utils<TReal>::fabs(A(1, 2)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 2;
            c = 1;
        }
        if (math_utils<TReal>::fabs(A(2, 0)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 0;
            c = 2;
        }
        if (math_utils<TReal>::fabs(A(2, 1)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 1;
            c = 2;
        }
        if (math_utils<TReal>::fabs(A(2, 2)) > math_utils<TReal>::fabs(A(c, r)))
        {
            r = 2;
            c = 2;
        }

        if (r > 0)
        {
            swap_rows(AA, 0, r);
            dd = -1;
        }
        if (c > 0)
        {
            swap_columns(AA, 0, c);
            dd = -dd;
        }
    }

    vector<TReal, 3> U;
    U(0) = AA(0, 0);

    matrix<TReal, 2, 2> AA_;
    // In case the whole matrix is 0, we add a small value.
    const TReal m1 = AA(1, 0) / (AA(0, 0) + (std::numeric_limits<TReal>::min)());
    const TReal m2 = AA(2, 0) / (AA(0, 0) + (std::numeric_limits<TReal>::min)());
    AA_(0, 0) = AA(1, 1) - AA(0, 1) * m1;
    AA_(1, 0) = AA(2, 1) - AA(0, 1) * m2;
    AA_(0, 1) = AA(1, 2) - AA(0, 2) * m1;
    AA_(1, 1) = AA(2, 2) - AA(0, 2) * m2;

    {
        int r = 0;
        int c = 0;
        if (math_utils<TReal>::fabs(AA_(0, 1)) > math_utils<TReal>::fabs(AA_(0, 0)))
        {
            r = 1;
        }
        if (math_utils<TReal>::fabs(AA_(1, 0)) > math_utils<TReal>::fabs(AA_(c, r)))
        {
            r = 0;
            c = 1;
        }
        if (math_utils<TReal>::fabs(AA_(1, 1)) > math_utils<TReal>::fabs(AA_(c, r)))
        {
            r = 1;
            c = 1;
        }

        if (r > 0)
        {
            dd = -dd;
        }
        if (c > 0)
        {
            dd = -dd;
        }
        U(1) = AA_(c, r);
        if (U(1) == 0)
        {
            U(2) = 0;
        }
        else
        {
            U(2) = AA_(1 - c, 1 - r) - AA_(1 - c, r) * AA_(c, 1 - r) / U(1);
        }
    }

    d = dd;
    dd = dd * U(0) * U(1) * U(2);
    if (U(0) < 0)
        d = -d;
    if (U(1) < 0)
        d = -d;
    if (U(2) < 0)
        d = -d;

    u22 = U(1);

    assert(((d < 0) && (dd < 0)) || (dd >= 0));
    assert((d == 1) || (d == -1));

    return dd;
}

// This method is based on the Matlab implementation.
// It computes an LDL^T factorization by block LDL^T factorization with Bunch-Parlett pivoting.
// This method modifies the Bs matrix.
// Note: The method does not fill the whole L matrix, just the lower left part.
// The caller should assume:
//       - L(i,i) == 1
//       - L(i,j) == 0 for i > j
//       - L(2,3) == 0
// FIXME: Not sure we actually need a 4x4 matrix for D.
template <typename TReal>
inline void compute_ldlt_factorization_bunch_parlett(matrix<TReal, 4, 4> &L, matrix<TReal, 4, 4> &D, vector<int, 4> &p, matrix<TReal, 4, 4> &Bs)
{
    p(0) = 0;
    p(1) = 1;
    p(2) = 2;
    p(3) = 3;

    // ANSME: Should we compare absolute values to pick which row to pivot?
    //        This whole code could be refactored and improved to be more coherent
    //        from one step to another.

    // First step.
    {
        int r = 3;
        if (Bs(3, 3) < Bs(2, 2))
            r = 2;
        if (Bs(r, r) < Bs(1, 1))
            r = 1;

        if (Bs(r, r) > Bs(0, 0))
        {
            std::swap(p(0), p(r));
            swap_rows(Bs, 0, r);
            swap_columns(Bs, 0, r);
        }

        D(0, 0) = Bs(0, 0);
        L(0, 1) = Bs(0, 1) / D(0, 0);
        L(0, 2) = Bs(0, 2) / D(0, 0);
        L(0, 3) = Bs(0, 3) / D(0, 0);

        Bs(1, 1) = Bs(1, 1) - L(0, 1) * Bs(1, 0);
        Bs(1, 2) = Bs(1, 2) - L(0, 1) * Bs(2, 0);
        Bs(2, 1) = Bs(1, 2);
        Bs(1, 3) = Bs(1, 3) - L(0, 1) * Bs(3, 0);
        Bs(3, 1) = Bs(1, 3);

        Bs(2, 2) = Bs(2, 2) - L(0, 2) * Bs(2, 0);
        Bs(2, 3) = Bs(2, 3) - L(0, 2) * Bs(3, 0);
        Bs(3, 2) = Bs(2, 3);

        Bs(3, 3) = Bs(3, 3) - L(0, 3) * Bs(3, 0);
    }

    // Second step.
    {
        int r = 2;
        if (Bs(2, 2) < Bs(1, 1))
            r = 1;

        if (Bs(r, r) > Bs(0, 0))
        {
            std::swap(p(1), p(r));
            swap_rows(Bs, 1, r);
            swap_columns(Bs, 1, r);
#if 0
                    swap_rows(L, 1, r);
                    swap_columns(L, 1, r);
#else
            // Here, only the first column has been written, so we can swap just that.
            //                 swap(1,2)      swap(1,3)
            // | 1 0 0 0 |    | 1 0 0 0 |    | 1 0 0 0 |
            // | a 1 0 0 |    | b 1 0 0 |    | c 1 0 0 |
            // | b 0 1 0 |    | a 0 1 0 |    | b 0 1 0 |
            // | c 0 0 1 |    | c 0 0 1 |    | a 0 0 1 |
            std::swap(L(0, 1), L(0, r));
#endif
        }

        D(1, 1) = Bs(1, 1);
        L(1, 2) = Bs(1, 2) / D(1, 1);
        L(1, 3) = Bs(1, 3) / D(1, 1);

        D(2, 2) = Bs(2, 2) - L(1, 2) * Bs(2, 1);
        D(2, 3) = Bs(2, 3) - L(1, 2) * Bs(3, 1);
        D(3, 2) = D(2, 3);

        D(3, 3) = Bs(3, 3) - L(1, 3) * Bs(3, 1);
    }
}

// This computes the null space of the matrix, knowing that it is symmetric:
//
// | a    b |
// | b    c |
template <typename TReal>
inline vector<TReal, 2> compute_null_space(const TReal a, const TReal b, const TReal c)
{
    // We assume that the matrix determinant is 0.
    assert(0 == a * c - b * b);

    vector<TReal, 2> nullSpace;
    if (a != 0)
    {
        // Transform: R1 = R1 / a
        //
        // | 1    b / a |
        // | b    c     |
        //
        // R2 = R2 - b R1
        //
        // | 1    b / a         |
        // | 0    c - b * b / a |
        //
        // Because a != 0 and a * c - b * b == 0, we have
        // 1 * x1 + b / a * x2 = 0
        // a x1 + b x2 = 0
        //
        // solution: [ b    -a ]^T
        nullSpace(0) = b;
        nullSpace(1) = -a;
    }
    else
    {
        // Since a == 0 and the determinant is null we have:
        //
        // 0 * c - b * b = 0
        // b = 0
        //
        // If b == 0
        //
        // We have:
        //
        // | 0    0 |
        // | b    c |
        //
        // b x1 + c x2 = 0
        // solution: [ c    -b ]^T
        //
        // (With b == 0, [ c    0 ]^T).
        //
        assert(a == 0);
        assert(b == 0);
        assert(c != 0);
        nullSpace(0) = c;
        nullSpace(1) = -b;
    }

    return nullSpace;
}

// These methods are used for optimized operations in reverse iteration with LDL^T.
template <typename TReal>
inline vector<TReal, 4> multiply_il_v(
    const TReal IL01, const TReal IL02, const TReal IL03, const TReal IL12, const TReal IL13,
    const vector<TReal, 4> &v)
{
    // {{1,0,0,0},{m_12,1,0,0},{m_13,m_23,1,0},{m_14,m_24,0,1}} * {{v_1},{v_2},{v_3},{v_4}}
    //
    // | 1       0       0       0 |       | v1 |       | v1                         |
    // | IL01    1       0       0 |   *   | v2 |   =   | v1 * IL01 + v              |
    // | IL02    IL12    1       0 |       | v3 |       | v1 * IL02 + v2 * IL12 + v3 |
    // | IL03    IL13    0       1 |       | v4 |       | v1 * IL03 + v2 * IL13 + v4 |
    vector<TReal, 4> result;
    result(0) = v(0);
    result(1) = v(0) * IL01 + v(1);
    result(2) = v(0) * IL02 + v(1) * IL12 + v(2);
    result(3) = v(0) * IL03 + v(1) * IL13 + v(3);
    return result;
}

template <typename TReal>
inline vector<TReal, 4> multiply_id_v(
    const TReal ID00, const TReal ID11, const matrix<TReal, 2, 2> &ID,
    const vector<TReal, 4> &v)
{
    // {{a_11,0,0,0},{0,a_22,0,0},{0,0,b_11,b_21},{0,0,b_12,b_22}} * {{v_1},{v_2},{v_3},{v_4}}
    //
    // | ID00    0       0          0       |       | v1 |       | v1 * ID00                   |
    // | 0000    ID11    0          0       |   *   | v2 |   =   | v2 * ID11                   |
    // | 0       0       ID(0,0)    ID(1,0) |       | v3 |       | v3 * ID(0,0) + v4 * ID(1,0) |
    // | 0       0       ID(0,1)    ID(1,1) |       | v4 |       | v3 * ID(0,1) + v4 * ID(1,1) |
    vector<TReal, 4> result;
    result(0) = v(0) * ID00;
    result(1) = v(1) * ID11;
    result(2) = v(2) * ID(0, 0) + v(3) * ID(1, 0);
    result(3) = v(2) * ID(0, 1) + v(3) * ID(1, 1);
    return result;
}

template <typename TReal>
inline vector<TReal, 4> multiply_v_il(
    const vector<TReal, 4> &v,
    const TReal IL01, const TReal IL02, const TReal IL03, const TReal IL12, const TReal IL13)
{
    // Transpose[{{1,0,0,0},{m_12,1,0,0},{m_13,m_23,1,0},{m_14,m_24,0,1}}] * {{v_1},{v_2},{v_3},{v_4}}
    // {v_1,v_2,v_3,v_4} * {{1,0,0,0},{m_12,1,0,0},{m_13,m_23,1,0},{m_14,m_24,0,1}}
    //
    // | v1 |^T   *   | 1       0       0       0 |       | v1 + v2 * IL01 + v3 * IL02 + v4 * IL03 |^T
    // | v2 |         | IL01    1       0       0 |   =   | v2 + v3 * IL12 + v4 * IL13             |
    // | v3 |         | IL02    IL12    1       0 |       | v3                                     |
    // | v4 |         | IL03    IL13    0       1 |       | v4                                     |
    vector<TReal, 4> result;
    result(0) = v(0) + v(1) * IL01 + v(2) * IL02 + v(3) * IL03;
    result(1) = v(1) + v(2) * IL12 + v(3) * IL13;
    result(2) = v(2);
    result(3) = v(3);
    return result;
}

template <typename TReal>
inline vector<TReal, 4> multiply_minus_v_d(
    const vector<TReal, 4> &v,
    const matrix<TReal, 4, 4> &D)
{
    // -{v1,v2,v3,v4} * {{D_11,0,0,0},{0,D_22,0,0},{0,0,D_33,D_43},{0,0,D_34,D_44}}
    //
    // - | v1 |^T   *   | D(0,0)    0         0            0      |       | -v1 * D(0,0)               |^T
    //   | v2 |         | 0         D(1,1)    0            0      |   =   | -v2 * D(1,1)               |
    //   | v3 |         | 0         0         D(2,2)       D(3,2) |       | -v3 * D(2,2) - v4 * D(2,3) |
    //   | v4 |         | 0         0         D(2,3)       D(3,3) |       | -v4 * D(3,2) - v4 * D(3,3) |
    vector<TReal, 4> result;
    result(0) = -v(0) * D(0, 0);
    result(1) = -v(1) * D(1, 1);
    result(2) = -v(2) * D(2, 2) - v(3) * D(2, 3);
    result(3) = -v(2) * D(3, 2) - v(3) * D(3, 3);
    return result;
}

// Orthonormalization of matrices of the type:
//
// | v00    v10 |
// | v01    v11 |
// | 1      0   |
// | 0      1   |
template <typename TReal>
inline void orthonormalize_v_with_qr(
    vector<TReal, 4> &v0,
    vector<TReal, 4> &v1,
    const TReal v00, const TReal v10, const TReal v01, const TReal v11)
{
    // The factorization was obtained symbolically by WolframAlpha
    // by running the following query:
    //
    // QRDecomposition[{{v_11,v_21},{v_12,v_22},{1,0},{0,1}}]
    v0(0) = v00;
    v0(1) = v01;
    v0(2) = 1;
    v0(3) = 0;
    normalize(v0);

    // OPTME: We could reuse some of the multiplications.
    v1(0) = v10 + v01 * v01 * v10 - v00 * v01 * v11;
    v1(1) = v11 - v00 * v01 * v10 + v00 * v00 * v11;
    v1(2) = -v00 * v10 - v01 * v11;
    v1(3) = v00 * v00 + v01 * v01 + 1;
    normalize(v1);
}

template <typename TReal>
inline void orthonormalize_v_with_qr(
    vector<TReal, 4> &v0,
    vector<TReal, 4> &v1)
{
    // The factorization was obtained symbolically by WolframAlpha
    // by running the following query:
    //
    // QRDecomposition[{{a,e},{b,f},{c,g},{d,h}}]

    normalize(v0);

    // To avoid numerical stability issues when multiplying too big values in the solution,
    // we scale down the second vector.
    TReal factor = v1(0);
    if (factor < math_utils<TReal>::fabs(v1(1)))
        factor = math_utils<TReal>::fabs(v1(1));
    if (factor < math_utils<TReal>::fabs(v1(2)))
        factor = math_utils<TReal>::fabs(v1(2));
    if (factor < math_utils<TReal>::fabs(v1(3)))
        factor = math_utils<TReal>::fabs(v1(3));
    factor = 1 / (factor + (std::numeric_limits<TReal>::min)());

    const TReal a = v0(0);
    const TReal b = v0(1);
    const TReal c = v0(2);
    const TReal d = v0(3);
    const TReal e = v1(0) * factor;
    const TReal f = v1(1) * factor;
    const TReal g = v1(2) * factor;
    const TReal h = v1(3) * factor;

    // OPTME: We could reuse some of the multiplications.
    // ANSME: Is there a more numerically stable way to do this?
    v1(0) = b * b * e + c * c * e + d * d * e - a * b * f - a * c * g - a * d * h;
    v1(1) = -a * b * e + a * a * f + c * c * f + d * d * f - b * c * g - b * d * h;
    v1(2) = -a * c * e - b * c * f + a * a * g + b * b * g + d * d * g - c * d * h;
    v1(3) = -a * d * e - b * d * f - c * d * g + a * a * h + b * b * h + c * c * h;
    normalize(v1);
}

}; // End of namespace detail.

namespace detail
{

// "Hard-coded" constants used by this algorithm.
//
// For the time being, they are same whether single or double precision
// is used.
template <typename TReal>
struct constants
{
    // Tolerance for determinant of matrix B.
    static inline TReal get_tau2();
    // Tolerance for third of determinant of matrix B.
    static inline TReal get_tau1();
    // Tolerance for Newton iterations.
    static inline TReal get_newton_tolerance();
    // Threshold for sub-space iterations.
    static inline TReal get_subspace_threshold();
};

template <typename TReal>
inline TReal constants<TReal>::get_tau2()
{
    return static_cast<TReal>(1.0e-4);
}

template <typename TReal>
inline TReal constants<TReal>::get_tau1()
{
    return static_cast<TReal>(1.0e-4);
}

template <typename TReal>
inline TReal constants<TReal>::get_newton_tolerance()
{
    // The paper mentions 10^-15, but the Matlab implementation uses 10^-12.
    return static_cast<TReal>(1.0e-12);
}

template <typename TReal>
inline TReal constants<TReal>::get_subspace_threshold()
{
    // This constant is used in the test:
    // if log10 |u22| > -7.18
    // This implies that u22 > 10^-7.18 ~= 6.607e-8
    // However, this constant was determined using the machine error on
    // floating-point representation, so it should probably be different
    // between the single and double precision versions.
    return static_cast<TReal>(6.607e-8);
}

}; // End of namespace detail.

// Implementation of 3x3 polar decomposition based on:
//
// "An algorithm to compute the polar decomposition of a 3x3 matrix"
// Nicholas J Higham and Vanni Noferini
// July 2015
//
// The paper is available at:
// http://eprints.ma.man.ac.uk/2352/01/covered/MIMS_ep2015_66.pdf
//
// It is now available at:
// http://link.springer.com/article/10.1007%2Fs11075-016-0098-7
//
// This C++ implementation is also extensively based on the Matlab
// implementation of this algorithm available at:
// https://github.com/higham/polar-decomp-3by3
//
// It is worth noting that Matlab uses a (row,column) matrix indexing,
// but this implementation uses (column,row) indexing.  Also, Matlab is
// 1-index based, but this implementation is 0-index based.
//
// This implementation has very specific goals that drive implementation
// choices:
// - It tries to avoid dependencies to third-party libraries.  Therefore,
//   it reimplements basic operations (matrix operations such as multiplication,
//   transposition, etc.).  These are straight-forward to implement and
//   are kept separate to the actual algorithm so that using an actual
//   linear algebra library would make the implementation more straight-forward.
// - It tries to be as efficient as possible.  Therefore it potentially combines
//   multiple operations into one (for instance, combine matrix transposition
//   with multiplication) to minimize runtime cost.  While this might reduce
//   code simplicity, we favor runtime efficiency while trying to make those
//   optimizations as easy to read as possible.
//
// It is also worth noting that this implementation relies on two major sources:
// - The algorithm as described in the paper
// - The algorithm as implemented in the Matlab implementation.
//
// This implementation tries to highlight its references to both the paper and
// the source code.
//
// The algorithms described in the paper are referred to using ### comments.
// For instance, specific lines such as the beginning of algorithm 3.5 will
// be highlighted:
//
// ### Algorithm 3.5
//
// So that references to the paper are obvious.  References to the Matlab
// implementation of the algorithm will be more textual.
namespace detail
{

template <typename TReal>
inline TReal run_algorithm_3_3(const TReal absDetA, const TReal detB);

template <typename TReal>
inline TReal run_algorithm_3_4(const TReal absDetA, const TReal detB);

// Implementation of algorithm 3.2.
template <typename TReal>
inline void run_algorithm_3_2(
    vector<TReal, 4> &v,
    vector<int, 4> &p,
    const matrix<TReal, 3, 3> &A,
    const matrix<TReal, 4, 4> &B,
    const TReal detB)
{
    // ### 1. Form B in R4? from A via (2.5).

    // ### 2. Compute b = det B from an LU factorization with partial pivoting.
    // Already done.
    const TReal b = detB;

    // ### 3. Compute d = det A from an LU factorization with partial pivoting.
    // Sign of the determinant.
    TReal d;
    // Determinant.
    TReal dd = compute_determinant_lu_partial(A, d);
    assert((d == 1) || (d == -1));

    // ### 4. if d < 0, B = -B, d = -d, end
    // We use the Bs matrix since we will need it anyways.  Bs matrix is formed from
    // minus B matrix.
    matrix<TReal, 4, 4> Bs = B;
    multiply(Bs, -d);
    dd *= d;

    // ### 5. Estimate lambda1, a dominant eigenvalue of B, via Algorithm 3.3.
    const TReal lambda1 = run_algorithm_3_3(dd, b);

    // ### 6. Bs = lambda1 I - B
    // Bs already holds -B.
    Bs(0, 0) += lambda1;
    Bs(1, 1) += lambda1;
    Bs(2, 2) += lambda1;
    Bs(3, 3) += lambda1;

    // ### 7. Compute an LDLT factorization with diagonal pivoting, P^T Bs P = L D L^T.
    matrix<TReal, 4, 4> L;
    vector<TReal, 4> D;
    compute_ldlt_factorization_diagonal(L, D, p, Bs);

    // ### 8. v = PL^-T e4 / ||L^-T e4||2
    // Normalization will be done in the common part.
    v(0) = L(0, 1) * L(1, 3) + L(0, 2) * L(2, 3) - L(0, 1) * L(2, 3) * L(1, 2) - L(0, 3);
    v(1) = L(2, 3) * L(1, 2) - L(1, 3);
    v(2) = -L(2, 3);
    v(3) = 1;

    // ### 9. Form the matrix Q using (2.7).
    // ### 10. Compute the upper triangle of H = Q^T A and set the lower triangle equal to
    //         the upper triangle.
    // Both are done at the end of algorithm 3.5, along with v normalization.
}

// Implementation of algorithm 3.3.
template <typename TReal>
inline TReal run_algorithm_3_3(
    const TReal absDetA,
    const TReal detB)
{
    TReal lambda1;

    const TReal &dd = absDetA;
    const TReal &b = detB;

    // ### 1. tau1 = 10^4 % Tolerance.
    static const TReal kTau1 = constants<TReal>::get_tau1();

    // ### 2. if b + 1/3 > 1
    if (b > kTau1 - 1 / static_cast<TReal>(3))
    {
        // ### 3. c = 8d
        const TReal c = 8 * dd;
        // ### 4. delta0 = 1 + 3b
        const TReal delta0 = 1 + 3 * b;
        // ### 5. delta1 = -1 + (27/16)c^2 + 9b
        const TReal delta1 = -1 + (27 / static_cast<TReal>(16)) * c * c + 9 * b;
        // ### 6. phi = delta1/delta0^(3/2)
        TReal phi = delta1 / (delta0 * math_utils<TReal>::sqrt(delta0));
        // This was not in the original algorithm, but clamp to [-1,1] in case of rounding errors.
        phi = math_utils<TReal>::clamp(phi, -1, 1);
        // ### 7. z = (4/3)(1 + delta0^(1/2)cos(arccos(alpha)/3))
        const TReal z = (4 / static_cast<TReal>(3)) * (1 + math_utils<TReal>::sqrt(delta0) * math_utils<TReal>::cos(math_utils<TReal>::acos(phi) / 3));
        // ### 8. s = z^0.5/2
        const TReal s = math_utils<TReal>::sqrt(z) / 2;
        // ### 9. lambda1 = s + (max(0, 4 - z + c/s))^(1/2)/2.
        lambda1 = s;
        const TReal temp = 4 - z + c / s;
        if (temp > 0)
            lambda1 += math_utils<TReal>::sqrt(temp) / 2;
    }
    // ### 10. else
    else
    {
        // ### 11. Use Newton's method (Algorithm 3.4) to approximate
        lambda1 = run_algorithm_3_4(dd, b);
    }
    // ### 12. end

    return lambda1;
}

// Implementation of algorithm 3.4.
template <typename TReal>
inline TReal run_algorithm_3_4(
    const TReal absDetA,
    const TReal detB)
{
    const TReal &dd = absDetA;
    const TReal &b = detB;

    // ### 1. x = sqrt(3)
    TReal x = math_utils<TReal>::sqrt(3);
    // ### 2. xold = 3
    TReal xold = 3;
    // ### 3. while xold - x > 10^-15
    static const TReal kNewtonTolerance = constants<TReal>::get_newton_tolerance();
    while (xold - x > kNewtonTolerance)
    {
        // ### 4. xold = x
        xold = x;
        // ### 5. Evaluate p = p(x) = det(xI - B) by Horner's method.
        const TReal c = 8 * dd;
        const TReal px = x * (x * (x * x - 2) - c) + b;
        // ### 6. Evaluate pd = p0(x) by Horner's method.
        const TReal dpx = x * (4 * x * x - 4) - c;
        // ### 7. x = x - p/pd
        x = x - px / dpx;
    }
    // ### 8. end
    const TReal lambda1 = x;

    return lambda1;
}

// Implementation of algorithm 3.5.
template <typename TReal>
inline void run_algorithm_3_5(
    matrix<TReal, 3, 3> &paramQ,
    matrix<TReal, 3, 3> &paramH,
    const matrix<TReal, 3, 3> &paramA)
{
    // First make sure the input matrix is normalized.
    matrix<TReal, 3, 3> A = paramA;
    normalize(A);

    // ### Algorithm 3.5

    // ### 1. tau2 = 10^-4 % Tolerance.
    static const TReal kTau2 = constants<TReal>::get_tau2();

    // ### 2. Form B in R4? from A via (2.5).
    matrix<TReal, 4, 4> B;
    // Computation of the matrix as described in the paper's formula.
    // Note: In the author's Matlab implementation, the computation
    // is slightly different.  It first computes the trace as the
    // sum of the diagonal and then computes the diagonal elements
    // of B derived from this.  This can create numerical differences
    // for which the algorithm should be tolerant, but depending on
    // whether single or double precision is used, different paths
    // in the algorithm might be used.  However, they should all yield
    // satisfactory solutions as the algorithm should be numerically stable
    // in both cases.
    //
    // OPTME: B matrix is symmetric, we could store it in a different way.
    B(0, 0) = A(0, 0) + A(1, 1) + A(2, 2);
    B(1, 1) = A(0, 0) - A(1, 1) - A(2, 2);
    B(2, 2) = A(1, 1) - A(0, 0) - A(2, 2);
    B(3, 3) = A(2, 2) - A(0, 0) - A(1, 1);
    B(0, 1) = B(1, 0) = A(2, 1) - A(1, 2);
    B(1, 2) = B(2, 1) = A(1, 0) + A(0, 1);
    B(2, 3) = B(3, 2) = A(2, 1) + A(1, 2);
    B(0, 2) = B(2, 0) = A(0, 2) - A(2, 0);
    B(1, 3) = B(3, 1) = A(2, 0) + A(0, 2);
    B(0, 3) = B(3, 0) = A(1, 0) - A(0, 1);

    // ### 3. Compute b = det B from an LU factorization with partial pivoting.
    // It actually computes B determinant from A matrix is B is a function of A.
    const TReal b = compute_b_determinant_from_a_matrix_lu_partial(A);

    vector<TReal, 4> v;
    vector<int, 4> p;

    // ### 4. if b < 1 - tau2
    if (b < 1 - kTau2)
    {
        // ### % Dominant eigenvalue of B is well separated.
        // ### 5. Call Algorithm 3.2.
        run_algorithm_3_2(v, p, A, B, b);
    }
    // ### 6. else
    else
    {
        // ### 7. Compute d = detA using an LU factorization with complete pivoting.
        // Also keep the second U value of the factorization.
        TReal u22;
        // Sign of the determinant.
        TReal d;
        // Determinant.
        TReal dd = compute_determinant_lu_complete(A, d, u22);
        assert((d == 1) || (d == -1));

        // ### 8. If d < 0, B = -B, end
        // We use the Bs matrix since we will need it anyways.  Bs matrix is formed from
        // minus B matrix.
        matrix<TReal, 4, 4> Bs = B;
        multiply(Bs, -d);
        dd *= d;

        // ### 9. Estimate lambda1 using Algorithm 3.3.
        const TReal lambda1 = run_algorithm_3_3(dd, b);

        // ### 10. Bs = lambda1 I - B
        // Bs already holds -B.
        Bs(0, 0) += lambda1;
        Bs(1, 1) += lambda1;
        Bs(2, 2) += lambda1;
        Bs(3, 3) += lambda1;

        // Compute Bs = LDL^T by block LDL^T factorization with Bunch-Parlett pivoting
        // This will be used in both following cases.
        matrix<TReal, 4, 4> L;
        matrix<TReal, 4, 4> D;
        compute_ldlt_factorization_bunch_parlett(L, D, p, Bs);
        assert(D(2, 3) == D(3, 2));

        const TReal DD = D(2, 2) * D(3, 3) - D(3, 2) * D(3, 2);
        if (DD == 0)
        {
            // Treat this case specially.  It is not really mentioned in the paper's algorithm,
            // but it is part of the Matlab implementation.
            const bool allZero = (D(2, 2) == 0) && (D(3, 3) == 0) && (D(3, 2) == 0);
            if (allZero)
            {
                // This is the equivalent of choosing a null space of (0,1) and do
                // the same calculation as the other case.
                v(0) = L(0, 1) * L(1, 3) - L(0, 3);
                v(1) = -L(1, 3);
                v(2) = 0;
                v(3) = 1;
            }
            else
            {
                // ANSME: A more robust way might be to get into this case for determinant close to 0,
                // instead of exactly equal to zero, and then use something more robust such as taking
                // the vector associated with the smallest singular value (0 if there is actually a
                // null space), which could probably be computed efficiently for 2x2 matrices.
                const vector<TReal, 2> nullSpace = compute_null_space(D(2, 2), D(2, 3), D(3, 3));

                // Since L is diagonal, we can solve v = L^-T * [0 0 a b]^T.
                // We also know, from compute_ldlt_factorization_bunch_parlett, that L(2,3) is 0.
                // So v can be computed symbolically by WolframAlpha by running the following query:
                // Inverse[Transpose[{{1,0,0,0},{l_12,1,0,0},{l_13,l_23,1,0},{l_14,l_24,0,1}}]] * {{0},{0},{a},{b}}
                //
                // We have L^-T = | 1    -L(0,1)    L(0,1)*L(1,2) - L(0,2)    -L(0,3) + L(0,2)*L(2,3) + L(0,1)*(L(1,3) - L(1,2)*L(2,3)) |
                //                | 0    1          -L(1,2)                   L(1,2)*L(2,3) - L(1,3)                                    |
                //                | 0    0          1                         -L(2,3)                                                   |
                //                | 0    0          0                         1                                                         |
                //
                // Using L(2,3) == 0 we get:
                //
                // We have L^-T = | 1    -L(0,1)    L(0,1)*L(1,2) - L(0,2)    L(0,1)*L(1,3) -L(0,3) |
                //                | 0    1          -L(1,2)                   -L(1,3)               |
                //                | 0    0          1                         0                     |
                //                | 0    0          0                         1                     |
                //
                // Multiplied by [0 0 a b]^T:
                //
                // | a * (L(0,1)*L(1,2) - L(0,2)) + b * (L(0,1)*L(1,3) - L(0,3)) |
                // | -a * L(1,2) - b * L(1,3)                                    |
                // | a                                                           |
                // | b                                                           |

                v(0) = nullSpace(0) * (L(0, 1) * L(1, 2) - L(0, 2)) + nullSpace(1) * (L(0, 1) * L(1, 3) - L(0, 3));
                v(1) = -nullSpace(0) * L(1, 2) - nullSpace(1) * L(1, 3);
                v(2) = nullSpace(0);
                v(3) = nullSpace(1);
            }
        }
        else
        {
            // Compute inverse of L.
            // See above for explanation of L^-1 computation (and assumptions about which values are 0 and 1).
            const TReal IL01 = -L(0, 1);
            const TReal IL02 = L(0, 1) * L(1, 2) - L(0, 2);
            const TReal IL12 = -L(1, 2);
            const TReal IL03 = L(0, 1) * L(1, 3) - L(0, 3);
            const TReal IL13 = -L(1, 3);

            // Compute inverse of D.
            // Inverse[{{d_11,0,0,0},{0,d_22,0,0},{0,0,d_33,d_43},{0,0,d_43,d_44}}]
            const TReal ID00 = 1 / D(0, 0);
            const TReal ID11 = 1 / D(1, 1);
            matrix<TReal, 2, 2> ID;
            ID(0, 0) = D(3, 3);
            ID(0, 1) = -D(3, 2);
            ID(1, 0) = -D(3, 2);
            ID(1, 1) = D(2, 2);
            multiply(ID, 1 / DD);

            // ### 11. if log10 |u22| > -7.18
            // This implies that u22 > 10^-7.18 ~= 6.607e-8
            static const TReal kSubspaceThreshold = constants<TReal>::get_subspace_threshold();
            const TReal AU = math_utils<TReal>::fabs(u22);
            if (AU > kSubspaceThreshold)
            {
                // ### 12. nit = ceil(15/(16.86 + 2 log10 |u22|))
                const int nit = static_cast<int>(
                    math_utils<TReal>::ceil(15 / (static_cast<TReal>(16.8) + 2 * math_utils<TReal>::log10(AU))));

                // ### 13. Compute Bs = LDL^T by block LDL^T factorization with Bunch-Parlett pivoting
                // Already done.

                // ### 14. v = L^-T e4 / ||L^-T e4|| % Initial guess.
                //
                // | 1       IL01    IL02    IL03 |       | 0 |       | IL03 |
                // | 0       1       IL12    IL13 |   *   | 0 |   =   | IL13 |
                // | 0       0       1       0    |       | 0 |       | 0    |
                // | 0       0       0       1    |       | 1 |       | 1    |
                v(0) = IL03;
                v(1) = IL13;
                v(2) = 0;
                v(3) = 1;
                // Normalization will happen in the loop.

                // ### 15. for i = 1 : nit
                for (int i = 0; i < nit; ++i)
                {
                    normalize(v);

                    // ### 16. Update v using one step of inverse iteration with LDL^T
                    // OPTME: Maybe some of these operations could be combined to be optimized...?

                    // v = L^-1 * v = IL * v;
                    v = multiply_il_v(IL01, IL02, IL03, IL12, IL13, v);

                    // v = D^-1 = ID * v;
                    v = multiply_id_v(ID00, ID11, ID, v);

                    // v = L^-T * v = IL^T * v = v * IL;
                    v = multiply_v_il(v, IL01, IL02, IL03, IL12, IL13);
                }
                // ### 17. end
                // The last normalization of v will be done at the end.
            }
            // ### 18. else
            else
            {
                // ### 19. Compute Bs = LDL^T by block LDL^T factorization with Bunch-Parlett pivoting
                // Already done.

                // ### 20. V = L^-T [e3 e4] % Initial guess.
                //
                // | 1       IL01    IL02    IL03 |       | 0    0 |       | IL02    IL03 |
                // | 0       1       IL12    IL13 |   *   | 0    0 |   =   | IL12    IL13 |
                // | 0       0       1       0    |       | 1    0 |       | 1       0    |
                // | 0       0       0       1    |       | 0    1 |       | 0       1    |
                const TReal v00 = IL02;
                const TReal v10 = IL03;
                const TReal v01 = IL12;
                const TReal v11 = IL13;

                // ### 21. for i = 1:2
                // ### 22. Orthonormalize V via QR factorization.
                // ### 23. Update V using one step of inverse subspace iteration with LDL^T.
                vector<TReal, 4> v0, v1;
                orthonormalize_v_with_qr(v0, v1, v00, v10, v01, v11);

                for (int i = 0; i < 2; ++i)
                {
                    v0 = multiply_il_v(IL01, IL02, IL03, IL12, IL13, v0);
                    v1 = multiply_il_v(IL01, IL02, IL03, IL12, IL13, v1);

                    v0 = multiply_id_v(ID00, ID11, ID, v0);
                    v1 = multiply_id_v(ID00, ID11, ID, v1);

                    v0 = multiply_v_il(v0, IL01, IL02, IL03, IL12, IL13);
                    v1 = multiply_v_il(v1, IL01, IL02, IL03, IL12, IL13);
                }
                // ### 24. end

                // ### 25. Orthonormalize V via QR factorization.
                orthonormalize_v_with_qr(v0, v1);

                // ### 26. Bp = V^T Bs V in R2x2
                // ### 27. Find w, eigenvector of smallest eigenvalue of Bp, by analytic formula.
                // ### 28. v = V w
                // L has the same form as IL, so we can use the same function to multiply them.
                const vector<TReal, 4> v0_temp = multiply_v_il(v0, L(0, 1), L(0, 2), L(0, 3), L(1, 2), L(1, 3));
                const vector<TReal, 4> v1_temp = multiply_v_il(v1, L(0, 1), L(0, 2), L(0, 3), L(1, 2), L(1, 3));
                const vector<TReal, 4> H0 = multiply_minus_v_d(v0_temp, D);
                const vector<TReal, 4> H1 = multiply_minus_v_d(v1_temp, D);
                const TReal H00 = dot(H0, v0_temp);
                const TReal H10 = dot(H0, v1_temp);
                const TReal H11 = dot(H1, v1_temp);
                if (math_utils<TReal>::fabs(H10) < static_cast<TReal>(1.0e-15))
                {
                    if (H00 > H10)
                        v = v0;
                    else
                        v = v1;
                }
                else
                {
                    const TReal r = (H00 - H11) / (2 * H10);
                    const int s = (H10 < 0 ? -1 : 1);
                    const TReal f = r + s * math_utils<TReal>::sqrt(1 + r * r);
                    v(0) = v0(0) * f + v1(0);
                    v(1) = v0(1) * f + v1(1);
                    v(2) = v0(2) * f + v1(2);
                    v(3) = v0(3) * f + v1(3);
                }
            }
            // ### 29. end
        }

        // ### 30. Form the matrix Q from v as in Theorem 2.5.
        // ### 31. Compute H = Q^T A.
        // Both are done at the end of the function.
    }

    // Compute rotation from dominant eigen vector v.
    normalize(v);
    vector<TReal, 4> vtemp = v;
    v(p(0)) = vtemp(0);
    v(p(1)) = vtemp(1);
    v(p(2)) = vtemp(2);
    v(p(3)) = vtemp(3);

    const TReal v12 = 2 * v(0) * v(1);
    const TReal v13 = 2 * v(0) * v(2);
    const TReal v14 = 2 * v(0) * v(3);
    const TReal v22 = 2 * v(1) * v(1);
    const TReal v23 = 2 * v(1) * v(2);
    const TReal v24 = 2 * v(1) * v(3);
    const TReal v33 = 2 * v(2) * v(2);
    const TReal v34 = 2 * v(2) * v(3);
    const TReal v44 = 2 * v(3) * v(3);

    paramQ(0, 0) = 1 - (v33 + v44);
    paramQ(0, 1) = v23 - v14;
    paramQ(0, 2) = v24 + v13;
    paramQ(1, 0) = v23 + v14;
    paramQ(1, 1) = 1 - (v22 + v44);
    paramQ(1, 2) = v34 - v12;
    paramQ(2, 0) = v24 - v13;
    paramQ(2, 1) = v34 + v12;
    paramQ(2, 2) = 1 - (v22 + v33);

    // The Matlab implementation returns the opposite of the matrix if det A < 0.
    // We don't do that because we want a right-handed rotation.

    // Compute scale.
    transpose_multiply(paramH, paramQ, paramA);

    // The Matlab implementation suggests averaging the top and lower part of the
    // matrix to ensure symmetry, but we don't do it.
}

}; // End of namespace detail.

}; // End of namespace polar.

#endif // __POLAR_DECOMPOSITION_3X3_IMPL_H__
