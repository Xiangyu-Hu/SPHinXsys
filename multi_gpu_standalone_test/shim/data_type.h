// Minimal stand-in for SPHinXsys data_type.h, sufficient to compile and exercise
// DomainPartition without Eigen/TBB/Simbody. Types mirror the real API surface used.
#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace SPH
{
using Real = double;
using UnsignedInt = size_t;
template <typename T> using StdVec = std::vector<T>;
constexpr int Dimensions = 3;

template <typename T, int N>
struct FixedVec
{
    T v[N]{};
    T &operator[](int i) { return v[i]; }
    const T &operator[](int i) const { return v[i]; }
    static FixedVec Zero() { return FixedVec(); }
    static FixedVec Constant(T c) { FixedVec r; for (int i = 0; i < N; ++i) r.v[i] = c; return r; }
    FixedVec operator+(const FixedVec &o) const { FixedVec r; for (int i = 0; i < N; ++i) r.v[i] = v[i] + o.v[i]; return r; }
    FixedVec operator-(const FixedVec &o) const { FixedVec r; for (int i = 0; i < N; ++i) r.v[i] = v[i] - o.v[i]; return r; }
};

using Vecd = FixedVec<Real, Dimensions>;
using Arrayi = FixedVec<int, Dimensions>;

struct BoundingBoxd
{
    Vecd lower_, upper_;
    BoundingBoxd() = default;
    BoundingBoxd(const Vecd &lower, const Vecd &upper) : lower_(lower), upper_(upper) {}
    bool checkContain(const Vecd &p) const
    {
        for (int i = 0; i < Dimensions; ++i)
            if (p[i] < lower_[i] || p[i] > upper_[i]) return false;
        return true;
    }
};

template <class T> inline T SMAX(T a, T b) { return a > b ? a : b; }
template <class T> inline T SMIN(T a, T b) { return a < b ? a : b; }
} // namespace SPH
