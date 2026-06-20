#ifndef ELECTROMAGNETIC_OPHELIE_PHI_RHS_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_PHI_RHS_DIAGNOSTICS_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_parameters.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OpheliePhiRhsFingerprint
{
    Real volume_sum = 0.0;
    Real vol_weighted_l2 = 0.0;
    Real min_value = 0.0;
    Real max_value = 0.0;
    std::uint64_t xor_checksum = 0;
};

inline std::uint64_t ophelieScalarXorChecksum(Real value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(Real) == sizeof(std::uint64_t) || sizeof(Real) == sizeof(std::uint32_t),
                  "Unexpected Real size for RHS checksum");
    if constexpr (sizeof(Real) == sizeof(std::uint64_t))
    {
        std::memcpy(&bits, &value, sizeof(Real));
    }
    else
    {
        std::uint32_t bits32 = 0;
        std::memcpy(&bits32, &value, sizeof(Real));
        bits = static_cast<std::uint64_t>(bits32);
    }
    return bits;
}

inline OpheliePhiRhsFingerprint computeOpheliePhiRhsFingerprint(BaseParticles &particles,
                                                                const std::string &rhs_field_name,
                                                                size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, rhs_field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(rhs_field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OpheliePhiRhsFingerprint fingerprint;
    if (total_real_particles == 0)
    {
        return fingerprint;
    }

    fingerprint.min_value = rhs[0];
    fingerprint.max_value = rhs[0];
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        fingerprint.volume_sum += vol[i] * rhs[i];
        fingerprint.vol_weighted_l2 += vol[i] * rhs[i] * rhs[i];
        fingerprint.min_value = std::min(fingerprint.min_value, rhs[i]);
        fingerprint.max_value = std::max(fingerprint.max_value, rhs[i]);
        fingerprint.xor_checksum ^= ophelieScalarXorChecksum(rhs[i]);
    }
    fingerprint.vol_weighted_l2 = std::sqrt(std::max(fingerprint.vol_weighted_l2, Real(0)));
    return fingerprint;
}

inline bool opheliePhiRhsFingerprintsMatch(const OpheliePhiRhsFingerprint &lhs,
                                           const OpheliePhiRhsFingerprint &rhs, Real l2_tol = Real(1.0e-12),
                                           Real sum_tol = Real(1.0e-12))
{
    return lhs.xor_checksum == rhs.xor_checksum &&
           std::abs(lhs.volume_sum - rhs.volume_sum) <= sum_tol &&
           std::abs(lhs.vol_weighted_l2 - rhs.vol_weighted_l2) <= l2_tol;
}

inline void logOpheliePhiRhsFingerprint(const char *stage, const OpheliePhiRhsFingerprint &fingerprint)
{
    std::cout << "[ophelie] phi_rhs_fingerprint stage=" << stage << " sum=" << fingerprint.volume_sum
              << " l2=" << fingerprint.vol_weighted_l2 << " min=" << fingerprint.min_value
              << " max=" << fingerprint.max_value << " xor_checksum=" << fingerprint.xor_checksum << std::endl;
}

inline Real computeOphelieEdgeFluxInputScaleFromRhsL2(Real rhs_l2, Real safe_rhs_l2 = Real(1.0e4))
{
    if (!std::isfinite(rhs_l2) || rhs_l2 <= safe_rhs_l2)
    {
        return 1.0;
    }
    return safe_rhs_l2 / rhs_l2;
}

inline Real computeOphelieEdgeFluxInputScaleFromRhsMax(Real rhs_max_abs, Real safe_rhs_max_abs = Real(1.2e3))
{
    if (!std::isfinite(rhs_max_abs) || rhs_max_abs <= safe_rhs_max_abs)
    {
        return 1.0;
    }
    return safe_rhs_max_abs / rhs_max_abs;
}

/** Ignore float-overflow spikes when choosing phi solver-local RHS scale (p99.5 of |rhs|). */
inline Real hostRobustPhiRhsMaxAbs(BaseParticles &particles, const std::string &rhs_field_name, size_t n,
                                  Real percentile = Real(0.995))
{
    syncVariableToHost<Real>(particles, rhs_field_name);
    const Real *rhs = particles.getVariableDataByName<Real>(rhs_field_name);
    if (n == 0)
    {
        return 0.0;
    }
    StdVec<Real> abs_values;
    abs_values.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const Real value = std::abs(rhs[i]);
        if (std::isfinite(value))
        {
            abs_values.push_back(value);
        }
    }
    if (abs_values.empty())
    {
        return 0.0;
    }
    std::sort(abs_values.begin(), abs_values.end());
    const Real rank = percentile * static_cast<Real>(abs_values.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(rank));
    const size_t upper = static_cast<size_t>(std::ceil(rank));
    const Real weight = rank - static_cast<Real>(lower);
    return abs_values[lower] * (Real(1) - weight) + abs_values[upper] * weight;
}

inline void hostScaleScalarFieldInPlace(BaseParticles &particles, const std::string &field_name, Real scale, size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    Real *values = particles.getVariableDataByName<Real>(field_name);
    for (size_t i = 0; i < n; ++i)
    {
        values[i] *= scale;
    }
    syncVariableToDevice<Real>(particles, field_name);
}

/** Solver-local scale: multiply PCG workspace RHS by s, then unscale phi solution by 1/s. */
inline Real computeOpheliePhiSolverLocalRhsScale(BaseParticles &particles, const std::string &rhs_field_name, size_t n,
                                                 Real safe_rhs_l2, Real safe_rhs_max_abs)
{
    const OpheliePhiRhsFingerprint rhs_fp = computeOpheliePhiRhsFingerprint(particles, rhs_field_name, n);
    const Real rhs_max_abs = hostRobustPhiRhsMaxAbs(particles, rhs_field_name, n);
    const Real scale_l2 = computeOphelieEdgeFluxInputScaleFromRhsL2(rhs_fp.vol_weighted_l2, safe_rhs_l2);
    const Real scale_max = computeOphelieEdgeFluxInputScaleFromRhsMax(rhs_max_abs, safe_rhs_max_abs);
    return std::max(scale_l2, scale_max);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_RHS_DIAGNOSTICS_H
