/**
 * @file electromagnetic_ophelie_fixed_space_probe.h
 * @brief Doc232 — shared fixed-space temperature sampling (production + gates).
 *
 * Wendland-C2-like compact weight; samples Temperature at a fixed Eulerian point.
 * Must not track particle IDs. Used by Doc231 B1 gates and Doc232 production timeseries.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_FIXED_SPACE_PROBE_H
#define ELECTROMAGNETIC_OPHELIE_FIXED_SPACE_PROBE_H

#include "base_particles.h"
#include "electromagnetic_ophelie_device_sync.h"

#include <cmath>
#include <limits>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real ophelieWendlandProbeWeight(Real r, Real h)
{
    if (!(h > TinyReal))
        return Real(0);
    const Real q = r / h;
    if (q >= Real(2))
        return Real(0);
    const Real t = Real(1) - Real(0.5) * q;
    return t * t * t * t * (Real(1) + Real(2) * q);
}

struct OphelieFixedSpaceProbeSample
{
    Real T = std::numeric_limits<Real>::quiet_NaN();
    Real weight_sum = Real(0);
    int n_neighbors = 0;
    Real nearest_dist = Real(-1);
    bool valid = false;
};

/**
 * Host kernel sample of Temperature at probe (Eulerian).
 * Caller must sync Position/Temperature to host first if device-dirty.
 */
inline OphelieFixedSpaceProbeSample hostSampleFixedSpaceTemperature(BaseParticles &particles, const Vecd &probe,
                                                                    Real kernel_h,
                                                                    const std::string &temperature_field = "Temperature")
{
    OphelieFixedSpaceProbeSample s;
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    const size_t n = particles.TotalRealParticles();
    Real t_sum = Real(0);
    Real best = MaxReal;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(static_cast<double>(pos[i][0])))
            continue;
        const Real r = (pos[i] - probe).norm();
        best = std::min(best, r);
        const Real w = ophelieWendlandProbeWeight(r, kernel_h);
        if (w <= TinyReal)
            continue;
        s.weight_sum += w;
        t_sum += w * temp[i];
        ++s.n_neighbors;
    }
    s.nearest_dist = (best < MaxReal * Real(0.5)) ? best : Real(-1);
    if (s.weight_sum > TinyReal && s.n_neighbors > 0 && std::isfinite(static_cast<double>(t_sum)))
    {
        s.T = t_sum / s.weight_sum;
        s.valid = std::isfinite(static_cast<double>(s.T));
    }
    return s;
}

inline void syncFixedSpaceProbeFieldsToHost(BaseParticles &particles,
                                            const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, temperature_field);
}

struct OphelieFixedSpaceProbeDef
{
    std::string id;
    Vecd pos;
    std::string location_status = "development"; // or unconfirmed_sensor_location
    bool expect_interior = true;
    std::string note;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FIXED_SPACE_PROBE_H
