/**
 * @file electromagnetic_ophelie_french_skin_depth_diagnostic.h
 * @brief Round A — French skin-depth estimate + radial Joule power shells (host).
 *
 * Does not claim an analytic cylinder match. Reports:
 *   δ = sqrt(2/(ω μ0 σ)), δ/R, δ/dp,
 *   simulated power fraction in the outer shell of thickness min(δ,R),
 *   and a plane-skin exponential proxy for the same shell.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_SKIN_DEPTH_DIAGNOSTIC_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_SKIN_DEPTH_DIAGNOSTIC_H

#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real frenchAnalyticSkinDepthM(Real frequency_hz, Real sigma, Real mu0 = Real(4.0e-7) * Pi)
{
    if (frequency_hz <= TinyReal || sigma <= TinyReal || mu0 <= TinyReal)
    {
        return Real(0);
    }
    const Real omega = Real(2) * Pi * frequency_hz;
    return std::sqrt(Real(2) / (omega * mu0 * sigma));
}

/** Plane-skin proxy: q ∝ exp(-2 x / δ), x = depth from surface. Fraction in [0, s]. */
inline Real frenchPlaneSkinPowerFractionOuter(Real shell_depth_m, Real delta_m, Real radius_m)
{
    if (delta_m <= TinyReal || radius_m <= TinyReal)
    {
        return Real(0);
    }
    const Real s = std::min(std::max(shell_depth_m, Real(0)), radius_m);
    const Real denom = Real(1) - std::exp(Real(-2) * radius_m / delta_m);
    if (denom <= TinyReal)
    {
        return s / radius_m;
    }
    return (Real(1) - std::exp(Real(-2) * s / delta_m)) / denom;
}

struct FrenchSkinDepthRadialReport
{
    Real frequency_hz = 0.0;
    Real sigma = 0.0;
    Real mu0 = 0.0;
    Real delta_m = 0.0;
    Real radius_m = 0.0;
    Real dp_m = 0.0;
    Real delta_over_R = 0.0;
    Real delta_over_dp = 0.0;
    Real p_total_w = 0.0;
    Real p_outer_1delta_w = 0.0;
    Real p_outer_1delta_frac = 0.0;
    Real p_center_quarter_w = 0.0;
    Real p_center_quarter_frac = 0.0;
    Real p_outer_075R_w = 0.0;
    Real p_outer_075R_frac = 0.0;
    Real plane_skin_outer_1delta_frac = 0.0;
    Real q_mean_outer_075R = 0.0;
    Real q_mean_center_025R = 0.0;
    Real q_outer_over_center = 0.0;
    std::string regime = "undefined";
};

inline std::string frenchSkinPenetrationRegime(Real delta_over_R)
{
    if (!(delta_over_R > Real(0)))
    {
        return "undefined";
    }
    if (delta_over_R < Real(0.1))
    {
        return "thin_skin";
    }
    if (delta_over_R < Real(0.3))
    {
        return "intermediate";
    }
    return "deep_penetration";
}

inline FrenchSkinDepthRadialReport computeFrenchSkinDepthRadialReport(
    const Vecd *pos, const Real *q, const Real *vol, size_t n, const FrenchCylindricalFrame &frame, Real frequency_hz,
    Real sigma, Real mu0, Real dp_m)
{
    FrenchSkinDepthRadialReport r;
    r.frequency_hz = frequency_hz;
    r.sigma = sigma;
    r.mu0 = mu0;
    r.radius_m = frame.radius;
    r.dp_m = dp_m;
    r.delta_m = frenchAnalyticSkinDepthM(frequency_hz, sigma, mu0);
    r.delta_over_R = r.delta_m / (r.radius_m + TinyReal);
    r.delta_over_dp = r.delta_m / (dp_m + TinyReal);
    r.regime = frenchSkinPenetrationRegime(r.delta_over_R);

    const Real shell = std::min(r.delta_m, r.radius_m);
    const Real r_outer_skin = r.radius_m - shell;
    const Real r_outer_075 = Real(0.75) * r.radius_m;
    const Real r_center_025 = Real(0.25) * r.radius_m;

    double p_tot = 0.0;
    double p_skin = 0.0;
    double p_ctr = 0.0;
    double p_o75 = 0.0;
    double q_o_acc = 0.0;
    double q_o_w = 0.0;
    double q_c_acc = 0.0;
    double q_c_w = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        if (vol[i] <= TinyReal || !std::isfinite(q[i]) || q[i] < Real(0))
        {
            continue;
        }
        Real rr = 0.0;
        Real zz = 0.0;
        frenchCylindricalRZ(pos[i], frame, rr, zz);
        const double pi = static_cast<double>(q[i]) * static_cast<double>(vol[i]);
        const double wi = static_cast<double>(vol[i]);
        p_tot += pi;
        if (rr >= r_outer_skin)
        {
            p_skin += pi;
        }
        if (rr <= r_center_025)
        {
            p_ctr += pi;
            q_c_acc += static_cast<double>(q[i]) * wi;
            q_c_w += wi;
        }
        if (rr >= r_outer_075)
        {
            p_o75 += pi;
            q_o_acc += static_cast<double>(q[i]) * wi;
            q_o_w += wi;
        }
    }

    r.p_total_w = static_cast<Real>(p_tot);
    r.p_outer_1delta_w = static_cast<Real>(p_skin);
    r.p_center_quarter_w = static_cast<Real>(p_ctr);
    r.p_outer_075R_w = static_cast<Real>(p_o75);
    r.p_outer_1delta_frac = static_cast<Real>(p_skin / (p_tot + TinyReal));
    r.p_center_quarter_frac = static_cast<Real>(p_ctr / (p_tot + TinyReal));
    r.p_outer_075R_frac = static_cast<Real>(p_o75 / (p_tot + TinyReal));
    r.q_mean_outer_075R = static_cast<Real>(q_o_acc / (q_o_w + TinyReal));
    r.q_mean_center_025R = static_cast<Real>(q_c_acc / (q_c_w + TinyReal));
    r.q_outer_over_center = r.q_mean_outer_075R / (r.q_mean_center_025R + TinyReal);
    r.plane_skin_outer_1delta_frac = frenchPlaneSkinPowerFractionOuter(shell, r.delta_m, r.radius_m);
    return r;
}

inline void logFrenchSkinDepthRadialReport(const FrenchSkinDepthRadialReport &r, const std::string &tag = "skin-depth")
{
    std::cout << "[ophelie][" << tag << "] f_hz=" << r.frequency_hz << " sigma=" << r.sigma << " mu0=" << r.mu0
              << " delta_m=" << r.delta_m << " R=" << r.radius_m << " dp=" << r.dp_m
              << " delta/R=" << r.delta_over_R << " delta/dp=" << r.delta_over_dp << " regime=" << r.regime
              << std::endl;
    std::cout << "[ophelie][" << tag << "] P_total_W=" << r.p_total_w
              << " P_outer_1delta_frac=" << r.p_outer_1delta_frac
              << " plane_skin_proxy_frac=" << r.plane_skin_outer_1delta_frac
              << " P_center_r<0.25R_frac=" << r.p_center_quarter_frac
              << " P_outer_r>0.75R_frac=" << r.p_outer_075R_frac << " Q_outer/center=" << r.q_outer_over_center
              << std::endl;
    if (r.regime == "deep_penetration")
    {
        std::cout << "[ophelie][" << tag
                  << "] note: delta~R — thin-skin plane proxy is weak; treat shell fractions as SPH labels, "
                     "not Jacoutot match."
                  << std::endl;
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_SKIN_DEPTH_DIAGNOSTIC_H
