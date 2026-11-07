/**
 * @file electromagnetic_ophelie_analytic_induction_weak_feedback.h
 * @brief Batch E — weak-feedback uniform axial B: E_θ ∝ r, Q ∝ r² (host analytics).
 *
 * Prescribed A = ½ B₀ (−y, x, 0) ⇒ B = B₀ ê_z, div A = 0.
 * With φ ≈ 0: E = −i ω A, |E| = ω B₀ r / 2, Q = ½ σ |E|² ∝ r².
 * Volume-mean Q ratio (r>0.75R)/(r<0.25R) = 25 for a right circular cylinder
 * of equal height coverage (doc 23 §14.1).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_ANALYTIC_INDUCTION_WEAK_FEEDBACK_H
#define ELECTROMAGNETIC_OPHELIE_ANALYTIC_INDUCTION_WEAK_FEEDBACK_H

#include "electromagnetic_ophelie_device_sync.h"
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

/** Rotational gauge for uniform B = B₀ ê_z. */
inline Vecd ophelieRotationalAUniformBz(const Vecd &pos, const Vecd &axis_origin, Real b0)
{
    const Real x = pos[0] - axis_origin[0];
    const Real y = pos[1] - axis_origin[1];
    return Vecd(-Real(0.5) * b0 * y, Real(0.5) * b0 * x, Real(0));
}

inline Real ophelieWeakFeedbackExactQ(Real r_cyl, Real omega, Real b0, Real sigma)
{
    const Real e_mag = Real(0.5) * omega * b0 * r_cyl;
    return Real(0.5) * sigma * e_mag * e_mag;
}

/** Analytic volume-mean ⟨Q⟩ ratio outer (r>0.75R) / center (r<0.25R) = 25. */
inline Real ophelieWeakFeedbackAnalyticQOuterOverCenter()
{
    return Real(25);
}

struct OphelieWeakFeedbackQr2ShellMetrics
{
    Real radius_m = 0.0;
    Real p_total_w = 0.0;
    Real p_analytic_w = 0.0;
    Real q_mean_outer_075R = 0.0;
    Real q_mean_center_025R = 0.0;
    Real q_outer_over_center = 0.0;
    Real analytic_outer_over_center = 25.0;
    Real corr_q_r2 = 0.0; // vol-weighted correlation of Q with r²
    Real rel_power_err = 0.0;
    size_t n_outer = 0;
    size_t n_center = 0;
};

inline OphelieWeakFeedbackQr2ShellMetrics hostOphelieWeakFeedbackQr2ShellMetrics(
    const Vecd *pos, const Real *q, const Real *vol, size_t n, const Vecd &axis_origin, Real radius_m, Real omega,
    Real b0, Real sigma)
{
    OphelieWeakFeedbackQr2ShellMetrics m;
    m.radius_m = radius_m;
    m.analytic_outer_over_center = ophelieWeakFeedbackAnalyticQOuterOverCenter();
    const Real r_o = Real(0.75) * radius_m;
    const Real r_c = Real(0.25) * radius_m;

    double p_tot = 0.0;
    double p_exact = 0.0;
    double q_o = 0.0, v_o = 0.0;
    double q_c = 0.0, v_c = 0.0;
    double sum_q = 0.0, sum_r2 = 0.0, sum_q_r2 = 0.0, sum_q2 = 0.0, sum_r4 = 0.0, sum_v = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real r2 = r * r;
        const Real qi = q[i];
        const Real vi = vol[i];
        p_tot += static_cast<double>(qi) * static_cast<double>(vi);
        p_exact += static_cast<double>(ophelieWeakFeedbackExactQ(r, omega, b0, sigma)) * static_cast<double>(vi);
        sum_v += vi;
        sum_q += qi * vi;
        sum_r2 += r2 * vi;
        sum_q_r2 += qi * r2 * vi;
        sum_q2 += qi * qi * vi;
        sum_r4 += r2 * r2 * vi;
        if (r > r_o)
        {
            q_o += qi * vi;
            v_o += vi;
            ++m.n_outer;
        }
        if (r < r_c)
        {
            q_c += qi * vi;
            v_c += vi;
            ++m.n_center;
        }
    }

    m.p_total_w = static_cast<Real>(p_tot);
    m.p_analytic_w = static_cast<Real>(p_exact);
    m.rel_power_err = std::abs(m.p_total_w - m.p_analytic_w) / (std::abs(m.p_analytic_w) + TinyReal);
    m.q_mean_outer_075R = v_o > TinyReal ? static_cast<Real>(q_o / v_o) : Real(0);
    m.q_mean_center_025R = v_c > TinyReal ? static_cast<Real>(q_c / v_c) : Real(0);
    m.q_outer_over_center = m.q_mean_outer_075R / (m.q_mean_center_025R + TinyReal);

    const double inv_v = sum_v > TinyReal ? 1.0 / sum_v : 0.0;
    const double mean_q = sum_q * inv_v;
    const double mean_r2 = sum_r2 * inv_v;
    const double cov = sum_q_r2 * inv_v - mean_q * mean_r2;
    const double var_q = sum_q2 * inv_v - mean_q * mean_q;
    const double var_r2 = sum_r4 * inv_v - mean_r2 * mean_r2;
    m.corr_q_r2 = (var_q > TinyReal && var_r2 > TinyReal)
                      ? static_cast<Real>(cov / std::sqrt(var_q * var_r2))
                      : Real(0);
    return m;
}

inline bool verifyOphelieWeakFeedbackAnalyticIdentities(std::string &message)
{
    // ⟨r²⟩_{r>0.75R} / ⟨r²⟩_{r<0.25R} with area measure r dr
    const Real R = Real(1);
    const Real num = Real(0.5) * (std::pow(R, 4) - std::pow(Real(0.75) * R, 4)) / (std::pow(R, 2) - std::pow(Real(0.75) * R, 2));
    const Real den = Real(0.5) * std::pow(Real(0.25) * R, 2);
    const Real ratio = num / (den + TinyReal);
    if (!(std::abs(ratio - Real(25)) <= Real(1.0e-9) * Real(25)))
    {
        message = "analytic_q_ratio_not_25";
        return false;
    }
    if (!(std::abs(ophelieWeakFeedbackAnalyticQOuterOverCenter() - Real(25)) <= TinyReal))
    {
        message = "helper_ratio_constant_wrong";
        return false;
    }
    message = "weak_feedback_analytic_ok";
    return true;
}

inline void logOphelieWeakFeedbackQr2ShellMetrics(const OphelieWeakFeedbackQr2ShellMetrics &m, const char *tag)
{
    std::cout << "[ophelie][weak-feedback-Qr2] " << tag << " R=" << m.radius_m
              << " Q_o/c=" << m.q_outer_over_center << " analytic=25"
              << " corr(Q,r2)=" << m.corr_q_r2 << " P_sph=" << m.p_total_w << " P_analytic=" << m.p_analytic_w
              << " rel_P=" << m.rel_power_err << " n_o/c=" << m.n_outer << "/" << m.n_center << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_ANALYTIC_INDUCTION_WEAK_FEEDBACK_H
