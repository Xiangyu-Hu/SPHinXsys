/**
 * @file electromagnetic_ophelie_analytic_cylinder_self_induction.h
 * @brief Doc 23 §14.1 — infinite conducting cylinder in uniform axial AC B₀ (self-induction).
 *
 * Convention: exp(+iωt), E = −iωA (matches SPH / GetDP metadata). Non-magnetic cylinder
 * (μ=μ₀) with impressed far-field / surface H_z(R)=H₀=B₀/μ₀:
 *   α² = −i ω μ σ ,  α = (1−i)/δ ,  δ = √(2/(ω μ σ))
 *   H_z(r) = H₀ J₀(α r) / J₀(α R)     (r ≤ R)
 *   E_θ(r) = (α/σ) H₀ J₁(α r) / J₀(α R)
 *   J_θ = σ E_θ ,  Q = ½ σ |E|²
 *
 * Weak-feedback limit |α R|→0 recovers E_θ = −i ω B₀ r / 2 and Q ∝ r².
 * (Doc 40 R09: previous α=(1+i)/δ was the complex conjugate under this convention.)
 * Finite-height SPH compares mid-plane radial profiles + P/H ≈ P' (W/m).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_SELF_INDUCTION_H
#define ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_SELF_INDUCTION_H

#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

using OphelieComplex = std::complex<Real>;

inline Real ophelieSkinDepthM(Real omega, Real mu0, Real sigma)
{
    return std::sqrt(Real(2) / (omega * mu0 * sigma + TinyReal));
}

/** Propagation constant α = √(−i ω μ σ) = (1−i)/δ under exp(+iωt). */
inline OphelieComplex ophelieCylinderAlpha(Real omega, Real mu0, Real sigma)
{
    const Real delta = ophelieSkinDepthM(omega, mu0, sigma);
    return OphelieComplex(Real(1) / delta, Real(-1) / delta);
}

/** Power-series J₀(z) for complex z (enough for |z| ≲ 20). */
inline OphelieComplex ophelieBesselJ0(OphelieComplex z)
{
    const OphelieComplex z2 = z * z;
    OphelieComplex term(1, 0);
    OphelieComplex sum = term;
    for (int k = 1; k <= 60; ++k)
    {
        term *= -z2 / (Real(4) * Real(k) * Real(k));
        sum += term;
        if (std::abs(term) < Real(1.0e-15) * (std::abs(sum) + Real(1)))
        {
            break;
        }
    }
    return sum;
}

/** Power-series J₁(z). */
inline OphelieComplex ophelieBesselJ1(OphelieComplex z)
{
    const OphelieComplex z2 = z * z;
    OphelieComplex term = z / Real(2); // k=0
    OphelieComplex sum = term;
    for (int k = 1; k <= 60; ++k)
    {
        term *= -z2 / (Real(4) * Real(k) * Real(k + 1));
        sum += term;
        if (std::abs(term) < Real(1.0e-15) * (std::abs(sum) + Real(1)))
        {
            break;
        }
    }
    return sum;
}

struct OphelieCylinderAnalyticSample
{
    Real r_m = 0.0;
    Real Hz_abs = 0.0;
    Real Bz_abs = 0.0;
    Real Hz_re = 0.0;
    Real Hz_im = 0.0;
    Real Jtheta_abs = 0.0;
    Real Jtheta_re = 0.0;
    Real Jtheta_im = 0.0;
    Real Etheta_abs = 0.0;
    Real Etheta_re = 0.0;
    Real Etheta_im = 0.0;
    Real Q_Wpm3 = 0.0;
};

struct OphelieCylinderAnalyticSolution
{
    Real R_m = 0.0;
    Real B0_T = 0.0;
    Real H0_Apm = 0.0;
    Real omega = 0.0;
    Real mu0 = 0.0;
    Real sigma = 0.0;
    Real delta_m = 0.0;
    Real delta_over_R = 0.0;
    Real P_per_length_Wpm = 0.0; // ∫ Q 2π r dr
    OphelieComplex alpha{};
    OphelieComplex j0_aR{};
};

inline OphelieCylinderAnalyticSolution makeOphelieCylinderAnalyticSolution(Real R, Real B0, Real omega, Real mu0,
                                                                            Real sigma)
{
    OphelieCylinderAnalyticSolution s;
    s.R_m = R;
    s.B0_T = B0;
    s.H0_Apm = B0 / (mu0 + TinyReal);
    s.omega = omega;
    s.mu0 = mu0;
    s.sigma = sigma;
    s.delta_m = ophelieSkinDepthM(omega, mu0, sigma);
    s.delta_over_R = s.delta_m / (R + TinyReal);
    s.alpha = ophelieCylinderAlpha(omega, mu0, sigma);
    s.j0_aR = ophelieBesselJ0(s.alpha * R);
    return s;
}

inline OphelieCylinderAnalyticSample sampleOphelieCylinderAnalytic(const OphelieCylinderAnalyticSolution &s, Real r)
{
    OphelieCylinderAnalyticSample o;
    o.r_m = r;
    if (r < Real(0) || r > s.R_m * Real(1.000001) || std::abs(s.j0_aR) < TinyReal)
    {
        return o;
    }
    const OphelieComplex j0 = ophelieBesselJ0(s.alpha * r);
    const OphelieComplex j1 = ophelieBesselJ1(s.alpha * r);
    const OphelieComplex hz = (s.H0_Apm / s.j0_aR) * j0;
    const OphelieComplex eth = (s.alpha / s.sigma) * (s.H0_Apm / s.j0_aR) * j1;
    const OphelieComplex jth = s.sigma * eth;
    o.Hz_re = hz.real();
    o.Hz_im = hz.imag();
    o.Hz_abs = std::abs(hz);
    o.Bz_abs = s.mu0 * o.Hz_abs;
    o.Jtheta_re = jth.real();
    o.Jtheta_im = jth.imag();
    o.Jtheta_abs = std::abs(jth);
    o.Etheta_re = eth.real();
    o.Etheta_im = eth.imag();
    o.Etheta_abs = std::abs(eth);
    o.Q_Wpm3 = Real(0.5) * s.sigma * o.Etheta_abs * o.Etheta_abs;
    return o;
}

/** Trapezoid ∫_0^R Q(r) 2π r dr. */
inline Real ophelieCylinderAnalyticPowerPerLength(const OphelieCylinderAnalyticSolution &s, int n_quad = 400)
{
    Real acc = 0.0;
    const Real dr = s.R_m / Real(n_quad);
    for (int i = 0; i <= n_quad; ++i)
    {
        const Real r = Real(i) * dr;
        const Real w = (i == 0 || i == n_quad) ? Real(0.5) : Real(1);
        const auto samp = sampleOphelieCylinderAnalytic(s, std::min(r, s.R_m));
        acc += w * samp.Q_Wpm3 * (Real(2) * Pi * r) * dr;
    }
    return acc;
}

inline void finalizeOphelieCylinderAnalyticPower(OphelieCylinderAnalyticSolution &s)
{
    s.P_per_length_Wpm = ophelieCylinderAnalyticPowerPerLength(s);
}

/** Host identity: |αR|→0 ⇒ E_θ → −i ω B₀ r / 2 (weak feedback, exp(+iωt)). */
inline bool verifyOphelieCylinderAnalyticWeakLimit(std::string &message)
{
    const Real R = Real(0.2);
    const Real B0 = Real(0.01);
    const Real mu0 = Real(1.2566370614e-6);
    const Real sigma = Real(16);
    // Very low frequency → deep penetration
    const Real omega = Real(2) * Pi * Real(50); // 50 Hz
    auto sol = makeOphelieCylinderAnalyticSolution(R, B0, omega, mu0, sigma);
    if (!(sol.delta_over_R > Real(5)))
    {
        message = "weak_limit_delta_over_R_too_small";
        return false;
    }
    // α must be (1−i)/δ under exp(+iωt) (doc 40 R09).
    if (!(sol.alpha.real() > TinyReal && sol.alpha.imag() < Real(0) &&
          std::abs(sol.alpha.real() + sol.alpha.imag()) < Real(1.0e-5) * sol.alpha.real()))
    {
        message = "alpha_phase_convention_failed";
        return false;
    }
    const Real r = Real(0.5) * R;
    const auto samp = sampleOphelieCylinderAnalytic(sol, r);
    const Real e_weak_mag = Real(0.5) * omega * B0 * r;
    const Real e_weak_im = -e_weak_mag; // E = −i ω B₀ r / 2
    const Real rel_abs = std::abs(samp.Etheta_abs - e_weak_mag) / (e_weak_mag + TinyReal);
    const Real rel_im = std::abs(samp.Etheta_im - e_weak_im) / (e_weak_mag + TinyReal);
    const Real rel_re = std::abs(samp.Etheta_re) / (e_weak_mag + TinyReal);
    if (!(rel_abs < Real(0.02) && rel_im < Real(0.02) && rel_re < Real(0.05)))
    {
        message = "weak_limit_E_phase_mismatch rel_abs=" + std::to_string(rel_abs) +
                  " rel_im=" + std::to_string(rel_im) + " rel_re=" + std::to_string(rel_re);
        return false;
    }
    // Bessel recurrence check at a moderate argument (float-safe step)
    const OphelieComplex z = OphelieComplex(Real(1.5), Real(-1.5));
    const OphelieComplex j1 = ophelieBesselJ1(z);
    const Real h = Real(1.0e-3);
    const OphelieComplex j0p = (ophelieBesselJ0(z + h) - ophelieBesselJ0(z - h)) / (Real(2) * h);
    const Real deriv_err = std::abs(j0p + j1) / (std::abs(j1) + Real(1));
    if (!(deriv_err < Real(5.0e-3)))
    {
        message = "bessel_J0_deriv_check_failed err=" + std::to_string(deriv_err);
        return false;
    }
    message = "cylinder_analytic_weak_limit_ok phase=exp(+iwt) E=-i*omega*B0*r/2";
    return true;
}

/** Doc198 E0: which analytic reference the mid-plane compare uses. */
enum class OphelieCylinderReferenceMode
{
    FullSkin = 0,           /**< infinite-cylinder self-consistent skin (Bessel) */
    PrescribedANoSelf = 1   /**< A=½B0(−y,x,0), φ=0, J=−iωσA (no reaction) */
};

inline const char *ophelieCylinderReferenceModeName(OphelieCylinderReferenceMode mode)
{
    return mode == OphelieCylinderReferenceMode::PrescribedANoSelf ? "prescribed_A_no_self" : "full_skin";
}

/**
 * Prescribed uniform Bz → A=½B0(−y,x,0); no self-induction.
 * E=−iωA ⇒ |J_θ|=σω B0 r/2, Q=½σ|E|² (exp(+iωt) contract).
 */
inline OphelieCylinderAnalyticSample sampleOphelieCylinderPrescribedANoSelf(Real r, Real R, Real B0, Real omega,
                                                                            Real sigma)
{
    OphelieCylinderAnalyticSample o;
    o.r_m = r;
    if (r < Real(0) || r > R * Real(1.000001))
    {
        return o;
    }
    const Real a_th = Real(0.5) * B0 * r;
    const Real e_abs = omega * a_th;
    o.Etheta_re = Real(0);
    o.Etheta_im = -e_abs;
    o.Etheta_abs = e_abs;
    o.Jtheta_re = Real(0);
    o.Jtheta_im = -sigma * e_abs;
    o.Jtheta_abs = sigma * e_abs;
    o.Q_Wpm3 = Real(0.5) * sigma * e_abs * e_abs;
    o.Hz_re = B0 / (Real(1.2566370614e-6)); // unused for J compare; leave rough
    o.Hz_im = Real(0);
    o.Hz_abs = o.Hz_re;
    o.Bz_abs = B0;
    return o;
}

inline Real ophelieCylinderPrescribedANoSelfPowerPerLength(Real R, Real B0, Real omega, Real sigma)
{
    // ∫_0^R ½ σ (ω B0 r / 2)^2 2π r dr = (π/16) σ ω² B0² R^4
    return (Pi / Real(16)) * sigma * omega * omega * B0 * B0 * R * R * R * R;
}

inline OphelieCylinderAnalyticSample sampleOphelieCylinderByReferenceMode(
    OphelieCylinderReferenceMode mode, const OphelieCylinderAnalyticSolution &full_skin_sol, Real r, Real R,
    Real B0, Real omega, Real sigma)
{
    if (mode == OphelieCylinderReferenceMode::PrescribedANoSelf)
    {
        return sampleOphelieCylinderPrescribedANoSelf(r, R, B0, omega, sigma);
    }
    return sampleOphelieCylinderAnalytic(full_skin_sol, r);
}

struct OphelieCylinderMidplaneCompareReport
{
    Real delta_over_R = 0.0;
    Real P_sph_W = 0.0;
    Real P_analytic_finite_W = 0.0; // P' * height
    Real rel_P = 0.0;
    Real corr_abs_Jtheta = 0.0;
    Real corr_Q_shape = 0.0;
    Real mean_rel_err_abs_J = 0.0;
    Real median_rel_err_abs_J = 0.0;
    /** Doc197: volume-weighted L2 of |J_θ| on mid-plane band (≠ bin mean_rel_J). */
    Real l2_rel_abs_J_mid = 0.0;
    /** Same L2 restricted to strong field |J_an| > strong_floor * max|J_an|_mid. */
    Real l2_rel_abs_J_mid_strong = 0.0;
    Real mean_rel_err_abs_J_end = 0.0; /**< end slabs; finite-height vs infinite ref */
    Real q_outer_over_center_sph = 0.0;
    Real q_outer_over_center_analytic = 0.0;
    size_t n_bins = 0;
    size_t n_midplane = 0;
    size_t n_endband = 0;
    size_t n_mid_strong = 0;
    OphelieCylinderReferenceMode reference_mode = OphelieCylinderReferenceMode::FullSkin;
};

inline Real ophelieMedianInplace(std::vector<Real> &xs)
{
    if (xs.empty())
    {
        return Real(0);
    }
    const size_t m = xs.size() / 2;
    std::nth_element(xs.begin(), xs.begin() + static_cast<std::ptrdiff_t>(m), xs.end());
    return xs[m];
}

/**
 * Volume-weighted mid-plane (|z-zc|<z_half_band) radial bins vs analytic.
 * J_θ from SPH: azimuthal component of j_edge_recon (real/imag).
 */
/**
 * Per-particle mid-plane samples for lattice-anisotropy audit (Doc197).
 * Columns: x,y,z,r,theta_deg,vol,J_sph,J_an,signed_rel,abs_rel
 */
inline bool hostOphelieCylinderDumpMidplaneParticleCsv(
    const Vecd *pos, const Real *vol, const Vecd *j_re, const Vecd *j_im, size_t n, const Vecd &axis_origin,
    Real radius_m, Real z_half_band, Real B0, Real omega, Real mu0, Real sigma, const std::string &csv_path,
    OphelieCylinderReferenceMode reference_mode = OphelieCylinderReferenceMode::FullSkin)
{
    if (csv_path.empty())
    {
        return false;
    }
    auto sol = makeOphelieCylinderAnalyticSolution(radius_m, B0, omega, mu0, sigma);
    finalizeOphelieCylinderAnalyticPower(sol);
    const Real zc = axis_origin[2];
    std::ofstream ofs(csv_path);
    if (!ofs)
    {
        std::cerr << "[ophelie][cyl-analytic] failed to write particle CSV: " << csv_path << std::endl;
        return false;
    }
    ofs << "x,y,z,r,theta_deg,vol,J_sph,J_an,signed_rel,abs_rel\n";
    size_t n_out = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        if (r >= radius_m || std::abs(pos[i][2] - zc) > z_half_band)
        {
            continue;
        }
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const Real jt_re = inv_r * (-j_re[i][0] * dy + j_re[i][1] * dx);
        const Real jt_im = inv_r * (-j_im[i][0] * dy + j_im[i][1] * dx);
        const Real jt = std::hypot(jt_re, jt_im);
        const auto samp =
            sampleOphelieCylinderByReferenceMode(reference_mode, sol, r, radius_m, B0, omega, sigma);
        const Real signed_rel =
            (samp.Jtheta_abs > Real(1.0e-9)) ? (jt - samp.Jtheta_abs) / samp.Jtheta_abs : Real(0);
        const Real th_deg = static_cast<Real>(std::atan2(dy, dx) * Real(180) / Pi);
        ofs << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << "," << th_deg << "," << vol[i]
            << "," << jt << "," << samp.Jtheta_abs << "," << signed_rel << "," << std::abs(signed_rel) << "\n";
        ++n_out;
    }
    std::cout << "[ophelie][cyl-analytic] midplane_particle_csv=" << csv_path << " n=" << n_out << std::endl;
    return n_out > 0;
}

/**
 * Same mid-plane bins as hostOphelieCylinderMidplaneCompare; optional CSV for Doc197 root-cause.
 * Columns: bin,r_mid_m,w_vol,J_sph,J_an,rel_err_J,Q_sph,Q_an
 */
inline bool hostOphelieCylinderDumpMidplaneRadialCsv(
    const Vecd *pos, const Real *vol, const Vecd *j_re, const Vecd *j_im, const Real *q, size_t n,
    const Vecd &axis_origin, Real radius_m, Real z_half_band, Real B0, Real omega, Real mu0, Real sigma,
    const std::string &csv_path, int n_bins = 12,
    OphelieCylinderReferenceMode reference_mode = OphelieCylinderReferenceMode::FullSkin)
{
    if (csv_path.empty())
    {
        return false;
    }
    auto sol = makeOphelieCylinderAnalyticSolution(radius_m, B0, omega, mu0, sigma);
    finalizeOphelieCylinderAnalyticPower(sol);
    std::vector<double> w_bin(static_cast<size_t>(n_bins), 0.0);
    std::vector<double> j_bin(static_cast<size_t>(n_bins), 0.0);
    std::vector<double> q_bin(static_cast<size_t>(n_bins), 0.0);
    const Real zc = axis_origin[2];
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        if (r >= radius_m || std::abs(pos[i][2] - zc) > z_half_band)
        {
            continue;
        }
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const Real jt_re = inv_r * (-j_re[i][0] * dy + j_re[i][1] * dx);
        const Real jt_im = inv_r * (-j_im[i][0] * dy + j_im[i][1] * dx);
        const Real jt = std::hypot(jt_re, jt_im);
        const int b = std::min(n_bins - 1, static_cast<int>(r / radius_m * Real(n_bins)));
        w_bin[static_cast<size_t>(b)] += static_cast<double>(vol[i]);
        j_bin[static_cast<size_t>(b)] += static_cast<double>(vol[i] * jt);
        q_bin[static_cast<size_t>(b)] += static_cast<double>(vol[i] * q[i]);
    }
    std::ofstream ofs(csv_path);
    if (!ofs)
    {
        std::cerr << "[ophelie][cyl-analytic] failed to write radial CSV: " << csv_path << std::endl;
        return false;
    }
    ofs << "bin,r_mid_m,w_vol,J_sph,J_an,rel_err_J,Q_sph,Q_an\n";
    for (int b = 0; b < n_bins; ++b)
    {
        if (w_bin[static_cast<size_t>(b)] <= TinyReal)
        {
            continue;
        }
        const Real r_mid = radius_m * (Real(b) + Real(0.5)) / Real(n_bins);
        const auto samp =
            sampleOphelieCylinderByReferenceMode(reference_mode, sol, r_mid, radius_m, B0, omega, sigma);
        const Real js = static_cast<Real>(j_bin[static_cast<size_t>(b)] / w_bin[static_cast<size_t>(b)]);
        const Real qs = static_cast<Real>(q_bin[static_cast<size_t>(b)] / w_bin[static_cast<size_t>(b)]);
        const Real rel =
            (samp.Jtheta_abs > Real(1.0e-9)) ? std::abs(js - samp.Jtheta_abs) / samp.Jtheta_abs : Real(0);
        ofs << b << "," << r_mid << "," << w_bin[static_cast<size_t>(b)] << "," << js << "," << samp.Jtheta_abs
            << "," << rel << "," << qs << "," << samp.Q_Wpm3 << "\n";
    }
    std::cout << "[ophelie][cyl-analytic] radial_csv=" << csv_path
              << " reference_mode=" << ophelieCylinderReferenceModeName(reference_mode) << std::endl;
    return true;
}

inline OphelieCylinderMidplaneCompareReport hostOphelieCylinderMidplaneCompare(
    const Vecd *pos, const Real *vol, const Vecd *j_re, const Vecd *j_im, const Real *q, size_t n,
    const Vecd &axis_origin, Real radius_m, Real half_height, Real z_half_band, Real B0, Real omega, Real mu0,
    Real sigma, int n_bins = 12,
    OphelieCylinderReferenceMode reference_mode = OphelieCylinderReferenceMode::FullSkin)
{
    OphelieCylinderMidplaneCompareReport rep;
    rep.reference_mode = reference_mode;
    auto sol = makeOphelieCylinderAnalyticSolution(radius_m, B0, omega, mu0, sigma);
    finalizeOphelieCylinderAnalyticPower(sol);
    rep.delta_over_R = sol.delta_over_R;
    const Real height = Real(2) * half_height;
    if (reference_mode == OphelieCylinderReferenceMode::PrescribedANoSelf)
    {
        // Weak-skin limit: δ/R → ∞; still report skin sol.delta for context.
        rep.P_analytic_finite_W =
            ophelieCylinderPrescribedANoSelfPowerPerLength(radius_m, B0, omega, sigma) * height;
    }
    else
    {
        rep.P_analytic_finite_W = sol.P_per_length_Wpm * height;
    }

    std::vector<double> w_bin(static_cast<size_t>(n_bins), 0.0);
    std::vector<double> j_bin(static_cast<size_t>(n_bins), 0.0);
    std::vector<double> q_bin(static_cast<size_t>(n_bins), 0.0);
    double p_tot = 0.0;
    double q_o = 0.0, v_o = 0.0, q_c = 0.0, v_c = 0.0;
    double l2_num = 0.0, l2_den = 0.0;
    double l2s_num = 0.0, l2s_den = 0.0;
    double end_rel_sum = 0.0;
    size_t n_end_rel = 0;
    Real j_an_max_mid = Real(0);
    const Real r_o = Real(0.75) * radius_m;
    const Real r_c = Real(0.25) * radius_m;
    const Real zc = axis_origin[2];
    // End slabs: outer 20% of half-height each side (finite-cylinder pollution vs infinite ref).
    const Real z_end_inner = Real(0.80) * half_height;

    struct MidSamp
    {
        Real vi, jt, j_an;
    };
    std::vector<MidSamp> mid_samps;
    mid_samps.reserve(n / 4);

    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real vi = vol[i];
        p_tot += static_cast<double>(q[i]) * static_cast<double>(vi);
        if (r > r_o)
        {
            q_o += q[i] * vi;
            v_o += vi;
        }
        if (r < r_c)
        {
            q_c += q[i] * vi;
            v_c += vi;
        }
        if (r >= radius_m)
        {
            continue;
        }
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const Real jt_re = inv_r * (-j_re[i][0] * dy + j_re[i][1] * dx);
        const Real jt_im = inv_r * (-j_im[i][0] * dy + j_im[i][1] * dx);
        const Real jt = std::hypot(jt_re, jt_im);
        const auto samp =
            sampleOphelieCylinderByReferenceMode(reference_mode, sol, r, radius_m, B0, omega, sigma);
        const Real abs_z = std::abs(pos[i][2] - zc);
        if (abs_z > z_end_inner)
        {
            ++rep.n_endband;
            if (samp.Jtheta_abs > Real(1.0e-9))
            {
                end_rel_sum += static_cast<double>(std::abs(jt - samp.Jtheta_abs) / samp.Jtheta_abs);
                ++n_end_rel;
            }
        }
        if (abs_z > z_half_band)
        {
            continue;
        }
        ++rep.n_midplane;
        mid_samps.push_back(MidSamp{vi, jt, samp.Jtheta_abs});
        j_an_max_mid = std::max(j_an_max_mid, samp.Jtheta_abs);
        const int b = std::min(n_bins - 1, static_cast<int>(r / radius_m * Real(n_bins)));
        w_bin[static_cast<size_t>(b)] += vi;
        j_bin[static_cast<size_t>(b)] += vi * jt;
        q_bin[static_cast<size_t>(b)] += vi * q[i];
    }
    const Real strong_floor = Real(0.10) * j_an_max_mid; // pre-declared: ignore weak 10% of peak
    for (const MidSamp &m : mid_samps)
    {
        const double d = static_cast<double>(m.jt - m.j_an);
        const double ja = static_cast<double>(m.j_an);
        const double w = static_cast<double>(m.vi);
        l2_num += w * d * d;
        l2_den += w * ja * ja;
        if (m.j_an > strong_floor)
        {
            ++rep.n_mid_strong;
            l2s_num += w * d * d;
            l2s_den += w * ja * ja;
        }
    }
    rep.l2_rel_abs_J_mid =
        (l2_den > TinyReal) ? static_cast<Real>(std::sqrt(l2_num / l2_den)) : Real(0);
    rep.l2_rel_abs_J_mid_strong =
        (l2s_den > TinyReal) ? static_cast<Real>(std::sqrt(l2s_num / l2s_den)) : Real(0);
    if (n_end_rel > 0)
    {
        rep.mean_rel_err_abs_J_end = static_cast<Real>(end_rel_sum / static_cast<double>(n_end_rel));
    }
    rep.P_sph_W = static_cast<Real>(p_tot);
    rep.rel_P = std::abs(rep.P_sph_W - rep.P_analytic_finite_W) / (std::abs(rep.P_analytic_finite_W) + TinyReal);
    rep.q_outer_over_center_sph = (v_o > TinyReal && v_c > TinyReal) ? static_cast<Real>((q_o / v_o) / (q_c / v_c + TinyReal))
                                                                    : Real(0);

    // Analytic shell means with area weight r dr
    {
        double qo = 0, vo = 0, qc = 0, vc = 0;
        const int nq = 200;
        for (int i = 0; i < nq; ++i)
        {
            const Real r0 = radius_m * Real(i) / Real(nq);
            const Real r1 = radius_m * Real(i + 1) / Real(nq);
            const Real rm = Real(0.5) * (r0 + r1);
            const Real area = Pi * (r1 * r1 - r0 * r0);
            const auto samp =
                sampleOphelieCylinderByReferenceMode(reference_mode, sol, rm, radius_m, B0, omega, sigma);
            if (rm > r_o)
            {
                qo += samp.Q_Wpm3 * area;
                vo += area;
            }
            if (rm < r_c)
            {
                qc += samp.Q_Wpm3 * area;
                vc += area;
            }
        }
        rep.q_outer_over_center_analytic =
            (vo > TinyReal && vc > TinyReal) ? static_cast<Real>((qo / vo) / (qc / vc + TinyReal)) : Real(0);
    }

    std::vector<Real> js, ja, qs, qa, errj;
    for (int b = 0; b < n_bins; ++b)
    {
        if (w_bin[static_cast<size_t>(b)] <= TinyReal)
        {
            continue;
        }
        const Real r_mid = radius_m * (Real(b) + Real(0.5)) / Real(n_bins);
        const auto samp =
            sampleOphelieCylinderByReferenceMode(reference_mode, sol, r_mid, radius_m, B0, omega, sigma);
        js.push_back(static_cast<Real>(j_bin[static_cast<size_t>(b)] / w_bin[static_cast<size_t>(b)]));
        ja.push_back(samp.Jtheta_abs);
        qs.push_back(static_cast<Real>(q_bin[static_cast<size_t>(b)] / w_bin[static_cast<size_t>(b)]));
        qa.push_back(samp.Q_Wpm3);
        if (samp.Jtheta_abs > Real(1.0e-9))
        {
            errj.push_back(std::abs(js.back() - samp.Jtheta_abs) / samp.Jtheta_abs);
        }
    }
    rep.n_bins = js.size();

    auto pearson = [](const std::vector<Real> &a, const std::vector<Real> &b) -> Real {
        if (a.size() < 2)
        {
            return Real(0);
        }
        double ma = 0, mb = 0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            ma += a[i];
            mb += b[i];
        }
        ma /= static_cast<double>(a.size());
        mb /= static_cast<double>(b.size());
        double num = 0, da = 0, db = 0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            const double xa = a[i] - ma;
            const double xb = b[i] - mb;
            num += xa * xb;
            da += xa * xa;
            db += xb * xb;
        }
        const double den = std::sqrt(da * db);
        return den > TinyReal ? static_cast<Real>(num / den) : Real(0);
    };
    rep.corr_abs_Jtheta = pearson(js, ja);
    // Q shape: normalize by mean
    if (!qs.empty())
    {
        Real ms = 0, ma = 0;
        for (size_t i = 0; i < qs.size(); ++i)
        {
            ms += qs[i];
            ma += qa[i];
        }
        ms /= Real(qs.size());
        ma /= Real(qa.size());
        std::vector<Real> qsn, qan;
        for (size_t i = 0; i < qs.size(); ++i)
        {
            qsn.push_back(qs[i] / (ms + TinyReal));
            qan.push_back(qa[i] / (ma + TinyReal));
        }
        rep.corr_Q_shape = pearson(qsn, qan);
    }
    if (!errj.empty())
    {
        double sum = 0;
        for (Real e : errj)
        {
            sum += e;
        }
        rep.mean_rel_err_abs_J = static_cast<Real>(sum / static_cast<double>(errj.size()));
        rep.median_rel_err_abs_J = ophelieMedianInplace(errj);
    }
    return rep;
}

inline void logOphelieCylinderMidplaneCompare(const OphelieCylinderMidplaneCompareReport &r, const char *tag)
{
    std::cout << "[ophelie][cyl-analytic] " << tag
              << " reference_mode=" << ophelieCylinderReferenceModeName(r.reference_mode)
              << " delta/R=" << r.delta_over_R << " n_bins=" << r.n_bins
              << " n_mid=" << r.n_midplane << " n_end=" << r.n_endband << " n_mid_strong=" << r.n_mid_strong
              << " corr_J=" << r.corr_abs_Jtheta << " corr_Q=" << r.corr_Q_shape
              << " mean/med_rel_J=" << r.mean_rel_err_abs_J << "/" << r.median_rel_err_abs_J
              << " L2_rel_J_mid=" << r.l2_rel_abs_J_mid << " L2_rel_J_mid_strong=" << r.l2_rel_abs_J_mid_strong
              << " mean_rel_J_end=" << r.mean_rel_err_abs_J_end
              << " P_sph=" << r.P_sph_W << " P_analytic_H=" << r.P_analytic_finite_W << " rel_P=" << r.rel_P
              << " Q_o/c_sph=" << r.q_outer_over_center_sph << " Q_o/c_an=" << r.q_outer_over_center_analytic
              << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_SELF_INDUCTION_H
