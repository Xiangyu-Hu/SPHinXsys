/**
 * @file electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h
 * @brief Doc199 E2: A=∇χ, φ=-iωχ ⇒ E=J=0 manufactured cancel (exp(+iωt)).
 *
 * Phase 1: inject analytic φ, recon only (no phi solve).
 * Gate: vol L2 |J| / (σ ω |∇χ|) ≤ 1%. Do NOT divide by J_ref=0.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC199_E2_GRAD_A_ZERO_J_H
#define ELECTROMAGNETIC_OPHELIE_DOC199_E2_GRAD_A_ZERO_J_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OphelieDoc199E2ChiKind
{
    LinearX = 0,
    QuadraticZ = 1
};

inline const char *ophelieDoc199E2ChiKindName(OphelieDoc199E2ChiKind k)
{
    return k == OphelieDoc199E2ChiKind::QuadraticZ ? "quadratic_z" : "linear_x";
}

struct OphelieDoc199E2ChiField
{
    OphelieDoc199E2ChiKind kind = OphelieDoc199E2ChiKind::LinearX;
    Real g0 = 0.0; // |∇χ| for linear; also α for quadratic (A_z = α z_rel)
    Vecd center = Vecd::Zero();
};

/** Characteristic |∇χ| matched to cylinder rim |A_rot|=B0 R/2. */
inline Real ophelieDoc199E2GradChiScale(Real b0, Real radius)
{
    return Real(0.5) * b0 * radius;
}

inline Real ophelieDoc199E2ChiAt(const Vecd &pos, const OphelieDoc199E2ChiField &chi)
{
    const Real xr = pos[0] - chi.center[0];
    const Real yr = pos[1] - chi.center[1];
    const Real zr = pos[2] - chi.center[2];
    (void)yr;
    if (chi.kind == OphelieDoc199E2ChiKind::QuadraticZ)
        return Real(0.5) * chi.g0 * zr * zr;
    return chi.g0 * xr;
}

inline Vecd ophelieDoc199E2GradChiAt(const Vecd &pos, const OphelieDoc199E2ChiField &chi)
{
    const Real zr = pos[2] - chi.center[2];
    if (chi.kind == OphelieDoc199E2ChiKind::QuadraticZ)
        return Vecd(0.0, 0.0, chi.g0 * zr);
    return Vecd(chi.g0, 0.0, 0.0);
}

struct OphelieDoc199E2Region
{
    size_t n = 0;
    double vol = 0.0;
    double j_num = 0.0;   // Σ V |J|²
    double scale = 0.0;   // Σ V (σω|∇χ|)²
    double q_abs = 0.0;   // Σ V |Q|
    double edge_drop_sq = 0.0; // optional diagnostic sum
};

inline Real ophelieDoc199E2RelJL2(const OphelieDoc199E2Region &s)
{
    return (s.scale > TinyReal) ? static_cast<Real>(std::sqrt(s.j_num / s.scale)) : Real(0);
}

struct OphelieDoc199E2ZeroJReport
{
    std::string chi_name = "unknown";
    int phase = 1;
    Real j_rel_l2_all = 0.0;
    Real j_rel_l2_interior = 0.0;
    Real j_rel_l2_rim = 0.0;
    Real q_vol_abs_mean = 0.0;
    Real phi_imag_rms = 0.0;
    Real phi_real_rms = 0.0;
    /** Doc200 T1: ‖φ_i + ωχ − mean‖_vol / (ω χ_rms); gauge constant removed. */
    Real phi_imag_vs_minus_omega_chi_rel = 0.0;
    Real phi_real_vs_zero_rel = 0.0; // ‖φ_r‖_vol / (ω χ_rms)
    Real a_vs_gradchi_rel = 0.0;
    Real omega_gradchi_l2 = 0.0;
    bool gate_j_all = false;
    bool gate_j_interior = false;
    bool gate_phi_match = false;
    OphelieDoc199E2Region all;
    OphelieDoc199E2Region interior;
    OphelieDoc199E2Region rim;
};

inline const char *ophelieDoc199E2RegionTag(const Vecd &pos, const Vecd &center, Real radius, Real half_height)
{
    const Real dx = pos[0] - center[0];
    const Real dy = pos[1] - center[1];
    const Real r = std::sqrt(dx * dx + dy * dy);
    const Real az = std::abs(pos[2] - center[2]);
    if (az > Real(0.80) * half_height)
        return "end";
    if (r >= Real(0.75) * radius)
        return "rim";
    if (r < Real(0.70) * radius && az < Real(0.50) * half_height)
        return "interior";
    return "bulk";
}

inline void ophelieDoc199E2Acc(OphelieDoc199E2Region &s, Real v, Real j2, Real scale2, Real qabs)
{
    const double vv = static_cast<double>(v);
    s.vol += vv;
    s.j_num += vv * static_cast<double>(j2);
    s.scale += vv * static_cast<double>(scale2);
    s.q_abs += vv * static_cast<double>(qabs);
    ++s.n;
}

/**
 * Host audit after edge recon.
 * phase=1: injected analytic φ; phase=2 (Doc200 T1): solved φ vs −ωχ (gauge removed).
 * Expects j_edge_recon_* and joule_heat_edge_recon_complex already filled.
 * J* scale = σ ω |∇χ|_char (never shrink by residual J).
 */
inline OphelieDoc199E2ZeroJReport hostOphelieDoc199E2ZeroJAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieDoc199E2ChiField &chi, Real omega,
    Real radius, Real half_height, const std::string &particle_csv_path = {}, int phase = 1)
{
    OphelieDoc199E2ZeroJReport rep;
    rep.chi_name = ophelieDoc199E2ChiKindName(chi.kind);
    rep.phase = phase;
    // Non-vanishing characteristic |∇χ| (Doc199: do not normalize by J_ref=0 or midplane zero grad).
    const Real g_char = (chi.kind == OphelieDoc199E2ChiKind::QuadraticZ)
                            ? std::abs(chi.g0) * half_height
                            : std::abs(chi.g0);

    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Vecd *a = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);

    double phi2 = 0.0, phi_r2 = 0.0, vol_tot = 0.0;
    double a_err = 0.0, a_den = 0.0, og2 = 0.0;
    double psi_sum = 0.0, chi2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd g = ophelieDoc199E2GradChiAt(pos[i], chi);
        const Real j2 = jr[i].squaredNorm() + ji[i].squaredNorm();
        const Real scale = std::max(sigma[i] * omega * g_char, Real(TinyReal));
        const Real scale2 = scale * scale;
        ophelieDoc199E2Acc(rep.all, vol[i], j2, scale2, std::abs(q[i]));
        const char *tag = ophelieDoc199E2RegionTag(pos[i], chi.center, radius, half_height);
        if (std::strcmp(tag, "interior") == 0)
            ophelieDoc199E2Acc(rep.interior, vol[i], j2, scale2, std::abs(q[i]));
        if (std::strcmp(tag, "rim") == 0)
            ophelieDoc199E2Acc(rep.rim, vol[i], j2, scale2, std::abs(q[i]));

        const double vv = static_cast<double>(vol[i]);
        vol_tot += vv;
        phi2 += vv * static_cast<double>(phi[i] * phi[i]);
        phi_r2 += vv * static_cast<double>(phi_r[i] * phi_r[i]);
        const Real chi_i = ophelieDoc199E2ChiAt(pos[i], chi);
        psi_sum += vv * static_cast<double>(phi[i] + omega * chi_i);
        chi2 += vv * static_cast<double>(chi_i * chi_i);
        const Vecd da = a[i] - g;
        a_err += vv * static_cast<double>(da.squaredNorm());
        a_den += vv * static_cast<double>(std::max(g.squaredNorm(), g_char * g_char));
        og2 += vv * static_cast<double>((omega * g_char) * (omega * g_char));
    }

    const double mean_psi = (vol_tot > TinyReal) ? (psi_sum / vol_tot) : 0.0;
    double psi_err2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double vv = static_cast<double>(vol[i]);
        const Real chi_i = ophelieDoc199E2ChiAt(pos[i], chi);
        const double dpsi = static_cast<double>(phi[i] + omega * chi_i) - mean_psi;
        psi_err2 += vv * dpsi * dpsi;
    }
    const Real chi_rms = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(chi2 / vol_tot)) : Real(0);
    const Real phi_scale = std::max(omega * std::max(chi_rms, Real(TinyReal)), Real(TinyReal));

    rep.j_rel_l2_all = ophelieDoc199E2RelJL2(rep.all);
    rep.j_rel_l2_interior = ophelieDoc199E2RelJL2(rep.interior);
    rep.j_rel_l2_rim = ophelieDoc199E2RelJL2(rep.rim);
    rep.q_vol_abs_mean = (rep.all.vol > TinyReal) ? static_cast<Real>(rep.all.q_abs / rep.all.vol) : Real(0);
    rep.phi_imag_rms = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(phi2 / vol_tot)) : Real(0);
    rep.phi_real_rms = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(phi_r2 / vol_tot)) : Real(0);
    rep.phi_imag_vs_minus_omega_chi_rel =
        (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(psi_err2 / vol_tot) / phi_scale) : Real(0);
    rep.phi_real_vs_zero_rel = rep.phi_real_rms / phi_scale;
    rep.a_vs_gradchi_rel = (a_den > TinyReal) ? static_cast<Real>(std::sqrt(a_err / a_den)) : Real(0);
    rep.omega_gradchi_l2 = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(og2 / vol_tot)) : Real(0);
    // Contract: J_rel all ≤1%; T1 also φ match ≤1% (gauge-free).
    rep.gate_j_all = rep.j_rel_l2_all <= Real(0.01);
    rep.gate_j_interior = rep.j_rel_l2_interior <= Real(0.01);
    rep.gate_phi_match = (rep.phi_imag_vs_minus_omega_chi_rel <= Real(0.01)) &&
                         (rep.phi_real_vs_zero_rel <= Real(0.01));

    if (!particle_csv_path.empty())
    {
        std::ofstream csv(particle_csv_path);
        csv << "i,x,y,z,r,region,Vol,chi,phi_imag,phi_real,phi_i+w*chi,|gradchi|,|A|,|J|,Q,"
               "Jx_re,Jy_re,Jz_re,Jx_im,Jy_im,Jz_im\n";
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - chi.center[0];
            const Real dy = pos[i][1] - chi.center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Vecd g = ophelieDoc199E2GradChiAt(pos[i], chi);
            const Real chi_i = ophelieDoc199E2ChiAt(pos[i], chi);
            const Real jn = std::sqrt(jr[i].squaredNorm() + ji[i].squaredNorm());
            csv << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << ","
                << ophelieDoc199E2RegionTag(pos[i], chi.center, radius, half_height) << "," << vol[i] << "," << chi_i
                << "," << phi[i] << "," << phi_r[i] << "," << (phi[i] + omega * chi_i) << "," << g.norm() << ","
                << a[i].norm() << "," << jn << "," << q[i] << "," << jr[i][0] << "," << jr[i][1] << "," << jr[i][2]
                << "," << ji[i][0] << "," << ji[i][1] << "," << ji[i][2] << "\n";
        }
    }

    std::cout << "[ophelie][doc199][e2] phase=" << rep.phase << " chi=" << rep.chi_name
              << " a_vs_gradchi_rel=" << rep.a_vs_gradchi_rel << " omega|gradchi|_L2=" << rep.omega_gradchi_l2
              << " phi_imag_rms=" << rep.phi_imag_rms << " phi_real_rms=" << rep.phi_real_rms
              << " phi_i_vs_-wchi_rel=" << rep.phi_imag_vs_minus_omega_chi_rel
              << " phi_r_vs_0_rel=" << rep.phi_real_vs_zero_rel << " J_rel_L2 all=" << rep.j_rel_l2_all
              << " interior=" << rep.j_rel_l2_interior << " rim=" << rep.j_rel_l2_rim
              << " Q_abs_mean=" << rep.q_vol_abs_mean << std::endl;
    std::cout << "[ophelie][doc199][e2] gates J_rel_all<=1%=" << (rep.gate_j_all ? 1 : 0)
              << " J_rel_interior<=1%=" << (rep.gate_j_interior ? 1 : 0)
              << " phi_match<=1%=" << (rep.gate_phi_match ? 1 : 0)
              << " note: scale=sigma*omega*|grad chi|; phase1=inject, phase2=solve+Neumann" << std::endl;
    if (!particle_csv_path.empty())
        std::cout << "[ophelie][doc199][e2] particle_csv=" << particle_csv_path << std::endl;
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC199_E2_GRAD_A_ZERO_J_H
