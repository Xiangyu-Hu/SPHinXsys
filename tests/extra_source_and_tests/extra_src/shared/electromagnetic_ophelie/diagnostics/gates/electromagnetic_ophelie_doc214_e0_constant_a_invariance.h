/**
 * @file electromagnetic_ophelie_doc214_e0_constant_a_invariance.h
 * @brief Doc214 E0: same production profile; Arot / const A0 / Arot+A0 invariance.
 *
 * Gates (project engineering, not literature):
 *   E0-R vs E0-RC: vol-weighted complex J L2 ≤1%, Q-field L2 ≤2%
 *   E0-C: Jrms / (σ ω |A0|) ≤1% (all / interior / side=rim / end)
 * Phi transform check removes only a global constant C (no linear fit).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC214_E0_CONSTANT_A_INVARIANCE_H
#define ELECTROMAGNETIC_OPHELIE_DOC214_E0_CONSTANT_A_INVARIANCE_H

#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h"
#include "electromagnetic_ophelie_doc199_e3_gauge_transform.h"
#include "electromagnetic_ophelie_field_names.h"

#include <cmath>
#include <cstring>
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

struct OphelieDoc214E0Snapshot
{
    std::string tag;
    std::vector<Vecd> jr, ji, ar, ai;
    std::vector<Real> pr, pi, q, vol;
    Real phi_eq = 0.0;
    Real phi_solver_rel = 0.0;
};

struct OphelieDoc214E0RegionStats
{
    size_t n = 0;
    Real vol = 0.0;
    Real j_rms = 0.0;     // sqrt(Σ V |J|² / Σ V)
    Real j_rel_char = 0.0; // j_rms / Jchar
    Real q_l2 = 0.0;       // sqrt(Σ V Q² / Σ V)
};

struct OphelieDoc214E0Report
{
    Real a0_mag = 0.0;
    Real jchar = 0.0; // σ ω |A0|
    Real omega = 0.0;
    Real j_rel_rc_vs_r_all = 0.0;
    Real j_rel_rc_vs_r_interior = 0.0;
    Real j_rel_rc_vs_r_side = 0.0;
    Real j_rel_rc_vs_r_end = 0.0;
    Real q_rel_rc_vs_r_all = 0.0;
    Real q_rel_rc_vs_r_interior = 0.0;
    Real q_rel_rc_vs_r_side = 0.0;
    Real q_rel_rc_vs_r_end = 0.0;
    Real p_r = 0.0;
    Real p_rc = 0.0;
    Real p_c = 0.0;
    Real p_rel_rc_vs_r = 0.0;
    /** After removing global mean C only: ‖φ_RC − (φ_R − i ω A0·Δx) − C‖ / scale */
    Real phi_transform_rel = 0.0;
    Real phi_raw_rel_no_align = 0.0;
    OphelieDoc214E0RegionStats c_all, c_interior, c_side, c_end;
    bool gate_j_rc = false;
    bool gate_q_rc = false;
    bool gate_c = false;
    bool passed = false;
};

inline void ophelieDoc214E0FillRegionMasks(const Vecd *pos, size_t n, const Vecd &center, Real radius,
                                          Real half_height, std::vector<char> &all, std::vector<char> &interior,
                                          std::vector<char> &side, std::vector<char> &end)
{
    all.assign(n, 1);
    interior.assign(n, 0);
    side.assign(n, 0);
    end.assign(n, 0);
    for (size_t i = 0; i < n; ++i)
    {
        const char *tag = ophelieDoc199E2RegionTag(pos[i], center, radius, half_height);
        if (std::strcmp(tag, "interior") == 0)
            interior[i] = 1;
        else if (std::strcmp(tag, "rim") == 0)
            side[i] = 1; // lateral / side wall band
        else if (std::strcmp(tag, "end") == 0)
            end[i] = 1;
    }
}

inline Real ophelieDoc214E0QFieldRel(const std::vector<Real> &vol, const std::vector<Real> &q0,
                                    const std::vector<Real> &q1, const std::vector<char> &mask)
{
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < vol.size(); ++i)
    {
        if (!mask[i])
            continue;
        const double v = static_cast<double>(vol[i]);
        const double d = static_cast<double>(q1[i] - q0[i]);
        num += v * d * d;
        den += v * static_cast<double>(q0[i] * q0[i]);
    }
    return (den > TinyReal) ? static_cast<Real>(std::sqrt(num / den)) : Real(0);
}

inline OphelieDoc214E0RegionStats ophelieDoc214E0JRegion(const std::vector<Real> &vol, const std::vector<Vecd> &jr,
                                                        const std::vector<Vecd> &ji, const std::vector<Real> &q,
                                                        const std::vector<char> &mask, Real jchar)
{
    OphelieDoc214E0RegionStats s;
    double j2 = 0.0, q2 = 0.0, vsum = 0.0;
    for (size_t i = 0; i < vol.size(); ++i)
    {
        if (!mask[i])
            continue;
        const double v = static_cast<double>(vol[i]);
        vsum += v;
        j2 += v * static_cast<double>(jr[i].squaredNorm() + ji[i].squaredNorm());
        q2 += v * static_cast<double>(q[i] * q[i]);
        ++s.n;
    }
    s.vol = static_cast<Real>(vsum);
    s.j_rms = (vsum > TinyReal) ? static_cast<Real>(std::sqrt(j2 / vsum)) : Real(0);
    s.q_l2 = (vsum > TinyReal) ? static_cast<Real>(std::sqrt(q2 / vsum)) : Real(0);
    s.j_rel_char = (jchar > TinyReal) ? s.j_rms / jchar : Real(0);
    return s;
}

inline void hostOphelieDoc214E0CaptureSnapshot(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                              const std::string &tag, Real phi_eq, Real phi_solver_rel,
                                              OphelieDoc214E0Snapshot &out)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Vecd>(particles, names.a_src_imag);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *ar = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Vecd *ai = particles.getVariableDataByName<Vecd>(names.a_src_imag);
    const Real *pr = particles.getVariableDataByName<Real>(names.phi_real);
    const Real *pi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);

    out.tag = tag;
    out.phi_eq = phi_eq;
    out.phi_solver_rel = phi_solver_rel;
    out.vol.resize(n);
    out.ar.resize(n);
    out.ai.resize(n);
    out.pr.resize(n);
    out.pi.resize(n);
    out.jr.resize(n);
    out.ji.resize(n);
    out.q.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        out.vol[i] = vol[i];
        out.ar[i] = ar[i];
        out.ai[i] = ai[i];
        out.pr[i] = pr[i];
        out.pi[i] = pi[i];
        out.jr[i] = jr[i];
        out.ji[i] = ji[i];
        out.q[i] = q[i];
    }
}

inline void hostOphelieDoc214E0DumpSnapshotCsv(const std::string &path, const Vecd *pos, const Vecd &center,
                                              Real radius, Real half_height, const OphelieDoc214E0Snapshot &s)
{
    if (path.empty())
        return;
    std::ofstream csv(path);
    csv << "i,x,y,z,r,region,Vol,"
           "Arx,Ary,Arz,Aix,Aiy,Aiz,phi_r,phi_i,"
           "Jrx,Jry,Jrz,Jix,Jiy,Jiz,Q\n";
    for (size_t i = 0; i < s.vol.size(); ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        csv << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << ","
            << ophelieDoc199E2RegionTag(pos[i], center, radius, half_height) << "," << s.vol[i] << ","
            << s.ar[i][0] << "," << s.ar[i][1] << "," << s.ar[i][2] << "," << s.ai[i][0] << "," << s.ai[i][1]
            << "," << s.ai[i][2] << "," << s.pr[i] << "," << s.pi[i] << "," << s.jr[i][0] << "," << s.jr[i][1]
            << "," << s.jr[i][2] << "," << s.ji[i][0] << "," << s.ji[i][1] << "," << s.ji[i][2] << "," << s.q[i]
            << "\n";
    }
}

inline OphelieDoc214E0Report hostOphelieDoc214E0Audit(const Vecd *pos, const Vecd &center, Real radius,
                                                     Real half_height, Real omega, Real sigma, const Vecd &a0,
                                                     bool a0_times_i, const OphelieDoc214E0Snapshot &snap_r,
                                                     const OphelieDoc214E0Snapshot &snap_c,
                                                     const OphelieDoc214E0Snapshot &snap_rc,
                                                     const std::string &summary_csv_path = {})
{
    OphelieDoc214E0Report rep;
    const size_t n = snap_r.vol.size();
    rep.a0_mag = a0.norm();
    rep.omega = omega;
    rep.jchar = sigma * omega * rep.a0_mag;

    std::vector<char> mask_all, mask_int, mask_side, mask_end;
    ophelieDoc214E0FillRegionMasks(pos, n, center, radius, half_height, mask_all, mask_int, mask_side, mask_end);

    rep.j_rel_rc_vs_r_all =
        ophelieDoc199E3RegionRelDiff(snap_r.vol, snap_r.jr, snap_r.ji, snap_rc.jr, snap_rc.ji, mask_all);
    rep.j_rel_rc_vs_r_interior =
        ophelieDoc199E3RegionRelDiff(snap_r.vol, snap_r.jr, snap_r.ji, snap_rc.jr, snap_rc.ji, mask_int);
    rep.j_rel_rc_vs_r_side =
        ophelieDoc199E3RegionRelDiff(snap_r.vol, snap_r.jr, snap_r.ji, snap_rc.jr, snap_rc.ji, mask_side);
    rep.j_rel_rc_vs_r_end =
        ophelieDoc199E3RegionRelDiff(snap_r.vol, snap_r.jr, snap_r.ji, snap_rc.jr, snap_rc.ji, mask_end);

    rep.q_rel_rc_vs_r_all = ophelieDoc214E0QFieldRel(snap_r.vol, snap_r.q, snap_rc.q, mask_all);
    rep.q_rel_rc_vs_r_interior = ophelieDoc214E0QFieldRel(snap_r.vol, snap_r.q, snap_rc.q, mask_int);
    rep.q_rel_rc_vs_r_side = ophelieDoc214E0QFieldRel(snap_r.vol, snap_r.q, snap_rc.q, mask_side);
    rep.q_rel_rc_vs_r_end = ophelieDoc214E0QFieldRel(snap_r.vol, snap_r.q, snap_rc.q, mask_end);

    double p_r = 0.0, p_c = 0.0, p_rc = 0.0, vol_tot = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double v = static_cast<double>(snap_r.vol[i]);
        vol_tot += v;
        p_r += v * static_cast<double>(snap_r.q[i]);
        p_c += v * static_cast<double>(snap_c.q[i]);
        p_rc += v * static_cast<double>(snap_rc.q[i]);
    }
    rep.p_r = static_cast<Real>(p_r);
    rep.p_c = static_cast<Real>(p_c);
    rep.p_rc = static_cast<Real>(p_rc);
    rep.p_rel_rc_vs_r = std::abs(rep.p_rc - rep.p_r) / (std::abs(rep.p_r) + TinyReal);

    rep.c_all = ophelieDoc214E0JRegion(snap_c.vol, snap_c.jr, snap_c.ji, snap_c.q, mask_all, rep.jchar);
    rep.c_interior = ophelieDoc214E0JRegion(snap_c.vol, snap_c.jr, snap_c.ji, snap_c.q, mask_int, rep.jchar);
    rep.c_side = ophelieDoc214E0JRegion(snap_c.vol, snap_c.jr, snap_c.ji, snap_c.q, mask_side, rep.jchar);
    rep.c_end = ophelieDoc214E0JRegion(snap_c.vol, snap_c.jr, snap_c.ji, snap_c.q, mask_end, rep.jchar);

    // Expected: φ_RC ≈ φ_R − i ω A0·(x−xc) + C  (A0 real → imag shift −ω A0·Δx)
    // If a0_times_i: A0 → i A0, then −i ω (i A0)·Δx = +ω A0·Δx on real part.
    double raw_num = 0.0, raw_den = 0.0;
    double dpr_sum = 0.0, dpi_sum = 0.0;
    std::vector<Real> dpr(n), dpi(n);
    for (size_t i = 0; i < n; ++i)
    {
        const Real xr = pos[i][0] - center[0];
        const Real yr = pos[i][1] - center[1];
        const Real zr = pos[i][2] - center[2];
        const Real adot = a0[0] * xr + a0[1] * yr + a0[2] * zr;
        Real expect_pr = snap_r.pr[i];
        Real expect_pi = snap_r.pi[i];
        if (a0_times_i)
        {
            // A0_eff = i * A0_real ⇒ −i ω A0_eff ·Δx = + ω A0·Δx (real)
            expect_pr = snap_r.pr[i] + omega * adot;
            expect_pi = snap_r.pi[i];
        }
        else
        {
            expect_pr = snap_r.pr[i];
            expect_pi = snap_r.pi[i] - omega * adot;
        }
        dpr[i] = snap_rc.pr[i] - expect_pr;
        dpi[i] = snap_rc.pi[i] - expect_pi;
        const double v = static_cast<double>(snap_r.vol[i]);
        dpr_sum += v * static_cast<double>(dpr[i]);
        dpi_sum += v * static_cast<double>(dpi[i]);
        raw_num += v * static_cast<double>(dpr[i] * dpr[i] + dpi[i] * dpi[i]);
        raw_den += v * static_cast<double>(expect_pr * expect_pr + expect_pi * expect_pi);
    }
    const Real mean_r = (vol_tot > TinyReal) ? static_cast<Real>(dpr_sum / vol_tot) : Real(0);
    const Real mean_i = (vol_tot > TinyReal) ? static_cast<Real>(dpi_sum / vol_tot) : Real(0);
    double align_num = 0.0, align_den = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double v = static_cast<double>(snap_r.vol[i]);
        const Real er = dpr[i] - mean_r;
        const Real ei = dpi[i] - mean_i;
        align_num += v * static_cast<double>(er * er + ei * ei);
        // scale: ω |A0| L_char with L_char ~ R
        const Real scale = omega * rep.a0_mag * radius;
        align_den += v * static_cast<double>(scale * scale);
        (void)align_den;
    }
    // Prefer relative to expected φ magnitude when available; else ω|A0|R.
    rep.phi_raw_rel_no_align = (raw_den > TinyReal) ? static_cast<Real>(std::sqrt(raw_num / raw_den)) : Real(0);
    if (raw_den > TinyReal)
        rep.phi_transform_rel = static_cast<Real>(std::sqrt(align_num / raw_den));
    else
    {
        double den = 0.0;
        const Real scale = omega * rep.a0_mag * radius;
        for (size_t i = 0; i < n; ++i)
            den += static_cast<double>(snap_r.vol[i]) * static_cast<double>(scale * scale);
        rep.phi_transform_rel = (den > TinyReal) ? static_cast<Real>(std::sqrt(align_num / den)) : Real(0);
    }

    rep.gate_j_rc = rep.j_rel_rc_vs_r_all <= Real(0.01);
    rep.gate_q_rc = rep.q_rel_rc_vs_r_all <= Real(0.02);
    rep.gate_c = rep.c_all.j_rel_char <= Real(0.01);
    rep.passed = rep.gate_j_rc && rep.gate_q_rc && rep.gate_c;

    std::cout << "[ophelie][doc214][e0] a0_mag=" << rep.a0_mag << " jchar=" << rep.jchar
              << " a0_times_i=" << (a0_times_i ? 1 : 0) << std::endl;
    std::cout << "[ophelie][doc214][e0] RC_vs_R J_rel all=" << rep.j_rel_rc_vs_r_all
              << " interior=" << rep.j_rel_rc_vs_r_interior << " side=" << rep.j_rel_rc_vs_r_side
              << " end=" << rep.j_rel_rc_vs_r_end << std::endl;
    std::cout << "[ophelie][doc214][e0] RC_vs_R Q_rel all=" << rep.q_rel_rc_vs_r_all
              << " interior=" << rep.q_rel_rc_vs_r_interior << " side=" << rep.q_rel_rc_vs_r_side
              << " end=" << rep.q_rel_rc_vs_r_end << " P_rel=" << rep.p_rel_rc_vs_r << " P_R=" << rep.p_r
              << " P_RC=" << rep.p_rc << " P_C=" << rep.p_c << std::endl;
    std::cout << "[ophelie][doc214][e0] C Jrms/Jchar all=" << rep.c_all.j_rel_char
              << " interior=" << rep.c_interior.j_rel_char << " side=" << rep.c_side.j_rel_char
              << " end=" << rep.c_end.j_rel_char << " Jrms_all=" << rep.c_all.j_rms << std::endl;
    std::cout << "[ophelie][doc214][e0] phi_transform_rel(demean_C)=" << rep.phi_transform_rel
              << " phi_raw_rel_no_align=" << rep.phi_raw_rel_no_align << " mean_C_r=" << mean_r
              << " mean_C_i=" << mean_i << std::endl;
    std::cout << "[ophelie][doc214][e0] residuals R phi_eq=" << snap_r.phi_eq
              << " phi_solver_rel=" << snap_r.phi_solver_rel << " | C phi_eq=" << snap_c.phi_eq
              << " phi_solver_rel=" << snap_c.phi_solver_rel << " | RC phi_eq=" << snap_rc.phi_eq
              << " phi_solver_rel=" << snap_rc.phi_solver_rel << std::endl;
    std::cout << "[ophelie][doc214][e0] gates J_RC<=1%=" << (rep.gate_j_rc ? 1 : 0)
              << " Q_RC<=2%=" << (rep.gate_q_rc ? 1 : 0) << " C_Jchar<=1%=" << (rep.gate_c ? 1 : 0)
              << " passed=" << (rep.passed ? 1 : 0) << std::endl;

    if (!summary_csv_path.empty())
    {
        std::ofstream csv(summary_csv_path);
        csv << "metric,value\n";
        csv << "a0_mag," << rep.a0_mag << "\n";
        csv << "jchar," << rep.jchar << "\n";
        csv << "j_rel_rc_vs_r_all," << rep.j_rel_rc_vs_r_all << "\n";
        csv << "j_rel_rc_vs_r_interior," << rep.j_rel_rc_vs_r_interior << "\n";
        csv << "j_rel_rc_vs_r_side," << rep.j_rel_rc_vs_r_side << "\n";
        csv << "j_rel_rc_vs_r_end," << rep.j_rel_rc_vs_r_end << "\n";
        csv << "q_rel_rc_vs_r_all," << rep.q_rel_rc_vs_r_all << "\n";
        csv << "q_rel_rc_vs_r_interior," << rep.q_rel_rc_vs_r_interior << "\n";
        csv << "q_rel_rc_vs_r_side," << rep.q_rel_rc_vs_r_side << "\n";
        csv << "q_rel_rc_vs_r_end," << rep.q_rel_rc_vs_r_end << "\n";
        csv << "p_rel_rc_vs_r," << rep.p_rel_rc_vs_r << "\n";
        csv << "c_j_rel_char_all," << rep.c_all.j_rel_char << "\n";
        csv << "c_j_rel_char_interior," << rep.c_interior.j_rel_char << "\n";
        csv << "c_j_rel_char_side," << rep.c_side.j_rel_char << "\n";
        csv << "c_j_rel_char_end," << rep.c_end.j_rel_char << "\n";
        csv << "phi_transform_rel," << rep.phi_transform_rel << "\n";
        csv << "phi_raw_rel_no_align," << rep.phi_raw_rel_no_align << "\n";
        csv << "gate_j_rc," << (rep.gate_j_rc ? 1 : 0) << "\n";
        csv << "gate_q_rc," << (rep.gate_q_rc ? 1 : 0) << "\n";
        csv << "gate_c," << (rep.gate_c ? 1 : 0) << "\n";
        csv << "passed," << (rep.passed ? 1 : 0) << "\n";
        std::cout << "[ophelie][doc214][e0] summary_csv=" << summary_csv_path << std::endl;
    }
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC214_E0_CONSTANT_A_INVARIANCE_H
