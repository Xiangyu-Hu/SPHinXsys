/**
 * @file electromagnetic_ophelie_doc218_d2_manufactured.h
 * @brief Doc218 C: M1–M8 manufactured suite on the production D2 call chain.
 *
 * Injects F via A_src (F = a_sign ω A; imag uses a_sign=+1 on ASrcReal), then
 * execOphelieComplexEdgeFluxSolveReconAndPower (D2 solve + local-poly recon).
 * No separate "correct" solver.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC218_D2_MANUFACTURED_H
#define ELECTROMAGNETIC_OPHELIE_DOC218_D2_MANUFACTURED_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_phi_d2_local_poly.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi_component.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc218MfgCaseResult
{
    std::string id;
    Real eJ = Real(0);
    Real q_rel = Real(-1);
    Real jchar = Real(0);
    Real phi_rel = Real(0);
    Real rank_fail = Real(0);
    Real neu_max_over_Fchar = Real(-1);
    bool gate = false;
    bool zero_ref = false;
};

struct OphelieDoc218MfgSuiteResult
{
    std::vector<OphelieDoc218MfgCaseResult> cases;
    bool passed = false;
};

inline Vecd doc218RotA(const Vecd &pos, const Vecd &c, Real B0)
{
    // A = 0.5 B0 × r_cyl  (uniform Bz) → A=(-0.5 B0 y', 0.5 B0 x', 0)
    const Real x = pos[0] - c[0];
    const Real y = pos[1] - c[1];
    return Vecd(Real(-0.5) * B0 * y, Real(0.5) * B0 * x, Real(0));
}

inline Vecd doc218GradPsiHarmonic(const Vecd &pos, const Vecd &c, Real R, Real psi0)
{
    const Real x = pos[0] - c[0];
    const Real y = pos[1] - c[1];
    return Vecd(Real(2) * psi0 * x / (R * R), Real(-2) * psi0 * y / (R * R), Real(0));
}

inline Vecd doc218GradPsiM7(const Vecd &pos, const Vecd &c, Real R, Real psi0)
{
    const Real x = pos[0] - c[0];
    return Vecd(Real(2) * psi0 * x / (R * R), Real(0), Real(0));
}

inline Vecd doc218GradPsiM8(const Vecd &pos, const Vecd &c, Real H, Real psi0)
{
    const Real z = pos[2] - c[2];
    return Vecd(Real(0), Real(0), Real(2) * psi0 * z / (H * H));
}

inline Vecd doc218Jstar(const Vecd &pos, const Vecd &c, Real R, Real J0)
{
    const Real x = pos[0] - c[0];
    const Real y = pos[1] - c[1];
    return J0 * Vecd(-y / R, x / R, Real(0));
}

inline void doc218WriteWallCsv(const std::string &path, const Vecd *pos, const Real *Vol, const Real *sigma,
                               const Vecd *F, const Real *phi, const Vecd *J, const Real *Q, size_t n,
                               const Vecd &c, Real R, Real H, Real dp)
{
    std::ofstream os(path);
    os << "x,y,z,Vol,sigma,Fx,Fy,Fz,phi,Jx,Jy,Jz,Q,is_side,is_top,is_bot,wx,wy,wz,nx,ny,nz\n";
    const Real wall_dp = Real(0.75) * dp;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd x = pos[i] - c;
        const Real r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
        const Real z = x[2];
        const bool is_side = std::abs(r - R) <= wall_dp;
        const bool is_top = (!is_side) && (std::abs(z - H) <= wall_dp);
        const bool is_bot = (!is_side) && (std::abs(z + H) <= wall_dp);
        Vecd wall = pos[i];
        Vecd nvec = Vecd::Zero();
        if (is_side)
        {
            const Real rh = std::max(r, TinyReal);
            wall[0] = c[0] + R * x[0] / rh;
            wall[1] = c[1] + R * x[1] / rh;
            nvec[0] = x[0] / rh;
            nvec[1] = x[1] / rh;
        }
        else if (is_top)
        {
            wall[2] = c[2] + H;
            nvec[2] = Real(1);
        }
        else if (is_bot)
        {
            wall[2] = c[2] - H;
            nvec[2] = Real(-1);
        }
        os << pos[i][0] << ',' << pos[i][1] << ',' << pos[i][2] << ',' << Vol[i] << ',' << sigma[i] << ','
           << F[i][0] << ',' << F[i][1] << ',' << F[i][2] << ',' << phi[i] << ',' << J[i][0] << ',' << J[i][1]
           << ',' << J[i][2] << ',' << Q[i] << ',' << (is_side ? 1 : 0) << ',' << (is_top ? 1 : 0) << ','
           << (is_bot ? 1 : 0) << ',' << wall[0] << ',' << wall[1] << ',' << wall[2] << ',' << nvec[0] << ','
           << nvec[1] << ',' << nvec[2] << '\n';
    }
}

/**
 * Run one manufactured case. F is prescribed at particle centers; ASrcReal = F/ω.
 * Uses production D2 path (solve+recon). Imag-only (skip_real).
 */
template <class ExecutionPolicy>
inline OphelieDoc218MfgCaseResult runOphelieDoc218D2ManufacturedOne(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const Vecd &center, Real R, Real H, Real B0, const std::string &id, const std::string &csv_prefix)
{
    OphelieDoc218MfgCaseResult out;
    out.id = id;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real omega = params.omega();
    const Real A0 = B0 * R / Real(2);
    const Real Fchar = omega * A0;
    const Real psi0 = Fchar * R / Real(2);
    const Real sigma0 = params.sigma_glass_;

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, names.sigma);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Vecd *a_src_r = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Vecd *a_src_i = particles.getVariableDataByName<Vecd>(names.a_src_imag);
    Vecd *a_coil_r = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    Vecd *a_coil_i = particles.getVariableDataByName<Vecd>(names.a_coil_imag);

    std::vector<Vecd> F(n, Vecd::Zero());
    std::vector<Vecd> Jref(n, Vecd::Zero());
    bool zero_ref = false;
    bool use_var_sigma = false;

    if (id == "M1")
    {
        // constant F → linear φ cancel, J*=0
        zero_ref = true;
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            F[i] = Vecd(omega * A0, Real(0), Real(0));
        }
    }
    else if (id == "M2")
    {
        // rotational F = ω A_rot; φ*=0; J* = -σ F
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            const Vecd Arot = doc218RotA(pos[i], center, B0);
            F[i] = omega * Arot;
            Jref[i] = -sigma[i] * F[i];
        }
    }
    else if (id == "M3")
    {
        // R+C: F = ω(A_rot + A0 êx); J* same as M2 (C cancelled by φ)
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            const Vecd Arot = doc218RotA(pos[i], center, B0);
            F[i] = omega * (Arot + Vecd(A0, Real(0), Real(0)));
            Jref[i] = -sigma[i] * (omega * Arot);
        }
    }
    else if (id == "M4")
    {
        zero_ref = true;
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            F[i] = -doc218GradPsiHarmonic(pos[i], center, R, psi0);
        }
    }
    else if (id == "M5")
    {
        const Real J0 = sigma0 * Fchar;
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            const Vecd Js = doc218Jstar(pos[i], center, R, J0);
            Jref[i] = Js;
            F[i] = -doc218GradPsiHarmonic(pos[i], center, R, psi0) - Js / sigma[i];
        }
    }
    else if (id == "M6")
    {
        use_var_sigma = true;
        const Real J0 = sigma0 * Fchar;
        for (size_t i = 0; i < n; ++i)
        {
            const Real x = pos[i][0] - center[0];
            sigma[i] = sigma0 * (Real(1.5) + Real(0.5) * x / R);
            const Vecd Js = doc218Jstar(pos[i], center, R, J0);
            Jref[i] = Js;
            F[i] = -doc218GradPsiHarmonic(pos[i], center, R, psi0) - Js / sigma[i];
        }
    }
    else if (id == "M7")
    {
        zero_ref = true;
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            F[i] = -doc218GradPsiM7(pos[i], center, R, psi0);
        }
    }
    else if (id == "M8")
    {
        zero_ref = true;
        for (size_t i = 0; i < n; ++i)
        {
            sigma[i] = sigma0;
            F[i] = -doc218GradPsiM8(pos[i], center, H, psi0);
        }
    }
    else
    {
        std::cerr << "[ophelie][doc218][mfg] unknown id=" << id << std::endl;
        out.gate = false;
        return out;
    }
    out.zero_ref = zero_ref;
    (void)use_var_sigma;

    // F = (+1) ω A_src_real  ⇒  A = F/ω
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd A = F[i] / (omega + TinyReal);
        a_src_r[i] = A;
        a_coil_r[i] = A;
        a_src_i[i] = Vecd::Zero();
        a_coil_i[i] = Vecd::Zero();
    }
    syncVariableToDevice<Real>(particles, names.sigma);
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
    syncVariableToDevice<Vecd>(particles, names.a_src_imag);
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_coil_imag);

    OphelieComplexEdgeFluxSolveReport report;
    (void)execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(
        glass_body, inner, names, params, nullptr, params.doc217_d2_dp_, &report, /*skip_real=*/true);
    out.phi_rel = report.phi_imag_solver_rel_residual;

    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_imag);
    const Vecd *J = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *Q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);

    Real vol_sum = Real(0), j2 = Real(0), e2 = Real(0), qn2 = Real(0), qr2 = Real(0), f2 = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        vol_sum += Vol[i];
        f2 += Vol[i] * F[i].squaredNorm();
        if (zero_ref)
        {
            e2 += Vol[i] * J[i].squaredNorm();
        }
        else
        {
            const Vecd dJ = J[i] - Jref[i];
            e2 += Vol[i] * dJ.squaredNorm();
            j2 += Vol[i] * Jref[i].squaredNorm();
            const Real Qn = Q[i];
            const Real Qr = Real(0.5) * Jref[i].squaredNorm() / (sigma[i] + TinyReal);
            qn2 += Vol[i] * (Qn - Qr) * (Qn - Qr);
            qr2 += Vol[i] * Qr * Qr;
        }
    }
    const Real sigma_mean = sigma0; // jchar uses base σ even for M6 (matches Python)
    out.jchar = sigma_mean * std::sqrt(f2 / (vol_sum + TinyReal));
    if (zero_ref)
        out.eJ = std::sqrt(e2 / (out.jchar * out.jchar * vol_sum + TinyReal));
    else
    {
        out.eJ = std::sqrt(e2 / (j2 + TinyReal));
        out.q_rel = std::sqrt(qn2 / (qr2 + TinyReal));
    }

    // Gates (Doc218 §6.1–6.2)
    if (id == "M1")
        out.gate = out.eJ <= Real(0.01);
    else if (id == "M2" || id == "M3" || id == "M4" || id == "M5" || id == "M6" || id == "M7" || id == "M8")
        out.gate = out.eJ <= Real(0.05) && (out.q_rel < Real(0) || out.q_rel <= Real(0.05));
    else
        out.gate = false;

    if (report.phi_space != std::string("d2"))
    {
        std::cerr << "[ophelie][doc218][mfg] refuse: phi_space=" << report.phi_space << " (need d2)" << std::endl;
        out.gate = false;
    }

    std::cout << "[ophelie][doc218][mfg] id=" << id << " eJ=" << out.eJ << " q_rel=" << out.q_rel
              << " jchar=" << out.jchar << " phi_rel=" << out.phi_rel << " wall_F_mode="
              << (params.doc218_wall_F_mode_ == 1 ? "wall_reconstructed" : "center_legacy")
              << " compat_mode=" << (params.doc218_compat_mode_ == 1 ? "none" : "boundary_mean_legacy")
              << " gate=" << (out.gate ? 1 : 0) << std::endl;

    if (!csv_prefix.empty())
        doc218WriteWallCsv(csv_prefix + "_" + id + ".csv", pos, Vol, sigma, F.data(), phi, J, Q, n, center, R, H,
                           params.doc217_d2_dp_);

    return out;
}

template <class ExecutionPolicy>
inline OphelieDoc218MfgSuiteResult runOphelieDoc218D2ManufacturedSuite(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const Vecd &center, Real R, Real H, Real B0, const std::vector<std::string> &ids,
    const std::string &csv_prefix)
{
    OphelieDoc218MfgSuiteResult suite;
    if (!params.doc217_d2_local_poly_)
    {
        std::cerr << "[ophelie][doc218][mfg] refuse: require --doc217-d2-local-poly" << std::endl;
        suite.passed = false;
        return suite;
    }
    // Formal modes for acceptance (legacy only via explicit CLI for A diagnosis).
    std::cout << "[ophelie][doc218][mfg] suite n_cases=" << ids.size()
              << " wall_F_mode=" << params.doc218_wall_F_mode_
              << " compat_mode=" << params.doc218_compat_mode_ << std::endl;

    OphelieDoc218MfgCaseResult m2, m3;
    std::vector<Vecd> j_m2_store;
    std::vector<Real> q_m2_store;
    Real dJ_m3_m2 = Real(-1), dQ_m3_m2 = Real(-1);
    bool m3_invariance_ok = true;
    for (const auto &id : ids)
    {
        const auto r = runOphelieDoc218D2ManufacturedOne<ExecutionPolicy>(glass_body, inner, names, params, center, R,
                                                                            H, B0, id, csv_prefix);
        suite.cases.push_back(r);
        BaseParticles &particles = glass_body.getBaseParticles();
        const size_t n = particles.TotalRealParticles();
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        // Manufactured suite uses skip_real=true → complex Q is not written; use imag Q.
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_imag);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Vecd *J = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
        const Real *Qc = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);
        const Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        if (id == "M2")
        {
            m2 = r;
            j_m2_store.assign(J, J + n);
            q_m2_store.assign(Qc, Qc + n);
        }
        if (id == "M3")
        {
            m3 = r;
            // Doc219 B §5.4: full-vector ΔJ(M3,M2) ≤ 1%, Q field ≤ 2%.
            if (j_m2_store.size() != n || q_m2_store.size() != n || Qc == nullptr)
            {
                std::cerr << "[ophelie][doc219][mfg] M3-M2 invariance requires prior M2 in same suite"
                          << std::endl;
                m3_invariance_ok = false;
                suite.cases.back().gate = false;
            }
            else
            {
                Real e2 = Real(0), j2 = Real(0), qe2 = Real(0), q2 = Real(0);
                for (size_t i = 0; i < n; ++i)
                {
                    const Vecd dJ = J[i] - j_m2_store[i];
                    e2 += Vol[i] * dJ.squaredNorm();
                    j2 += Vol[i] * j_m2_store[i].squaredNorm();
                    const Real dQ = Qc[i] - q_m2_store[i];
                    qe2 += Vol[i] * dQ * dQ;
                    q2 += Vol[i] * q_m2_store[i] * q_m2_store[i];
                }
                dJ_m3_m2 = std::sqrt(e2 / (j2 + TinyReal));
                dQ_m3_m2 = std::sqrt(qe2 / (q2 + TinyReal));
                m3_invariance_ok = (dJ_m3_m2 <= Real(0.01)) && (dQ_m3_m2 <= Real(0.02));
                std::cout << "[ophelie][doc219][mfg] M3-M2 invariance dJ=" << dJ_m3_m2 << " dQ=" << dQ_m3_m2
                          << " gate_dJ_1pct_dQ_2pct=" << (m3_invariance_ok ? 1 : 0) << std::endl;
                if (!m3_invariance_ok)
                    suite.cases.back().gate = false;
            }
        }
    }
    suite.passed = true;
    for (const auto &c : suite.cases)
        if (!c.gate)
            suite.passed = false;
    if (!m3.id.empty() && !m3_invariance_ok)
        suite.passed = false;
    return suite;
}

inline std::vector<std::string> parseDoc218MfgIds(const std::string &spec)
{
    if (spec == "all" || spec == "M1-M8" || spec == "m1-m8")
        return {"M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"};
    if (spec == "M1-M4")
        return {"M1", "M2", "M3", "M4"};
    if (spec == "M7-M8")
        return {"M7", "M8"};
    std::vector<std::string> ids;
    std::stringstream ss(spec);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (!tok.empty())
            ids.push_back(tok);
    }
    if (ids.empty() && !spec.empty())
        ids.push_back(spec);
    return ids;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC218_D2_MANUFACTURED_H
