/**
 * @file electromagnetic_ophelie_doc201_a1_phase_rotate.h
 * @brief Doc201 A1: compare solve(A) vs solve(iA) through the same complex edge-flux entry.
 *
 * Expectation (gauge-adjusted): J(iA)=i J(A), phi(iA)=i phi(A), Q invariant.
 * Does not swap stored solutions — caller must re-prescribe A'=iA and re-solve.
 */
#pragma once

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"

#include <algorithm>
#include <cmath>
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

struct OphelieDoc201A1PhaseRotateReport
{
    Real j_rel_l2 = 1;
    Real q_rel = 1;
    Real phi_rel_l2 = 1;
    Real amp_j0 = 0;
    Real amp_q0 = 0;
    size_t n = 0;
    bool gate_j = false;
    bool gate_q = false;
    bool passed = false;
};

/** In-place A := iA and B := iB on coil fields (host). */
inline void hostOphelieMultiplyCoilFieldsByI(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                             size_t n)
{
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.a_coil_imag);
    syncVariableToHost<Vecd>(particles, names.b_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_coil_imag);
    Vecd *ar = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    Vecd *ai = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
    Vecd *br = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    Vecd *bi = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd ar0 = ar[i];
        const Vecd ai0 = ai[i];
        ar[i] = -ai0;
        ai[i] = ar0;
        const Vecd br0 = br[i];
        const Vecd bi0 = bi[i];
        br[i] = -bi0;
        bi[i] = br0;
    }
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
    syncVariableToDevice<Vecd>(particles, names.b_coil_real);
    syncVariableToDevice<Vecd>(particles, names.b_coil_imag);
}

/**
 * Compare baseline snapshots (jr0,ji0,pr0,pi0,q0) to current fields after solve(iA).
 * Expected: Jr'≈-Ji0, Ji'≈Jr0; phi similarly after removing complex mean; Q'≈Q0.
 */
inline OphelieDoc201A1PhaseRotateReport hostOphelieDoc201A1ComparePhaseRotate(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const std::vector<Vecd> &jr0,
    const std::vector<Vecd> &ji0, const std::vector<Real> &pr0, const std::vector<Real> &pi0,
    const std::vector<Real> &q0, Real j_tol = Real(0.01), Real q_tol = Real(0.02),
    const std::string &dump_csv = "")
{
    OphelieDoc201A1PhaseRotateReport rep;
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *pr = particles.getVariableDataByName<Real>(names.phi_real);
    const Real *pi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = jr0.size();
    rep.n = n;

    // Complex mean of phi0 and phi' for gauge removal.
    double vsum = 0, p0r = 0, p0i = 0, p1r = 0, p1i = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double w = static_cast<double>(vol[i]);
        vsum += w;
        p0r += w * static_cast<double>(pr0[i]);
        p0i += w * static_cast<double>(pi0[i]);
        p1r += w * static_cast<double>(pr[i]);
        p1i += w * static_cast<double>(pi[i]);
    }
    const Real m0r = vsum > TinyReal ? static_cast<Real>(p0r / vsum) : 0;
    const Real m0i = vsum > TinyReal ? static_cast<Real>(p0i / vsum) : 0;
    const Real m1r = vsum > TinyReal ? static_cast<Real>(p1r / vsum) : 0;
    const Real m1i = vsum > TinyReal ? static_cast<Real>(p1i / vsum) : 0;

    double err_j2 = 0, ref_j2 = 0, err_p2 = 0, ref_p2 = 0, err_q2 = 0, ref_q2 = 0;
    std::ofstream csv;
    if (!dump_csv.empty())
    {
        csv.open(dump_csv);
        csv << "i,vol,Jr0x,Jr0y,Jr0z,Ji0x,Ji0y,Ji0z,Jr1x,Jr1y,Jr1z,Ji1x,Ji1y,Ji1z,Q0,Q1,dJr_vs_mJi,dJi_vs_Jr\n";
    }
    for (size_t i = 0; i < n; ++i)
    {
        const Real w = vol[i];
        // i*J0 = (-Ji0, Jr0)
        const Vecd exp_jr = -ji0[i];
        const Vecd exp_ji = jr0[i];
        const Vecd djr = jr[i] - exp_jr;
        const Vecd dji = ji[i] - exp_ji;
        err_j2 += static_cast<double>(w) * static_cast<double>(djr.squaredNorm() + dji.squaredNorm());
        ref_j2 += static_cast<double>(w) * static_cast<double>(jr0[i].squaredNorm() + ji0[i].squaredNorm());

        const Real p0rc = pr0[i] - m0r;
        const Real p0ic = pi0[i] - m0i;
        const Real exp_pr = -p0ic;
        const Real exp_pi = p0rc;
        const Real p1rc = pr[i] - m1r;
        const Real p1ic = pi[i] - m1i;
        const Real dpr = p1rc - exp_pr;
        const Real dpi = p1ic - exp_pi;
        err_p2 += static_cast<double>(w) * static_cast<double>(dpr * dpr + dpi * dpi);
        ref_p2 += static_cast<double>(w) * static_cast<double>(p0rc * p0rc + p0ic * p0ic);

        const Real dq = q[i] - q0[i];
        err_q2 += static_cast<double>(w) * static_cast<double>(dq * dq);
        ref_q2 += static_cast<double>(w) * static_cast<double>(q0[i] * q0[i]);

        if (csv.is_open())
        {
            csv << i << "," << w << "," << jr0[i][0] << "," << jr0[i][1] << "," << jr0[i][2] << "," << ji0[i][0]
                << "," << ji0[i][1] << "," << ji0[i][2] << "," << jr[i][0] << "," << jr[i][1] << "," << jr[i][2]
                << "," << ji[i][0] << "," << ji[i][1] << "," << ji[i][2] << "," << q0[i] << "," << q[i] << ","
                << djr.norm() << "," << dji.norm() << "\n";
        }
    }
    rep.j_rel_l2 = static_cast<Real>(std::sqrt(err_j2) / std::max(std::sqrt(ref_j2), static_cast<double>(TinyReal)));
    rep.phi_rel_l2 =
        static_cast<Real>(std::sqrt(err_p2) / std::max(std::sqrt(ref_p2), static_cast<double>(TinyReal)));
    rep.q_rel = static_cast<Real>(std::sqrt(err_q2) / std::max(std::sqrt(ref_q2), static_cast<double>(TinyReal)));
    rep.amp_j0 = static_cast<Real>(std::sqrt(ref_j2 / std::max(vsum, 1e-30)));
    rep.amp_q0 = static_cast<Real>(std::sqrt(ref_q2 / std::max(vsum, 1e-30)));
    rep.gate_j = rep.j_rel_l2 <= j_tol;
    rep.gate_q = rep.q_rel <= q_tol;
    rep.passed = rep.gate_j && rep.gate_q;

    std::cout << "[ophelie][doc201][a1] phase_rotate J_rel_L2=" << rep.j_rel_l2 << " Q_rel=" << rep.q_rel
              << " phi_rel_L2=" << rep.phi_rel_l2 << " |J0|_rms=" << rep.amp_j0 << " gate_j=" << (rep.gate_j ? 1 : 0)
              << " gate_q=" << (rep.gate_q ? 1 : 0) << " passed=" << (rep.passed ? 1 : 0)
              << " tol_J=" << j_tol << " tol_Q=" << q_tol << std::endl;
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH
