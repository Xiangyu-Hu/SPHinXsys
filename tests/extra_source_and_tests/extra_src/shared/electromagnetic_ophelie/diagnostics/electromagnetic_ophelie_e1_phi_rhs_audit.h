/**
 * @file electromagnetic_ophelie_e1_phi_rhs_audit.h
 * @brief Doc198: prescribed-A φ-RHS audit (why discrete RHS ≠ 0 when div A = 0).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_E1_PHI_RHS_AUDIT_H
#define ELECTROMAGNETIC_OPHELIE_E1_PHI_RHS_AUDIT_H

#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_boundary_diagnostics.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "update_body_relation.h"

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

struct OphelieE1ScalarRegionStats
{
    size_t n = 0;
    double w = 0.0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double max_abs = 0.0;
};

inline void ophelieE1ScalarAccum(OphelieE1ScalarRegionStats &s, Real vol, Real v)
{
    const double w = static_cast<double>(vol);
    s.w += w;
    s.sum += w * static_cast<double>(v);
    s.sum_sq += w * static_cast<double>(v) * static_cast<double>(v);
    s.max_abs = std::max(s.max_abs, static_cast<double>(std::abs(v)));
    ++s.n;
}

inline Real ophelieE1ScalarRms(const OphelieE1ScalarRegionStats &s)
{
    return (s.w > TinyReal) ? static_cast<Real>(std::sqrt(s.sum_sq / s.w)) : Real(0);
}
inline Real ophelieE1ScalarMean(const OphelieE1ScalarRegionStats &s)
{
    return (s.w > TinyReal) ? static_cast<Real>(s.sum / s.w) : Real(0);
}

struct OphelieE1PhiRhsOperatorReport
{
    const char *operator_name = "";
    OphelieE1ScalarRegionStats all;
    OphelieE1ScalarRegionStats interior; // r<0.7R, |z|<0.5 Hh
    OphelieE1ScalarRegionStats rim;      // r>=0.75R
    OphelieE1ScalarRegionStats endband;  // |z|>0.8 Hh
    Real vol_weighted_l2 = 0.0;          // fingerprint-style √(∑v b²/∑v) — same as rms if mean0
    Real omega_a_vol_l2 = 0.0;           // √〈|ωA|²〉 scale
    Real rhs_over_omega_a = 0.0;         // all.rms / omega_a_vol_l2 (dimensionless edge-scale mix)
};

struct OphelieE1PhiRhsAuditReport
{
    OphelieE1PhiRhsOperatorReport legacy_edge; // production EdgeFlux LegacyFlux path
    OphelieE1PhiRhsOperatorReport div_sigma_a;
    Real cosine_legacy_vs_div = 0.0; // vol-weighted cosine of two RHS fields
    Real rel_l2_legacy_vs_div = 0.0;
};

inline void hostOphelieE1FillScalarRegions(BaseParticles &particles, const Real *field, const Vecd &axis_origin,
                                          Real radius_m, Real half_height, OphelieE1PhiRhsOperatorReport &rep)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real zc = axis_origin[2];
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real abs_z = std::abs(pos[i][2] - zc);
        const Real v = field[i];
        ophelieE1ScalarAccum(rep.all, vol[i], v);
        if (r < Real(0.70) * radius_m && abs_z < Real(0.50) * half_height)
            ophelieE1ScalarAccum(rep.interior, vol[i], v);
        if (r >= Real(0.75) * radius_m)
            ophelieE1ScalarAccum(rep.rim, vol[i], v);
        if (abs_z > Real(0.80) * half_height)
            ophelieE1ScalarAccum(rep.endband, vol[i], v);
    }
    rep.vol_weighted_l2 = ophelieE1ScalarRms(rep.all);
}

inline void logOphelieE1PhiRhsOperatorReport(const OphelieE1PhiRhsOperatorReport &r)
{
    auto one = [](const char *nm, const OphelieE1ScalarRegionStats &s) {
        std::cout << " " << nm << "{n=" << s.n << " rms=" << ophelieE1ScalarRms(s)
                  << " mean=" << ophelieE1ScalarMean(s) << " max|b|=" << s.max_abs << "}";
    };
    std::cout << "[ophelie][e1][phi-rhs] op=" << r.operator_name << " omegaA_L2=" << r.omega_a_vol_l2
              << " rhs/omegaA=" << r.rhs_over_omega_a;
    one("all", r.all);
    one("interior", r.interior);
    one("rim", r.rim);
    one("end", r.endband);
    std::cout << std::endl;
}

/**
 * Build φ imag RHS twice (LegacyFlux/edge + DivSigmaA), report region pollution.
 * Restores params.phi_rhs_operator_kind_ and leaves phi_rhs_imag = production (legacy) field.
 */
template <class ExecutionPolicy>
inline OphelieE1PhiRhsAuditReport execOphelieE1PrescribedAPhiRhsAudit(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const Vecd &axis_origin, Real radius_m, Real half_height, Real b0, const std::string &csv_path = {})
{
    OphelieE1PhiRhsAuditReport out;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const OpheliePhiRhsOperatorKind saved = params.phi_rhs_operator_kind_;
    params.use_a_total_for_edge_flux_ = true;

    // ω|A| scale
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *a_coil = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    double wa2 = 0.0, wtot = 0.0;
    const Real omega = params.omega();
    for (size_t i = 0; i < n; ++i)
    {
        const double w = static_cast<double>(vol[i]);
        const Real an = (omega * a_coil[i]).norm();
        wa2 += w * static_cast<double>(an * an);
        wtot += w;
    }
    const Real omega_a_l2 = (wtot > TinyReal) ? static_cast<Real>(std::sqrt(wa2 / wtot)) : Real(0);

    StdVec<Real> rhs_legacy(n, Real(0));
    StdVec<Real> rhs_div(n, Real(0));

    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, nullptr, 0.0, nullptr);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_legacy.data(), n);
    out.legacy_edge.operator_name = "legacy-flux/edge";
    out.legacy_edge.omega_a_vol_l2 = omega_a_l2;
    hostOphelieE1FillScalarRegions(particles, rhs_legacy.data(), axis_origin, radius_m, half_height, out.legacy_edge);
    out.legacy_edge.rhs_over_omega_a =
        out.legacy_edge.vol_weighted_l2 / (omega_a_l2 + TinyReal);

    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, nullptr, 0.0, nullptr);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_div.data(), n);
    out.div_sigma_a.operator_name = "div-sigma-a";
    out.div_sigma_a.omega_a_vol_l2 = omega_a_l2;
    hostOphelieE1FillScalarRegions(particles, rhs_div.data(), axis_origin, radius_m, half_height, out.div_sigma_a);
    out.div_sigma_a.rhs_over_omega_a =
        out.div_sigma_a.vol_weighted_l2 / (omega_a_l2 + TinyReal);

    // Cosine / relative difference
    double dot = 0.0, n1 = 0.0, n2 = 0.0, d2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double w = static_cast<double>(vol[i]);
        const double a = static_cast<double>(rhs_legacy[i]);
        const double b = static_cast<double>(rhs_div[i]);
        dot += w * a * b;
        n1 += w * a * a;
        n2 += w * b * b;
        d2 += w * (a - b) * (a - b);
    }
    const double den = std::sqrt(n1 * n2);
    out.cosine_legacy_vs_div = (den > TinyReal) ? static_cast<Real>(dot / den) : Real(0);
    out.rel_l2_legacy_vs_div =
        (n1 > TinyReal) ? static_cast<Real>(std::sqrt(d2 / n1)) : Real(0);

    logOphelieE1PhiRhsOperatorReport(out.legacy_edge);
    logOphelieE1PhiRhsOperatorReport(out.div_sigma_a);
    std::cout << "[ophelie][e1][phi-rhs] cosine(legacy,div)=" << out.cosine_legacy_vs_div
              << " rel_L2(legacy-div)/legacy=" << out.rel_l2_legacy_vs_div
              << " B0=" << b0 << " note: analytic divA=0 ⇒ expect small interior RHS"
              << std::endl;

    if (!csv_path.empty())
    {
        syncVariableToHost<Vecd>(particles, "Position");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        std::ofstream ofs(csv_path);
        if (ofs)
        {
            ofs << "x,y,z,r,rhs_legacy,rhs_div,abs_legacy,abs_div\n";
            for (size_t i = 0; i < n; ++i)
            {
                const Real dx = pos[i][0] - axis_origin[0];
                const Real dy = pos[i][1] - axis_origin[1];
                const Real r = std::sqrt(dx * dx + dy * dy);
                ofs << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << "," << rhs_legacy[i] << ","
                    << rhs_div[i] << "," << std::abs(rhs_legacy[i]) << "," << std::abs(rhs_div[i]) << "\n";
            }
            std::cout << "[ophelie][e1][phi-rhs] csv=" << csv_path << std::endl;
        }
    }

    // Restore production RHS on particles for subsequent L3 solve
    params.phi_rhs_operator_kind_ = saved;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, nullptr, 0.0, nullptr);
    return out;
}

/** After φ solve: compare |∇φ| to |ωA| by region (expects ≪1 for prescribed A). */
inline void hostOphelieE1LogGradPhiVsOmegaA(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                           const Vecd &axis_origin, Real radius_m, Real half_height, Real omega)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Real>(particles, names.phi_imag);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Vecd *g = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);

    OphelieE1ScalarRegionStats ratio_all, ratio_int, ratio_rim;
    OphelieE1ScalarRegionStats phi_all, phi_int, phi_rim;
    const Real zc = axis_origin[2];
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real abs_z = std::abs(pos[i][2] - zc);
        const Real oa = (omega * a[i]).norm() + TinyReal;
        const Real ratio = g[i].norm() / oa;
        ophelieE1ScalarAccum(ratio_all, vol[i], ratio);
        ophelieE1ScalarAccum(phi_all, vol[i], phi[i]);
        if (r < Real(0.70) * radius_m && abs_z < Real(0.50) * half_height)
        {
            ophelieE1ScalarAccum(ratio_int, vol[i], ratio);
            ophelieE1ScalarAccum(phi_int, vol[i], phi[i]);
        }
        if (r >= Real(0.75) * radius_m)
        {
            ophelieE1ScalarAccum(ratio_rim, vol[i], ratio);
            ophelieE1ScalarAccum(phi_rim, vol[i], phi[i]);
        }
    }
    std::cout << "[ophelie][e1][phi-pollution] |gradPhi|/|omegaA| mean: all=" << ophelieE1ScalarMean(ratio_all)
              << " interior=" << ophelieE1ScalarMean(ratio_int) << " rim=" << ophelieE1ScalarMean(ratio_rim)
              << " | phi_rms all=" << ophelieE1ScalarRms(phi_all) << " int=" << ophelieE1ScalarRms(phi_int)
              << " rim=" << ophelieE1ScalarRms(phi_rim) << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_E1_PHI_RHS_AUDIT_H
