/**
 * @file electromagnetic_ophelie_doc199_p1_power_decomposition.h
 * @brief Doc199 P1: split power error into field / discrete-integral / continuous geometry.
 *
 * Prescribed rotational A, no self-induction. Peak phasor: Q=|J|²/(2σ).
 * Continuous: P = π H σ ω² B0² R⁴ / 16.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC199_P1_POWER_DECOMPOSITION_H
#define ELECTROMAGNETIC_OPHELIE_DOC199_P1_POWER_DECOMPOSITION_H

#include "electromagnetic_ophelie_analytic_cylinder_self_induction.h"
#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_field_names.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc199P1PowerReport
{
    Real p_num_particle = 0.0;      // Σ V Q_num (edge recon complex)
    Real p_ref_discrete = 0.0;      // Σ V |J_ref|²/(2σ) on same particles
    Real p_ref_continuous = 0.0;    // π H σ ω² B0² R⁴ / 16
    Real p_edge_graph = 0.0;
    Real p_edge_recon = 0.0;
    Real p_j_sq_over_2sigma = 0.0;  // from numerical J recon identity
    Real vol_particles = 0.0;
    Real vol_continuous = 0.0;
    Real sum_v_r2 = 0.0;
    Real moment_r2_continuous = 0.0; // ∫ r² dV = π H R⁴ / 2
    Real field_rel = 0.0;            // |P_num - P_ref_disc| / P_ref_disc
    Real geometry_rel = 0.0;         // |P_ref_disc - P_ref_cont| / P_ref_cont
    Real total_vs_cont_rel = 0.0;    // |P_num - P_ref_cont| / P_ref_cont
    Real edge_vs_particle_rel = 0.0;
    Real q_outer_over_core_num = 0.0;
    Real q_outer_over_core_ref_disc = 0.0;
    Real q_outer_over_core_cont = 0.0;
    Real j_rel_l2_vs_ref = 0.0;
};

inline OphelieDoc199P1PowerReport hostOphelieDoc199P1PowerDecomposition(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const Vecd &center, Real radius, Real half_height, Real b0, const std::string &particle_csv_path = {})
{
    OphelieDoc199P1PowerReport rep;
    const size_t n = particles.TotalRealParticles();
    const Real omega = params.omega();
    const Real sigma = params.sigma_glass_;
    const Real height = Real(2) * half_height;

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *q_num = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);

    const OphelieEdgeFluxPowerMetrics edge_m = computeHostEdgeFluxPowerMetrics(particles, names, n);
    rep.p_edge_graph = edge_m.p_graph_edge;
    rep.p_edge_recon = hostEdgeFluxReconPower(particles, names, n, params);

    const Real r_core = Real(0.5) * radius;
    const Real r_outer = Real(0.75) * radius;
    double qn_o = 0, vn_o = 0, qn_c = 0, vn_c = 0;
    double qr_o = 0, vr_o = 0, qr_c = 0, vr_c = 0;
    double j_err = 0, j_den = 0;
    double p_jsq = 0;

    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r2 = dx * dx + dy * dy;
        const Real r = std::sqrt(r2);
        const double v = static_cast<double>(vol[i]);
        rep.vol_particles += static_cast<Real>(v);
        rep.sum_v_r2 += static_cast<Real>(v * static_cast<double>(r2));

        const Vecd a_ref = ophelieRotationalAUniformBz(pos[i], center, b0);
        const Vecd j_ref_i = -sigma * omega * a_ref; // imag; real=0
        const Real j_ref2 = j_ref_i.squaredNorm();
        const Real q_ref = j_ref2 / (Real(2) * sigma + TinyReal);

        rep.p_num_particle += static_cast<Real>(v * static_cast<double>(q_num[i]));
        rep.p_ref_discrete += static_cast<Real>(v * static_cast<double>(q_ref));

        const Real jn2 = jr[i].squaredNorm() + ji[i].squaredNorm();
        p_jsq += v * static_cast<double>(jn2 / (Real(2) * sigma + TinyReal));

        const Vecd dj_r = jr[i];
        const Vecd dj_i = ji[i] - j_ref_i;
        j_err += v * static_cast<double>(dj_r.squaredNorm() + dj_i.squaredNorm());
        j_den += v * static_cast<double>(j_ref2);

        // mid-band radial shells for Q ratio (same cut as cylinder compare)
        const Real az = std::abs(pos[i][2] - center[2]);
        if (az <= Real(0.15) * half_height)
        {
            if (r >= r_outer)
            {
                qn_o += v * static_cast<double>(q_num[i]);
                vn_o += v;
                qr_o += v * static_cast<double>(q_ref);
                vr_o += v;
            }
            if (r <= r_core)
            {
                qn_c += v * static_cast<double>(q_num[i]);
                vn_c += v;
                qr_c += v * static_cast<double>(q_ref);
                vr_c += v;
            }
        }
        (void)a;
    }

    rep.p_j_sq_over_2sigma = static_cast<Real>(p_jsq);
    rep.vol_continuous = Pi * radius * radius * height;
    rep.moment_r2_continuous = (Pi / Real(2)) * height * radius * radius * radius * radius;
    rep.p_ref_continuous =
        ophelieCylinderPrescribedANoSelfPowerPerLength(radius, b0, omega, sigma) * height;

    rep.field_rel = std::abs(rep.p_num_particle - rep.p_ref_discrete) / (std::abs(rep.p_ref_discrete) + TinyReal);
    rep.geometry_rel =
        std::abs(rep.p_ref_discrete - rep.p_ref_continuous) / (std::abs(rep.p_ref_continuous) + TinyReal);
    rep.total_vs_cont_rel =
        std::abs(rep.p_num_particle - rep.p_ref_continuous) / (std::abs(rep.p_ref_continuous) + TinyReal);
    rep.edge_vs_particle_rel =
        std::abs(rep.p_edge_recon - rep.p_num_particle) / (std::abs(rep.p_num_particle) + TinyReal);
    rep.j_rel_l2_vs_ref = (j_den > TinyReal) ? static_cast<Real>(std::sqrt(j_err / j_den)) : Real(0);

    auto ratio = [](double qo, double vo, double qc, double vc) -> Real {
        if (vo <= TinyReal || vc <= TinyReal)
            return Real(0);
        return static_cast<Real>((qo / vo) / (qc / vc + TinyReal));
    };
    rep.q_outer_over_core_num = ratio(qn_o, vn_o, qn_c, vn_c);
    rep.q_outer_over_core_ref_disc = ratio(qr_o, vr_o, qr_c, vr_c);
    // continuous area-weighted: Q∝r² ⇒ outer/core mean Q ratio
    {
        double qo = 0, vo = 0, qc = 0, vc = 0;
        const int nq = 400;
        for (int i = 0; i < nq; ++i)
        {
            const Real r0 = radius * Real(i) / Real(nq);
            const Real r1 = radius * Real(i + 1) / Real(nq);
            const Real rm = Real(0.5) * (r0 + r1);
            const Real area = Pi * (r1 * r1 - r0 * r0);
            const auto s = sampleOphelieCylinderPrescribedANoSelf(rm, radius, b0, omega, sigma);
            if (rm >= r_outer)
            {
                qo += s.Q_Wpm3 * area;
                vo += area;
            }
            if (rm <= r_core)
            {
                qc += s.Q_Wpm3 * area;
                vc += area;
            }
        }
        rep.q_outer_over_core_cont = ratio(qo, vo, qc, vc);
    }

    if (!particle_csv_path.empty())
    {
        std::ofstream csv(particle_csv_path);
        csv << "i,x,y,z,r,Vol,Q_num,Q_ref,|J_num|,|J_ref|\n";
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - center[0];
            const Real dy = pos[i][1] - center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Vecd a_ref = ophelieRotationalAUniformBz(pos[i], center, b0);
            const Real j_ref_n = sigma * omega * a_ref.norm();
            const Real q_ref = (j_ref_n * j_ref_n) / (Real(2) * sigma + TinyReal);
            const Real jn = std::sqrt(jr[i].squaredNorm() + ji[i].squaredNorm());
            csv << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << "," << vol[i]
                << "," << q_num[i] << "," << q_ref << "," << jn << "," << j_ref_n << "\n";
        }
    }

    std::cout << "[ophelie][doc199][p1] P_num=" << rep.p_num_particle << " P_ref_disc=" << rep.p_ref_discrete
              << " P_ref_cont=" << rep.p_ref_continuous << " P_edge_graph=" << rep.p_edge_graph
              << " P_edge_recon=" << rep.p_edge_recon << " P_J2/2sig=" << rep.p_j_sq_over_2sigma << std::endl;
    std::cout << "[ophelie][doc199][p1] vol_part=" << rep.vol_particles << " vol_cont=" << rep.vol_continuous
              << " sum_V_r2=" << rep.sum_v_r2 << " moment_r2_cont=" << rep.moment_r2_continuous << std::endl;
    std::cout << "[ophelie][doc199][p1] split field_rel=" << rep.field_rel << " geometry_rel=" << rep.geometry_rel
              << " total_vs_cont_rel=" << rep.total_vs_cont_rel
              << " edge_vs_particle_rel=" << rep.edge_vs_particle_rel << " J_rel_L2=" << rep.j_rel_l2_vs_ref
              << std::endl;
    std::cout << "[ophelie][doc199][p1] Q_o/c num=" << rep.q_outer_over_core_num
              << " ref_disc=" << rep.q_outer_over_core_ref_disc << " cont=" << rep.q_outer_over_core_cont
              << std::endl;
    // Informational gates (do not redefine overall PASS as geometry-only).
    const bool field_ok = rep.field_rel <= Real(0.01);
    const bool geom_note = rep.geometry_rel <= Real(0.05);
    const bool total_old = rep.total_vs_cont_rel <= Real(0.05);
    std::cout << "[ophelie][doc199][p1] gates field_rel<=1%=" << (field_ok ? 1 : 0)
              << " geometry_rel<=5%=" << (geom_note ? 1 : 0) << " total_vs_cont<=5%=" << (total_old ? 1 : 0)
              << " note: split reported even if total fails" << std::endl;
    if (!particle_csv_path.empty())
        std::cout << "[ophelie][doc199][p1] particle_csv=" << particle_csv_path << std::endl;
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC199_P1_POWER_DECOMPOSITION_H
