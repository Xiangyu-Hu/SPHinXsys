/**
 * @file electromagnetic_ophelie_doc199_e3_gauge_transform.h
 * @brief Doc199 E3: A1=A0+∇χ, φ1=φ0-ωχ ⇒ J,Q invariant (exp(+iωt)).
 *
 * Baseline A0 = rotational uniform Bz. Compare stored baseline J0/Q0 to gauge run.
 * Gates: ΔJ_rel≤1%, ΔP_rel≤2%; also J0 vs analytic -iωσA0 ≤5% (informational if inject).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC199_E3_GAUGE_TRANSFORM_H
#define ELECTROMAGNETIC_OPHELIE_DOC199_E3_GAUGE_TRANSFORM_H

#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h"
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

struct OphelieDoc199E3GaugeReport
{
    std::string chi_name = "unknown";
    std::string mode = "unknown"; // inject | solve
    Real j0_vs_analytic_rel = 0.0;
    Real j1_vs_j0_rel = 0.0;
    Real j1_vs_j0_rel_interior = 0.0;
    Real j1_vs_j0_rel_rim = 0.0;
    Real p0 = 0.0;
    Real p1 = 0.0;
    Real p_rel_change = 0.0;
    Real phi1_plus_omega_chi_rms = 0.0; // should be ~const (gauge); report rms after mean subtract
    Real phi1_rms = 0.0;
    bool gate_dj = false;
    bool gate_dp = false;
    bool gate_j0_analytic = false;
};

inline Real ophelieDoc199E3RegionRelDiff(const std::vector<Real> &vol, const std::vector<Vecd> &jr0,
                                        const std::vector<Vecd> &ji0, const std::vector<Vecd> &jr1,
                                        const std::vector<Vecd> &ji1, const std::vector<char> &mask)
{
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < vol.size(); ++i)
    {
        if (!mask[i])
            continue;
        const double v = static_cast<double>(vol[i]);
        const Vecd djr = jr1[i] - jr0[i];
        const Vecd dji = ji1[i] - ji0[i];
        num += v * static_cast<double>(djr.squaredNorm() + dji.squaredNorm());
        den += v * static_cast<double>(jr0[i].squaredNorm() + ji0[i].squaredNorm());
    }
    return (den > TinyReal) ? static_cast<Real>(std::sqrt(num / den)) : Real(0);
}

inline OphelieDoc199E3GaugeReport hostOphelieDoc199E3GaugeAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieDoc199E2ChiField &chi, Real omega,
    Real b0, Real sigma, Real radius, Real half_height, const std::vector<Vecd> &jr0,
    const std::vector<Vecd> &ji0, const std::vector<Real> &q0, const std::string &mode,
    const std::string &particle_csv_path = {})
{
    OphelieDoc199E3GaugeReport rep;
    rep.chi_name = ophelieDoc199E2ChiKindName(chi.kind);
    rep.mode = mode;

    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol_ptr = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Vecd *jr1 = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji1 = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *q1 = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);

    std::vector<Real> vol(n);
    std::vector<Vecd> jr0v = jr0, ji0v = ji0, jr1v(n), ji1v(n);
    std::vector<char> mask_all(n, 1), mask_int(n, 0), mask_rim(n, 0);
    double p0 = 0.0, p1 = 0.0, vol_tot = 0.0;
    double j0_err = 0.0, j0_den = 0.0;
    double phi2 = 0.0, gauge2 = 0.0, gauge_sum = 0.0;
    std::vector<Real> gauge_phi(n);

    for (size_t i = 0; i < n; ++i)
    {
        vol[i] = vol_ptr[i];
        jr1v[i] = jr1[i];
        ji1v[i] = ji1[i];
        const double v = static_cast<double>(vol[i]);
        vol_tot += v;
        p0 += v * static_cast<double>(q0[i]);
        p1 += v * static_cast<double>(q1[i]);

        const Vecd a0 = ophelieRotationalAUniformBz(pos[i], chi.center, b0);
        const Vecd j_an_i = -sigma * omega * a0; // φ0=0 imag chain
        const Vecd dj0 = ji0[i] - j_an_i;
        const Vecd dj0r = jr0[i]; // analytic real J = 0
        j0_err += v * static_cast<double>(dj0r.squaredNorm() + dj0.squaredNorm());
        j0_den += v * static_cast<double>(j_an_i.squaredNorm());

        const Real chi_i = ophelieDoc199E2ChiAt(pos[i], chi);
        gauge_phi[i] = phi[i] + omega * chi_i; // φ1 + ωχ ~ φ0 + const
        gauge_sum += v * static_cast<double>(gauge_phi[i]);
        phi2 += v * static_cast<double>(phi[i] * phi[i]);

        const char *tag = ophelieDoc199E2RegionTag(pos[i], chi.center, radius, half_height);
        if (std::strcmp(tag, "interior") == 0)
            mask_int[i] = 1;
        if (std::strcmp(tag, "rim") == 0)
            mask_rim[i] = 1;

        (void)a;
    }

    const Real gauge_mean = (vol_tot > TinyReal) ? static_cast<Real>(gauge_sum / vol_tot) : Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        const double v = static_cast<double>(vol[i]);
        const Real d = gauge_phi[i] - gauge_mean;
        gauge2 += v * static_cast<double>(d * d);
    }

    rep.j0_vs_analytic_rel = (j0_den > TinyReal) ? static_cast<Real>(std::sqrt(j0_err / j0_den)) : Real(0);
    rep.j1_vs_j0_rel = ophelieDoc199E3RegionRelDiff(vol, jr0v, ji0v, jr1v, ji1v, mask_all);
    rep.j1_vs_j0_rel_interior = ophelieDoc199E3RegionRelDiff(vol, jr0v, ji0v, jr1v, ji1v, mask_int);
    rep.j1_vs_j0_rel_rim = ophelieDoc199E3RegionRelDiff(vol, jr0v, ji0v, jr1v, ji1v, mask_rim);
    rep.p0 = static_cast<Real>(p0);
    rep.p1 = static_cast<Real>(p1);
    rep.p_rel_change = std::abs(rep.p1 - rep.p0) / (std::abs(rep.p0) + TinyReal);
    rep.phi1_rms = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(phi2 / vol_tot)) : Real(0);
    rep.phi1_plus_omega_chi_rms =
        (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(gauge2 / vol_tot)) : Real(0);

    rep.gate_dj = rep.j1_vs_j0_rel <= Real(0.01);
    rep.gate_dp = rep.p_rel_change <= Real(0.02);
    rep.gate_j0_analytic = rep.j0_vs_analytic_rel <= Real(0.05);

    if (!particle_csv_path.empty())
    {
        std::ofstream csv(particle_csv_path);
        csv << "i,x,y,z,r,region,Vol,chi,phi1,gauge_phi,|J0|,|J1|,Q0,Q1,"
               "J0x_re,J0y_re,J0z_re,J0x_im,J0y_im,J0z_im,"
               "J1x_re,J1y_re,J1z_re,J1x_im,J1y_im,J1z_im\n";
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - chi.center[0];
            const Real dy = pos[i][1] - chi.center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real j0n = std::sqrt(jr0[i].squaredNorm() + ji0[i].squaredNorm());
            const Real j1n = std::sqrt(jr1[i].squaredNorm() + ji1[i].squaredNorm());
            csv << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << ","
                << ophelieDoc199E2RegionTag(pos[i], chi.center, radius, half_height) << "," << vol[i] << ","
                << ophelieDoc199E2ChiAt(pos[i], chi) << "," << phi[i] << "," << gauge_phi[i] << "," << j0n << ","
                << j1n << "," << q0[i] << "," << q1[i] << "," << jr0[i][0] << "," << jr0[i][1] << "," << jr0[i][2]
                << "," << ji0[i][0] << "," << ji0[i][1] << "," << ji0[i][2] << "," << jr1[i][0] << "," << jr1[i][1]
                << "," << jr1[i][2] << "," << ji1[i][0] << "," << ji1[i][1] << "," << ji1[i][2] << "\n";
        }
    }

    std::cout << "[ophelie][doc199][e3] mode=" << rep.mode << " chi=" << rep.chi_name
              << " J0_vs_analytic=" << rep.j0_vs_analytic_rel << " J1_vs_J0 all=" << rep.j1_vs_j0_rel
              << " interior=" << rep.j1_vs_j0_rel_interior << " rim=" << rep.j1_vs_j0_rel_rim << " P0=" << rep.p0
              << " P1=" << rep.p1 << " P_rel_change=" << rep.p_rel_change << " phi1_rms=" << rep.phi1_rms
              << " (phi1+omega*chi)_rms_demean=" << rep.phi1_plus_omega_chi_rms << std::endl;
    std::cout << "[ophelie][doc199][e3] gates dJ<=1%=" << (rep.gate_dj ? 1 : 0)
              << " dP<=2%=" << (rep.gate_dp ? 1 : 0) << " J0_analytic<=5%=" << (rep.gate_j0_analytic ? 1 : 0)
              << std::endl;
    if (!particle_csv_path.empty())
        std::cout << "[ophelie][doc199][e3] particle_csv=" << particle_csv_path << std::endl;
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC199_E3_GAUGE_TRANSFORM_H
