/**
 * @file electromagnetic_ophelie_doc215_e1d_zone_residual.h
 * @brief Doc215 E1d: per-particle / region residual of no-self J vs prescribed-A analytic,
 *        plus optional PhiGradLinearCorrection det (surface vs bulk diagnosis).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC215_E1D_ZONE_RESIDUAL_H
#define ELECTROMAGNETIC_OPHELIE_DOC215_E1D_ZONE_RESIDUAL_H

#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h"
#include "electromagnetic_ophelie_field_names.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc215E1dZoneAcc
{
    size_t n = 0;
    double w = 0;
    double se_j = 0;   // Σ V |J_sph−J_an|²
    double sn_j = 0;   // Σ V |J_an|²
    double sum_rel = 0;
    double max_rel = 0;
    size_t i_max = 0;
    double sum_det = 0;
    double min_det = 1e300;
    double max_det = -1e300;
    size_t n_det = 0;
};

inline void ophelieDoc215E1dAcc(OphelieDoc215E1dZoneAcc &z, size_t i, double vol, double dj2, double jan2,
                               double rel, double det, bool have_det)
{
    z.n++;
    z.w += vol;
    z.se_j += vol * dj2;
    z.sn_j += vol * jan2;
    z.sum_rel += rel;
    if (rel > z.max_rel)
    {
        z.max_rel = rel;
        z.i_max = i;
    }
    if (have_det)
    {
        z.n_det++;
        z.sum_det += det;
        z.min_det = std::min(z.min_det, det);
        z.max_det = std::max(z.max_det, det);
    }
}

inline Real ophelieDoc215E1dL2(const OphelieDoc215E1dZoneAcc &z)
{
    return (z.sn_j > TinyReal) ? static_cast<Real>(std::sqrt(z.se_j / z.sn_j)) : Real(0);
}

/**
 * Dump particles + zone summary. Analytic: prescribed_A_no_self J = (0, −σω A_rot).
 * Regions: Doc199 E2 tags interior|rim|end (+ all).
 * Optional radial shells by r/R and |z|/Hh mid vs endband.
 */
inline void hostOphelieDoc215E1dZoneResidualDump(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                 const Vecd &center, Real radius, Real half_height, Real b0,
                                                 Real omega, Real sigma, const std::string &csv_prefix,
                                                 bool expect_corrected_g)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Real>(particles, names.phi_imag);

    syncVariableToHost<Matd>(particles, names.phi_grad_linear_correction);
    Matd *B = particles.getVariableDataByName<Matd>(names.phi_grad_linear_correction);
    const bool have_B = (B != nullptr);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
    const Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);

    OphelieDoc215E1dZoneAcc z_all, z_int, z_rim, z_end;
    OphelieDoc215E1dZoneAcc z_r0, z_r1, z_r2; // r/R: <0.5, 0.5–0.75, ≥0.75
    OphelieDoc215E1dZoneAcc z_zmid, z_zend;   // |z|/Hh <0.5 vs ≥0.5

    const std::string part_path = csv_prefix + ".particles.csv";
    const std::string zone_path = csv_prefix + ".zones.csv";
    std::ofstream pofs(part_path);
    pofs << "i,x,y,z,r_over_R,abs_z_over_Hh,region,Vol,J_sph_abs,J_an_abs,rel_J,phi_r,phi_i,B_det\n";
    pofs << std::scientific << std::setprecision(10);

    for (size_t i = 0; i < n; ++i)
    {
        const Real xr = pos[i][0] - center[0];
        const Real yr = pos[i][1] - center[1];
        const Real zr = pos[i][2] - center[2];
        const Real r = std::sqrt(xr * xr + yr * yr);
        const Real rR = r / (radius + TinyReal);
        const Real zH = std::abs(zr) / (half_height + TinyReal);
        const char *tag = ophelieDoc199E2RegionTag(pos[i], center, radius, half_height);

        const Vecd a_an = ophelieRotationalAUniformBz(pos[i], center, b0);
        const Vecd j_an_r = Vecd::Zero();
        const Vecd j_an_i = -sigma * omega * a_an;
        const Vecd dj_r = jr[i] - j_an_r;
        const Vecd dj_i = ji[i] - j_an_i;
        const double jan2 = static_cast<double>(j_an_r.squaredNorm() + j_an_i.squaredNorm());
        const double jsp2 = static_cast<double>(jr[i].squaredNorm() + ji[i].squaredNorm());
        const double dj2 = static_cast<double>(dj_r.squaredNorm() + dj_i.squaredNorm());
        const double rel = std::sqrt(dj2) / (std::sqrt(jan2) + 1e-30);
        const double j_sph_abs = std::sqrt(jsp2);
        const double j_an_abs = std::sqrt(jan2);

        double det = 0.0;
        bool have_det = false;
        if (have_B)
        {
            det = static_cast<double>(B[i].determinant());
            have_det = std::isfinite(det);
        }

        const double v = static_cast<double>(vol[i]);
        ophelieDoc215E1dAcc(z_all, i, v, dj2, jan2, rel, det, have_det);
        if (std::strcmp(tag, "interior") == 0)
            ophelieDoc215E1dAcc(z_int, i, v, dj2, jan2, rel, det, have_det);
        else if (std::strcmp(tag, "rim") == 0)
            ophelieDoc215E1dAcc(z_rim, i, v, dj2, jan2, rel, det, have_det);
        else if (std::strcmp(tag, "end") == 0)
            ophelieDoc215E1dAcc(z_end, i, v, dj2, jan2, rel, det, have_det);

        if (rR < Real(0.5))
            ophelieDoc215E1dAcc(z_r0, i, v, dj2, jan2, rel, det, have_det);
        else if (rR < Real(0.75))
            ophelieDoc215E1dAcc(z_r1, i, v, dj2, jan2, rel, det, have_det);
        else
            ophelieDoc215E1dAcc(z_r2, i, v, dj2, jan2, rel, det, have_det);

        if (zH < Real(0.5))
            ophelieDoc215E1dAcc(z_zmid, i, v, dj2, jan2, rel, det, have_det);
        else
            ophelieDoc215E1dAcc(z_zend, i, v, dj2, jan2, rel, det, have_det);

        pofs << i << "," << double(pos[i][0]) << "," << double(pos[i][1]) << "," << double(pos[i][2]) << ","
             << double(rR) << "," << double(zH) << "," << tag << "," << double(vol[i]) << "," << j_sph_abs << ","
             << j_an_abs << "," << rel << "," << double(phi_r[i]) << "," << double(phi_i[i]) << ","
             << (have_det ? det : 0.0) << "\n";
    }
    pofs.close();

    std::ofstream zofs(zone_path);
    zofs << "zone,n,L2_rel_J,mean_rel_J,max_rel_J,i_max,B_det_mean,B_det_min,B_det_max,n_det\n";
    auto dumpz = [&](const char *name, const OphelieDoc215E1dZoneAcc &z) {
        const double mean_rel = (z.n > 0) ? (z.sum_rel / double(z.n)) : 0.0;
        const double det_mean = (z.n_det > 0) ? (z.sum_det / double(z.n_det)) : 0.0;
        const double det_min = (z.n_det > 0) ? z.min_det : 0.0;
        const double det_max = (z.n_det > 0) ? z.max_det : 0.0;
        zofs << name << "," << z.n << "," << double(ophelieDoc215E1dL2(z)) << "," << mean_rel << "," << z.max_rel
             << "," << z.i_max << "," << det_mean << "," << det_min << "," << det_max << "," << z.n_det << "\n";
        std::cout << std::scientific << std::setprecision(6) << "[ophelie][doc215][e1d] zone=" << name
                  << " n=" << z.n << " L2_rel_J=" << double(ophelieDoc215E1dL2(z)) << " mean_rel=" << mean_rel
                  << " max_rel=" << z.max_rel << " B_det_mean=" << det_mean << " B_det_min=" << det_min << "\n";
    };
    dumpz("all", z_all);
    dumpz("interior", z_int);
    dumpz("rim", z_rim);
    dumpz("end", z_end);
    dumpz("r_lt_0p5R", z_r0);
    dumpz("r_0p5_0p75R", z_r1);
    dumpz("r_ge_0p75R", z_r2);
    dumpz("z_midband", z_zmid);
    dumpz("z_endband", z_zend);

    std::cout << "[ophelie][doc215][e1d] expect_corrected_G=" << (expect_corrected_g ? 1 : 0)
              << " have_B=" << (have_B ? 1 : 0) << " particles=" << part_path << " zones=" << zone_path << "\n";
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC215_E1D_ZONE_RESIDUAL_H
