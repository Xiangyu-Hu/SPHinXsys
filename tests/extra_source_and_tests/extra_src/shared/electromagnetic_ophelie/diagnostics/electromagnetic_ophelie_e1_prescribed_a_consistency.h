/**
 * @file electromagnetic_ophelie_e1_prescribed_a_consistency.h
 * @brief Doc198 E1: prescribed rotating A, isolate edge-EMF / recon / phi layers.
 *
 * Phasor: exp(+iωt), E=-iωA-∇φ. With A real and φ=0: E_r=0, E_i=-ωA, J=σE.
 * Affine A=½B0(-y,x,0) → midpoint edge formula equals exact straight-line ∫A·dl.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_E1_PRESCRIBED_A_CONSISTENCY_H
#define ELECTROMAGNETIC_OPHELIE_E1_PRESCRIBED_A_CONSISTENCY_H

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

struct OphelieE1RegionStats
{
    double w = 0.0;
    double l2_num_j = 0.0;
    double l2_den_j = 0.0;
    double l2_num_e = 0.0;
    double l2_den_e = 0.0;
    double l2_num_j_nonaz = 0.0; // non-azimuthal |J|² weight
    size_t n = 0;
    size_t n_fallback = 0;
};

inline Real ophelieE1RelL2(const OphelieE1RegionStats &s)
{
    return (s.l2_den_j > TinyReal) ? static_cast<Real>(std::sqrt(s.l2_num_j / s.l2_den_j)) : Real(0);
}
inline Real ophelieE1RelL2E(const OphelieE1RegionStats &s)
{
    return (s.l2_den_e > TinyReal) ? static_cast<Real>(std::sqrt(s.l2_num_e / s.l2_den_e)) : Real(0);
}
inline Real ophelieE1NonAzFrac(const OphelieE1RegionStats &s)
{
    return (s.l2_den_j > TinyReal) ? static_cast<Real>(std::sqrt(s.l2_num_j_nonaz / s.l2_den_j)) : Real(0);
}

struct OphelieE1PrescribedAReport
{
    std::string layer_tag = "unknown";
    Real a_prescription_rel_l2 = 0.0;   // stored A vs analytic rotational A
    Real edge_affine_identity_max = 0.0; // |midpoint∫A - exact affine| / (|A|·|Δx|) proxy on sample
    Real phi_imag_vol_rms = 0.0;
    Real phi_real_vol_rms = 0.0;
    OphelieE1RegionStats all;
    OphelieE1RegionStats interior; // r<0.7R, |z-zc|<0.5 Hh
    OphelieE1RegionStats rim;      // r>=0.75R midband
    OphelieE1RegionStats axis_lat; // midband |θ| mod 90 <=15°
    OphelieE1RegionStats diag_lat; // midband |θ| mod 90 >=30°
    size_t n_particles = 0;
};

inline void ophelieE1Accumulate(OphelieE1RegionStats &s, Real vol, const Vecd &j_sph, const Vecd &j_an,
                                const Vecd &e_sph, const Vecd &e_an, const Vecd &j_nonaz, bool fallback)
{
    const double w = static_cast<double>(vol);
    const Vecd dj = j_sph - j_an;
    const Vecd de = e_sph - e_an;
    s.w += w;
    s.l2_num_j += w * static_cast<double>(dj.squaredNorm());
    s.l2_den_j += w * static_cast<double>(j_an.squaredNorm());
    s.l2_num_e += w * static_cast<double>(de.squaredNorm());
    s.l2_den_e += w * static_cast<double>(e_an.squaredNorm());
    s.l2_num_j_nonaz += w * static_cast<double>(j_nonaz.squaredNorm());
    ++s.n;
    if (fallback)
    {
        ++s.n_fallback;
    }
}

/**
 * Host E1 audit after edge-flux recon (complex). A_coil must be prescribed rotational field.
 * J_an = (0, -σω A) in (real,imag); E_an = (0, -ω A).
 */
inline OphelieE1PrescribedAReport hostOphelieE1PrescribedAConsistencyAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const Vecd &axis_origin, Real radius_m,
    Real half_height, Real b0, Real omega, Real sigma, const char *layer_tag)
{
    OphelieE1PrescribedAReport rep;
    rep.layer_tag = layer_tag;
    const size_t n = particles.TotalRealParticles();
    rep.n_particles = n;

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Real>(particles, names.edge_recon_fallback);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a_coil = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *er = particles.getVariableDataByName<Vecd>(names.e_edge_recon_real);
    const Vecd *ei = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
    const Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
    const Real *fb = particles.getVariableDataByName<Real>(names.edge_recon_fallback);

    double a_num = 0.0, a_den = 0.0;
    double phi_i2 = 0.0, phi_r2 = 0.0, w_tot = 0.0;
    const Real zc = axis_origin[2];
    const Real z_int = Real(0.5) * half_height;
    const Real z_mid = Real(0.15) * half_height;

    // Affine edge identity: sample consecutive particles is weak; use particle self-pair with
    // displaced synthetic neighbor offset (1e-3 R) as formula check only.
    Real edge_id_max = 0.0;
    for (size_t i = 0; i < n; i += std::max<size_t>(1, n / 200))
    {
        const Vecd a_i = ophelieRotationalAUniformBz(pos[i], axis_origin, b0);
        const Vecd dx(radius_m * Real(1.0e-3), radius_m * Real(2.0e-3), Real(0));
        const Vecd pj = pos[i] + dx;
        const Vecd a_j = ophelieRotationalAUniformBz(pj, axis_origin, b0);
        const Vecd a_avg = Real(0.5) * (a_i + a_j);
        const Real mid = a_avg.dot(dx);
        // Exact ∫ A·dl for A=½B0(-y,x,0) relative to axis: ½ B0 (x_i y_j - y_i x_j) in plane
        // (translation-invariant curl form). With axis origin subtracted:
        const Real xi = pos[i][0] - axis_origin[0];
        const Real yi = pos[i][1] - axis_origin[1];
        const Real xj = pj[0] - axis_origin[0];
        const Real yj = pj[1] - axis_origin[1];
        const Real exact = Real(0.5) * b0 * (xi * yj - yi * xj);
        const Real scale = (a_i.norm() + TinyReal) * dx.norm();
        edge_id_max = std::max(edge_id_max, std::abs(mid - exact) / (scale + TinyReal));
    }
    rep.edge_affine_identity_max = edge_id_max;

    for (size_t i = 0; i < n; ++i)
    {
        const Real w = vol[i];
        w_tot += w;
        phi_i2 += w * static_cast<double>(phi_i[i] * phi_i[i]);
        phi_r2 += w * static_cast<double>(phi_r[i] * phi_r[i]);

        const Vecd a_an = ophelieRotationalAUniformBz(pos[i], axis_origin, b0);
        const Vecd da = a_coil[i] - a_an;
        a_num += w * static_cast<double>(da.squaredNorm());
        a_den += w * static_cast<double>(a_an.squaredNorm());

        // Analytic E/J (φ=0 contract)
        const Vecd e_an_r = Vecd::Zero();
        const Vecd e_an_i = -omega * a_an;
        const Vecd j_an_r = Vecd::Zero();
        const Vecd j_an_i = sigma * e_an_i;

        const Vecd j_sph_r = jr[i];
        const Vecd j_sph_i = ji[i];
        // Pack complex error via stacked real+imag into one Vecd metric: use sum of squares
        const Vecd dj_r = j_sph_r - j_an_r;
        const Vecd dj_i = j_sph_i - j_an_i;
        const Vecd de_r = er[i] - e_an_r;
        const Vecd de_i = ei[i] - e_an_i;

        // Fake "vector" norms by accumulating into region with combined weights
        OphelieE1RegionStats tmp{};
        const double ww = static_cast<double>(w);
        tmp.l2_num_j = ww * (static_cast<double>(dj_r.squaredNorm()) + static_cast<double>(dj_i.squaredNorm()));
        tmp.l2_den_j = ww * (static_cast<double>(j_an_r.squaredNorm()) + static_cast<double>(j_an_i.squaredNorm()));
        tmp.l2_num_e = ww * (static_cast<double>(de_r.squaredNorm()) + static_cast<double>(de_i.squaredNorm()));
        tmp.l2_den_e = ww * (static_cast<double>(e_an_r.squaredNorm()) + static_cast<double>(e_an_i.squaredNorm()));

        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        // Azimuthal unit e_θ = (-y/r, x/r, 0); radial e_r=(x/r,y/r,0)
        const Real er_x = inv_r * dx;
        const Real er_y = inv_r * dy;
        const Real eth_x = -inv_r * dy;
        const Real eth_y = inv_r * dx;
        auto nonaz = [&](const Vecd &jv) {
            const Real jth = jv[0] * eth_x + jv[1] * eth_y;
            const Real jr = jv[0] * er_x + jv[1] * er_y;
            const Real jz = jv[2];
            // residual non-azimuthal = full - azimuthal-only reconstruction of imag/real separately handled below
            (void)jth;
            return Vecd(jr * er_x, jr * er_y, jz); // radial + z parts
        };
        // Non-az of SPH J (real+imag): measure |J_r| and |J_z|
        const Vecd naz_r = nonaz(j_sph_r);
        const Vecd naz_i = nonaz(j_sph_i);
        tmp.l2_num_j_nonaz = ww * (static_cast<double>(naz_r.squaredNorm()) + static_cast<double>(naz_i.squaredNorm()));

        auto add = [&](OphelieE1RegionStats &s) {
            s.w += ww;
            s.l2_num_j += tmp.l2_num_j;
            s.l2_den_j += tmp.l2_den_j;
            s.l2_num_e += tmp.l2_num_e;
            s.l2_den_e += tmp.l2_den_e;
            s.l2_num_j_nonaz += tmp.l2_num_j_nonaz;
            ++s.n;
            if (fb[i] > Real(0.5))
                ++s.n_fallback;
        };
        add(rep.all);

        const Real abs_z = std::abs(pos[i][2] - zc);
        if (r < Real(0.70) * radius_m && abs_z < z_int)
        {
            add(rep.interior);
        }
        if (r >= Real(0.75) * radius_m && abs_z < z_mid)
        {
            add(rep.rim);
        }
        if (abs_z < z_mid && r < radius_m)
        {
            Real th = static_cast<Real>(std::atan2(dy, dx) * Real(180) / Pi);
            th = std::abs(th);
            th = std::fmod(th, Real(90));
            th = std::min(th, Real(90) - th);
            if (th <= Real(15))
                add(rep.axis_lat);
            if (th >= Real(30))
                add(rep.diag_lat);
        }
    }

    rep.a_prescription_rel_l2 =
        (a_den > TinyReal) ? static_cast<Real>(std::sqrt(a_num / a_den)) : Real(0);
    rep.phi_imag_vol_rms = (w_tot > TinyReal) ? static_cast<Real>(std::sqrt(phi_i2 / w_tot)) : Real(0);
    rep.phi_real_vol_rms = (w_tot > TinyReal) ? static_cast<Real>(std::sqrt(phi_r2 / w_tot)) : Real(0);
    return rep;
}

inline void logOphelieE1PrescribedAReport(const OphelieE1PrescribedAReport &r)
{
    auto one = [](const char *name, const OphelieE1RegionStats &s) {
        std::cout << " " << name << "{n=" << s.n << " fb=" << s.n_fallback << " L2_J=" << ophelieE1RelL2(s)
                  << " L2_E=" << ophelieE1RelL2E(s) << " nonaz/J=" << ophelieE1NonAzFrac(s) << "}";
    };
    std::cout << "[ophelie][e1] layer=" << r.layer_tag << " n=" << r.n_particles
              << " A_prescription_L2=" << r.a_prescription_rel_l2
              << " edge_affine_id_max=" << r.edge_affine_identity_max
              << " phi_imag_rms=" << r.phi_imag_vol_rms << " phi_real_rms=" << r.phi_real_vol_rms;
    one("all", r.all);
    one("interior", r.interior);
    one("rim", r.rim);
    one("axis", r.axis_lat);
    one("diag", r.diag_lat);
    std::cout << std::endl;
}

inline bool hostOphelieE1DumpParticleCsv(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                         const Vecd &axis_origin, Real radius_m, Real half_height, Real b0,
                                         Real omega, Real sigma, const std::string &csv_path)
{
    if (csv_path.empty())
        return false;
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Real>(particles, names.edge_recon_fallback);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
    const Real *fb = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
    std::ofstream ofs(csv_path);
    if (!ofs)
        return false;
    ofs << "x,y,z,r,theta_deg,vol,Jth_sph,Jth_an,signed_rel_Jth,phi_imag,fallback\n";
    const Real zc = axis_origin[2];
    const Real z_mid = Real(0.15) * half_height;
    size_t n_out = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        if (r >= radius_m || std::abs(pos[i][2] - zc) > z_mid)
            continue;
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const Real jt_re = inv_r * (-jr[i][0] * dy + jr[i][1] * dx);
        const Real jt_im = inv_r * (-ji[i][0] * dy + ji[i][1] * dx);
        const Real jt = std::hypot(jt_re, jt_im);
        const Vecd a_an = ophelieRotationalAUniformBz(pos[i], axis_origin, b0);
        const Real a_th = inv_r * (-a_an[0] * dy + a_an[1] * dx); // = ½ B0 r
        const Real j_an = sigma * omega * std::abs(a_th);
        const Real signed_rel = (j_an > Real(1e-9)) ? (jt - j_an) / j_an : Real(0);
        const Real th = static_cast<Real>(std::atan2(dy, dx) * Real(180) / Pi);
        ofs << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << "," << th << "," << vol[i] << ","
            << jt << "," << j_an << "," << signed_rel << "," << phi_i[i] << "," << fb[i] << "\n";
        ++n_out;
    }
    std::cout << "[ophelie][e1] particle_csv=" << csv_path << " n=" << n_out << std::endl;
    return n_out > 0;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_E1_PRESCRIBED_A_CONSISTENCY_H
