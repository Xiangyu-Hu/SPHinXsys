/**
 * @file electromagnetic_ophelie_doc199_edge_continuity_audit.h
 * @brief Doc199 E1: after phi solve, separate linear eq residual vs edge-current continuity.
 *
 * Edge current density (imag chain): q_ij = -C_ij * [(φ_j-φ_i) + a_sign ω A_avg·(x_j-x_i)].
 * Continuity residual at particle i: r_edge_i = Σ_j q_ij.
 * Linear residual: r_lin_i = (Lφ)_i - b_i with production LHS/RHS operators.
 *
 * Does NOT zero-mean or project residuals. Global Σ V r_edge ≈ 0 is NOT a pass criterion.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC199_EDGE_CONTINUITY_AUDIT_H
#define ELECTROMAGNETIC_OPHELIE_DOC199_EDGE_CONTINUITY_AUDIT_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_gradient.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_phi_solvability.h"

#include "base_body_relation.h"
#include "inner_body_relation.h"
#include "update_body_relation.h"

#include <algorithm>
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

struct OphelieDoc199RegionRes
{
    size_t n = 0;
    double vol = 0.0;
    double l2_num = 0.0; // Σ V r²
    double scale_num = 0.0; // Σ V s²  (local scale)
    double max_abs = 0.0;
};

inline void ophelieDoc199AccRegion(OphelieDoc199RegionRes &s, Real v, Real r, Real scale)
{
    const double vv = static_cast<double>(v);
    s.vol += vv;
    s.l2_num += vv * static_cast<double>(r) * static_cast<double>(r);
    s.scale_num += vv * static_cast<double>(scale) * static_cast<double>(scale);
    s.max_abs = std::max(s.max_abs, static_cast<double>(std::abs(r)));
    ++s.n;
}

inline Real ophelieDoc199RelL2(const OphelieDoc199RegionRes &s)
{
    return (s.scale_num > TinyReal) ? static_cast<Real>(std::sqrt(s.l2_num / s.scale_num)) : Real(0);
}

inline Real ophelieDoc199AbsL2(const OphelieDoc199RegionRes &s)
{
    return (s.vol > TinyReal) ? static_cast<Real>(std::sqrt(s.l2_num / s.vol)) : Real(0);
}

struct OphelieDoc199EdgeContinuityReport
{
    std::string phi_rhs_name = "unknown";
    std::string phi_lhs_name = "unknown";
    Real phi_eq_res_vol = 0.0;          // ||Lφ-b||_vol / ||b||_vol
    Real edge_res_l2 = 0.0;            // abs L2 of Σq
    Real edge_res_rel_l2 = 0.0;        // vs local |Σ C ω A·Δx| scale
    Real edge_global_integral = 0.0;   // Σ V Σq (diagnostic only)
    Real antisym_I_rel_rms = 0.0;      // rms(|I_ij+I_ji|) / rms(|I_ij|)
    Real rhs_div_vs_edge_src_rel = 0.0; // DivSigmaA b vs φ=0 edge-source Σ(-C ω A·Δx)
    Real rhs_legacy_vs_edge_src_rel = 0.0;
    Real div_j_recon_abs_l2 = 0.0; // ||div J_edge_recon_imag||_vol (legacy abs)
    Real div_j_recon_rel_l2 = 0.0; // LEGACY unit-bearing: ||div||/||σω|A||| (Doc200: do not treat as %)
    // Doc200 P1 dimensionless: ε = L* * sqrt(Σ V |div J_c|² / Σ V) / J*
    Real L_star = 0.0;
    Real J_star = 0.0;
    Real epsilon_div_complex = 0.0; // preferred continuity metric (edge-recon J)
    Real div_j_complex_abs_l2 = 0.0; // sqrt(Σ V (|div Jr|²+|div Ji|²)/Σ V)
    /**
     * Doc201 A4 dual column: J_σE from uncorrected Gradφ + A (SPH family matching DivΣA identity),
     * then same uncorrected div + (L*,J*) as epsilon_div_complex. Does NOT replace edge-recon gate.
     */
    Real epsilon_div_dual_gradphi = 0.0;
    Real j_dual_vs_recon_rel = 0.0; // vol-weighted |J_dual - J_recon| / |J_recon|
    // Doc212: region split of the same ε_div measure (edge-recon complex |div J|).
    Real epsilon_div_interior = 0.0;
    Real epsilon_div_rim = 0.0;
    Real recon_fallback_fraction = 0.0;
    Real recon_fallback_fraction_rim = 0.0;
    /**
     * Doc212 pairing column: same (L*,J*) nondim as epsilon_div_complex, but on edge Kirchhoff
     * residual r = Σq (imag+real). NOT a substitute gate — proves SPH-div vs edge-Σq split.
     */
    Real epsilon_div_edge_kirchhoff = 0.0;
    Real sph_vs_edge_eps_ratio = 0.0; // epsilon_div_complex / max(epsilon_div_edge_kirchhoff, Tiny)
    Real corr_edge_res_vs_sph_div_im = 0.0; // vol-weighted corr(|r_edge_im|, |div J_im|)
    OphelieDoc199RegionRes all_edge;
    OphelieDoc199RegionRes interior_edge;
    OphelieDoc199RegionRes rim_edge;
    size_t n_edges_unique = 0;
    size_t n_particles = 0;
};

inline const char *ophelieDoc199RegionTag(const Vecd &pos, const Vecd &center, Real radius, Real half_height)
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
    return "other";
}

/**
 * Doc199 E1 continuity audit. Requires post-phi fields already on particles.
 * Recomputes production RHS + LHS at current φ; evaluates edge Σq; dumps CSVs.
 */
template <class ExecutionPolicy>
inline OphelieDoc199EdgeContinuityReport execOphelieDoc199EdgeContinuityAudit(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters params,
    const Vecd &center, Real radius, Real half_height, const std::string &particle_csv_path,
    const std::string &edge_csv_path)
{
    OphelieDoc199EdgeContinuityReport rep;
    rep.phi_rhs_name = phiRhsOperatorKindName(params.phi_rhs_operator_kind_);
    rep.phi_lhs_name =
        (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::LegacyPairwise) ? "legacy-pairwise" : "div-sigma-grad";

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    rep.n_particles = n;

    // --- Linear residual at production operators ---
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    rep.phi_eq_res_vol = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);

    std::vector<Real> rhs_prod(n, Real(0));
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_prod.data(), n);
    std::vector<Real> lhs_prod(n, Real(0));
    hostReadScalarField(particles, names.phi_lhs_imag, lhs_prod.data(), n);

    // --- Edge residual Σ q_ij (uses final φ) ---
    const OphelieEdgeFluxComponent imag = makeOphelieEdgeFluxImagComponent(names, params);
    const OphelieEdgeFluxResidualMetrics edge_m =
        evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, inner, names, imag, params);
    rep.edge_res_l2 = edge_m.edge_res_l2;
    rep.edge_global_integral = edge_m.edge_res_volume_integral;

    syncVariableToHost<Real>(particles, imag.edge_residual_field);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, imag.active_a_field);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);

    const Real *edge_res = particles.getVariableDataByName<Real>(imag.edge_residual_field);
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Vecd *a = particles.getVariableDataByName<Vecd>(imag.active_a_field);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *j_im = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *j_re = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);

    // Doc212: keep imag Σq, then evaluate real-chain Σq for complex edge-Kirchhoff ε.
    std::vector<Real> edge_res_im_host(n, Real(0));
    std::vector<Real> edge_res_re_host(n, Real(0));
    for (size_t i = 0; i < n; ++i)
        edge_res_im_host[i] = edge_res[i];
    {
        const OphelieEdgeFluxComponent real_c = makeOphelieEdgeFluxRealComponent(names, params);
        (void)evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, inner, names, real_c, params);
        syncVariableToHost<Real>(particles, real_c.edge_residual_field);
        const Real *edge_res_re = particles.getVariableDataByName<Real>(real_c.edge_residual_field);
        for (size_t i = 0; i < n; ++i)
            edge_res_re_host[i] = edge_res_re[i];
        // Restore imag residual field pointer target for later CSV (real eval wrote real field only).
        syncVariableToHost<Real>(particles, imag.edge_residual_field);
        edge_res = particles.getVariableDataByName<Real>(imag.edge_residual_field);
    }

    // Classic InnerRelation for unique-edge dump + local source scale
    InnerRelation classic(glass_body);
    classic.updateConfiguration();

    const Real omega = params.omega();
    const Real a_sign = imag.a_sign;
    const Real h = glass_body.getSPHAdaptation().ReferenceSmoothingLength();
    const Real pair_reg = params.pair_weight_regularization_;

    std::vector<Real> edge_src_scale(n, Real(0)); // |Σ -C a_sign ω A·Δx|  (φ=0 source)
    std::vector<Real> edge_src_signed(n, Real(0));
    double antisym_sq = 0.0, i_scale_sq = 0.0;
    size_t n_unique = 0;

    std::ofstream edge_csv;
    if (!edge_csv_path.empty())
    {
        edge_csv.open(edge_csv_path);
        edge_csv << "edge_id,i,j,xi,yi,zi,xj,yj,zj,r_ij,C_ij,dphi,emf_A,edge_drop,q_ij,I_ij,I_ji,antisym_I\n";
    }

    for (size_t i = 0; i < n; ++i)
    {
        const Neighborhood &nb = classic.inner_configuration_[i];
        Real src_i = 0.0;
        for (size_t k = 0; k < nb.current_size_; ++k)
        {
            const size_t j = nb.j_[k];
            const Real distance = nb.r_ij_[k];
            const Real distance_sq = distance * distance;
            const Vecd r_ij = nb.e_ij_[k] * distance; // e_ij points j→i? check: usually e_ij from i to j
            // SPHinXsys Neighborhood: e_ij unit vector; r_ij distance. Convention in CK: vec_r_ij = pos_i-pos_j.
            // classic e_ij_ is typically (pos_i - pos_j)/r = pointing from j to i.
            // Edge flux CK: r_ij_vec = vec_r_ij(i,j) = pos_i - pos_j; xj_minus_xi = -r_ij_vec = pos_j - pos_i.
            const Vecd r_ij_vec = r_ij; // if e from j to i, this is pos_i-pos_j — matches CK.
            const Vecd xj_minus_xi = -r_ij_vec;
            const Real dW_ijV_j = nb.dW_ij_[k] * vol[j];
            const Real c_ij = computeOphelieEdgeFluxPairConductance(sigma[i], sigma[j], dW_ijV_j, distance, distance_sq,
                                                                    pair_reg, h);
            const Real emf_A = a_sign * omega * (Real(0.5) * (a[i] + a[j])).dot(xj_minus_xi);
            const Real dphi = phi[j] - phi[i];
            const Real edge_drop = dphi + emf_A;
            const Real q_ij = -c_ij * edge_drop;
            src_i += -c_ij * emf_A;

            // Antisymmetry / unique edge (i < j): C_ji from j's neighborhood view.
            if (i < j)
            {
                Real c_ji = c_ij;
                const Neighborhood &nb_j = classic.inner_configuration_[j];
                for (size_t kj = 0; kj < nb_j.current_size_; ++kj)
                {
                    if (nb_j.j_[kj] != i)
                        continue;
                    const Real dW_jiV_i = nb_j.dW_ij_[kj] * vol[i];
                    c_ji = computeOphelieEdgeFluxPairConductance(sigma[j], sigma[i], dW_jiV_i, distance, distance_sq,
                                                                 pair_reg, h);
                    break;
                }
                const Real edge_drop_ji =
                    computeOphelieEdgeFluxEdgeDrop(phi[j], phi[i], a[j], a[i], r_ij_vec, omega, a_sign);
                const Real q_ji = -c_ji * edge_drop_ji;
                const Real I_ij = vol[i] * q_ij;
                const Real I_ji = vol[j] * q_ji;
                const Real antisym = I_ij + I_ji;
                antisym_sq += static_cast<double>(antisym) * static_cast<double>(antisym);
                i_scale_sq += static_cast<double>(I_ij) * static_cast<double>(I_ij);
                ++n_unique;
                if (edge_csv)
                {
                    edge_csv << n_unique << "," << i << "," << j << "," << pos[i][0] << "," << pos[i][1] << ","
                             << pos[i][2] << "," << pos[j][0] << "," << pos[j][1] << "," << pos[j][2] << "," << distance
                             << "," << c_ij << "," << dphi << "," << emf_A << "," << edge_drop << "," << q_ij << ","
                             << I_ij << "," << I_ji << "," << antisym << "\n";
                }
            }
        }
        edge_src_signed[i] = src_i;
        edge_src_scale[i] = std::abs(src_i);
    }
    rep.n_edges_unique = n_unique;
    rep.antisym_I_rel_rms =
        (i_scale_sq > TinyReal) ? static_cast<Real>(std::sqrt(antisym_sq / i_scale_sq)) : Real(0);

    // Region edge continuity
    for (size_t i = 0; i < n; ++i)
    {
        const char *tag = ophelieDoc199RegionTag(pos[i], center, radius, half_height);
        const Real scale = std::max(edge_src_scale[i], Real(TinyReal));
        ophelieDoc199AccRegion(rep.all_edge, vol[i], edge_res[i], scale);
        if (std::strcmp(tag, "interior") == 0)
            ophelieDoc199AccRegion(rep.interior_edge, vol[i], edge_res[i], scale);
        if (std::strcmp(tag, "rim") == 0)
            ophelieDoc199AccRegion(rep.rim_edge, vol[i], edge_res[i], scale);
    }
    rep.edge_res_rel_l2 = ophelieDoc199RelL2(rep.all_edge);

    // Compare production DivSigmaA / LegacyFlux RHS to φ=0 edge source
    {
        OphelieParameters p_div = params;
        p_div.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
        setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, p_div);
        std::vector<Real> rhs_div(n, Real(0));
        hostReadScalarField(particles, names.phi_rhs_imag, rhs_div.data(), n);

        OphelieParameters p_leg = params;
        p_leg.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
        setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, p_leg);
        std::vector<Real> rhs_leg(n, Real(0));
        hostReadScalarField(particles, names.phi_rhs_imag, rhs_leg.data(), n);

        rep.rhs_div_vs_edge_src_rel =
            hostVolWeightedRelativeDifference(particles, edge_src_signed.data(), rhs_div.data(), n);
        // Note: edge_src and RHS may differ by sign convention; also report flipped.
        {
            std::vector<Real> neg(n);
            for (size_t i = 0; i < n; ++i)
                neg[i] = -edge_src_signed[i];
            const Real flip = hostVolWeightedRelativeDifference(particles, neg.data(), rhs_div.data(), n);
            if (flip < rep.rhs_div_vs_edge_src_rel)
                rep.rhs_div_vs_edge_src_rel = flip;
        }
        rep.rhs_legacy_vs_edge_src_rel =
            hostVolWeightedRelativeDifference(particles, edge_src_signed.data(), rhs_leg.data(), n);
        {
            std::vector<Real> neg(n);
            for (size_t i = 0; i < n; ++i)
                neg[i] = -edge_src_signed[i];
            const Real flip = hostVolWeightedRelativeDifference(particles, neg.data(), rhs_leg.data(), n);
            if (flip < rep.rhs_legacy_vs_edge_src_rel)
                rep.rhs_legacy_vs_edge_src_rel = flip;
        }

        // Restore production RHS on particles for any later use
        setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    }

    // --- Particle div of edge-reconstructed J (≠ edge Kirchhoff r_edge) ---
    std::vector<Real> div_j_im_host(n, Real(0));
    std::vector<Real> div_j_re_host(n, Real(0));
    {
        UpdateCellLinkedList<ExecutionPolicy, RealBody> upd_cll(glass_body);
        UpdateRelation<ExecutionPolicy, Inner<>> upd_rel(inner);
        upd_cll.exec();
        upd_rel.exec();
        execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.j_edge_recon_imag, names.div_j_imag);
        syncVariableToHost<Real>(particles, names.div_j_imag);
        {
            const Real *div_j = particles.getVariableDataByName<Real>(names.div_j_imag);
            for (size_t i = 0; i < n; ++i)
                div_j_im_host[i] = div_j[i];
        }
        // Reuse DivJImag buffer for real-chain divergence (host copy already saved).
        execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.j_edge_recon_real, names.div_j_imag);
        syncVariableToHost<Real>(particles, names.div_j_imag);
        {
            const Real *div_j = particles.getVariableDataByName<Real>(names.div_j_imag);
            for (size_t i = 0; i < n; ++i)
                div_j_re_host[i] = div_j[i];
        }
        // Note: DivJImag particle field currently holds real-chain div after the second
        // exec above; CSV/metrics use div_j_{re,im}_host. Do not treat DivJImag as imag-only
        // after this audit without an explicit restore.

        OphelieDoc199RegionRes div_all_legacy;
        double vol_tot = 0.0, divc2 = 0.0, j2 = 0.0;
        double vol_int = 0.0, divc2_int = 0.0, vol_rim = 0.0, divc2_rim = 0.0;
        size_t n_fb = 0, n_fb_rim = 0, n_rim_part = 0;
        syncVariableToHost<Real>(particles, names.edge_recon_fallback);
        const Real *fb = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
        for (size_t i = 0; i < n; ++i)
        {
            const Real scale = std::max(sigma[i] * omega * a[i].norm(), Real(TinyReal));
            ophelieDoc199AccRegion(div_all_legacy, vol[i], div_j_im_host[i], scale);
            const double v = static_cast<double>(vol[i]);
            vol_tot += v;
            const double dri = static_cast<double>(div_j_re_host[i]);
            const double dii = static_cast<double>(div_j_im_host[i]);
            const double d2 = dri * dri + dii * dii;
            divc2 += v * d2;
            j2 += v * static_cast<double>(j_re[i].squaredNorm() + j_im[i].squaredNorm());
            const char *tag = ophelieDoc199RegionTag(pos[i], center, radius, half_height);
            const bool is_rim = (std::strcmp(tag, "rim") == 0);
            if (is_rim)
            {
                vol_rim += v;
                divc2_rim += v * d2;
                ++n_rim_part;
                if (fb[i] > Real(0.5))
                    ++n_fb_rim;
            }
            else if (std::strcmp(tag, "interior") == 0)
            {
                vol_int += v;
                divc2_int += v * d2;
            }
            if (fb[i] > Real(0.5))
                ++n_fb;
        }
        rep.div_j_recon_abs_l2 = ophelieDoc199AbsL2(div_all_legacy);
        rep.div_j_recon_rel_l2 = ophelieDoc199RelL2(div_all_legacy); // legacy unit-bearing
        rep.L_star = radius; // cylinder R (Doc200 contract)
        rep.J_star = (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(j2 / vol_tot)) : Real(0);
        // Doc201 A4 / plan §7: zero-J manufactured uses prescribed scale J★=RMS(σ ω |A|), not TinyReal.
        double jpres2 = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            const double v = static_cast<double>(vol[i]);
            const double jp = static_cast<double>(sigma[i] * omega * a[i].norm());
            jpres2 += v * jp * jp;
        }
        const Real j_prescribed =
            (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(jpres2 / vol_tot)) : Real(0);
        const Real j_recon = rep.J_star;
        const Real jstar = std::max(j_recon, std::max(j_prescribed, Real(TinyReal)));
        if (j_recon < Real(0.01) * std::max(j_prescribed, Real(TinyReal)))
        {
            std::cout << "[ophelie][doc199][e1-cont] Jstar_floor=prescribed_sigma_omega_|A| J_recon="
                      << j_recon << " J_prescribed=" << j_prescribed << " Jstar_used=" << jstar
                      << std::endl;
            rep.J_star = jstar; // report the scale actually used in epsilon_div
        }
        rep.div_j_complex_abs_l2 =
            (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(divc2 / vol_tot)) : Real(0);
        rep.epsilon_div_complex = rep.L_star * rep.div_j_complex_abs_l2 / jstar;
        const Real div_abs_int =
            (vol_int > TinyReal) ? static_cast<Real>(std::sqrt(divc2_int / vol_int)) : Real(0);
        const Real div_abs_rim =
            (vol_rim > TinyReal) ? static_cast<Real>(std::sqrt(divc2_rim / vol_rim)) : Real(0);
        rep.epsilon_div_interior = rep.L_star * div_abs_int / jstar;
        rep.epsilon_div_rim = rep.L_star * div_abs_rim / jstar;
        rep.recon_fallback_fraction = (n > 0) ? static_cast<Real>(n_fb) / static_cast<Real>(n) : Real(0);
        rep.recon_fallback_fraction_rim =
            (n_rim_part > 0) ? static_cast<Real>(n_fb_rim) / static_cast<Real>(n_rim_part) : Real(0);

        // Doc212: same (L*,J*) on edge Kirchhoff Σq (complex) vs SPH div(J_recon).
        {
            double edge_c2 = 0.0;
            double sum_w = 0.0, sum_e = 0.0, sum_d = 0.0, sum_ee = 0.0, sum_dd = 0.0, sum_ed = 0.0;
            for (size_t i = 0; i < n; ++i)
            {
                const double v = static_cast<double>(vol[i]);
                const double eri = static_cast<double>(edge_res_im_host[i]);
                const double err = static_cast<double>(edge_res_re_host[i]);
                edge_c2 += v * (eri * eri + err * err);
                const double ae = std::abs(eri);
                const double ad = std::abs(static_cast<double>(div_j_im_host[i]));
                sum_w += v;
                sum_e += v * ae;
                sum_d += v * ad;
                sum_ee += v * ae * ae;
                sum_dd += v * ad * ad;
                sum_ed += v * ae * ad;
            }
            const Real edge_abs =
                (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(edge_c2 / vol_tot)) : Real(0);
            rep.epsilon_div_edge_kirchhoff = rep.L_star * edge_abs / jstar;
            rep.sph_vs_edge_eps_ratio =
                rep.epsilon_div_edge_kirchhoff > TinyReal
                    ? rep.epsilon_div_complex / rep.epsilon_div_edge_kirchhoff
                    : Real(0);
            if (sum_w > TinyReal)
            {
                const double mean_e = sum_e / sum_w;
                const double mean_d = sum_d / sum_w;
                const double cov = sum_ed / sum_w - mean_e * mean_d;
                const double ve = sum_ee / sum_w - mean_e * mean_e;
                const double vd = sum_dd / sum_w - mean_d * mean_d;
                rep.corr_edge_res_vs_sph_div_im =
                    (ve > TinyReal && vd > TinyReal)
                        ? static_cast<Real>(cov / std::sqrt(ve * vd))
                        : Real(0);
            }
        }

        // --- Doc201 A4: dual Gradφ+A J (same uncorrected SPH family as DivΣA identity) ---
        {
            InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
                inner, names);
            compute_grad_phi.exec();
            // Real Gradφ: swap φ_real into φ_imag slot, run same kernel, stash GradPhiReal, restore.
            syncVariableToHost<Real>(particles, names.phi_imag);
            syncVariableToHost<Real>(particles, names.phi_real);
            syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
            Real *phi_im_w = particles.getVariableDataByName<Real>(names.phi_imag);
            Real *phi_re_w = particles.getVariableDataByName<Real>(names.phi_real);
            Vecd *gphi_im_w = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
            Vecd *gphi_re_w = particles.getVariableDataByName<Vecd>(names.grad_phi_real);
            std::vector<Real> phi_im_save(n);
            std::vector<Vecd> gphi_im_save(n);
            for (size_t i = 0; i < n; ++i)
            {
                phi_im_save[i] = phi_im_w[i];
                gphi_im_save[i] = gphi_im_w[i];
                phi_im_w[i] = phi_re_w[i];
            }
            syncVariableToDevice<Real>(particles, names.phi_imag);
            compute_grad_phi.exec();
            syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
            for (size_t i = 0; i < n; ++i)
            {
                gphi_re_w[i] = gphi_im_w[i];
                phi_im_w[i] = phi_im_save[i];
                gphi_im_w[i] = gphi_im_save[i];
            }
            syncVariableToDevice<Real>(particles, names.phi_imag);
            syncVariableToDevice<Vecd>(particles, names.grad_phi_imag);
            syncVariableToDevice<Vecd>(particles, names.grad_phi_real);

            syncVariableToHost<Vecd>(particles, names.a_src_real);
            syncVariableToHost<Vecd>(particles, names.a_src_imag);
            syncVariableToHost<Real>(particles, names.sigma);
            syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
            syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
            syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
            syncVariableToHost<Vecd>(particles, names.grad_phi_real);
            const Vecd *gphi_im = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
            const Vecd *gphi_re = particles.getVariableDataByName<Vecd>(names.grad_phi_real);
            const Vecd *a_re = particles.getVariableDataByName<Vecd>(names.a_src_real);
            const Vecd *a_im = particles.getVariableDataByName<Vecd>(names.a_src_imag);
            const Real *sig = particles.getVariableDataByName<Real>(names.sigma);
            const Vecd *jrec_im = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
            const Vecd *jrec_re = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
            Vecd *j_scratch_im = particles.getVariableDataByName<Vecd>(names.j_imag);
            Vecd *j_scratch_re = particles.getVariableDataByName<Vecd>(names.j_real);
            double jdiff2 = 0.0, jrec2 = 0.0;
            for (size_t i = 0; i < n; ++i)
            {
                const Vecd e_im = -params.omega() * a_re[i] - gphi_im[i];
                const Vecd e_re = params.omega() * a_im[i] - gphi_re[i];
                j_scratch_im[i] = sig[i] * e_im;
                j_scratch_re[i] = sig[i] * e_re;
                const double v = static_cast<double>(vol[i]);
                const Vecd d_im = j_scratch_im[i] - jrec_im[i];
                const Vecd d_re = j_scratch_re[i] - jrec_re[i];
                jdiff2 += v * static_cast<double>(d_im.squaredNorm() + d_re.squaredNorm());
                jrec2 += v * static_cast<double>(jrec_im[i].squaredNorm() + jrec_re[i].squaredNorm());
            }
            rep.j_dual_vs_recon_rel =
                (jrec2 > TinyReal) ? static_cast<Real>(std::sqrt(jdiff2 / jrec2)) : Real(0);
            syncVariableToDevice<Vecd>(particles, names.j_imag);
            syncVariableToDevice<Vecd>(particles, names.j_real);
            execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.j_imag, names.div_j_imag);
            syncVariableToHost<Real>(particles, names.div_j_imag);
            std::vector<Real> div_dual_im(n, Real(0));
            {
                const Real *div_j = particles.getVariableDataByName<Real>(names.div_j_imag);
                for (size_t i = 0; i < n; ++i)
                    div_dual_im[i] = div_j[i];
            }
            execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.j_real, names.div_j_imag);
            syncVariableToHost<Real>(particles, names.div_j_imag);
            double divc2_dual = 0.0;
            {
                const Real *div_re = particles.getVariableDataByName<Real>(names.div_j_imag);
                for (size_t i = 0; i < n; ++i)
                {
                    const double v = static_cast<double>(vol[i]);
                    const double dri = static_cast<double>(div_re[i]);
                    const double dii = static_cast<double>(div_dual_im[i]);
                    divc2_dual += v * (dri * dri + dii * dii);
                }
            }
            const Real div_dual_abs =
                (vol_tot > TinyReal) ? static_cast<Real>(std::sqrt(divc2_dual / vol_tot)) : Real(0);
            rep.epsilon_div_dual_gradphi = rep.L_star * div_dual_abs / jstar;
            for (size_t i = 0; i < n; ++i)
            {
                j_scratch_im[i] = jrec_im[i];
                j_scratch_re[i] = jrec_re[i];
            }
            syncVariableToDevice<Vecd>(particles, names.j_imag);
            syncVariableToDevice<Vecd>(particles, names.j_real);
        }
    }

    if (!particle_csv_path.empty())
    {
        syncVariableToHost<Real>(particles, names.edge_recon_fallback);
        const Real *fb_csv = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
        std::ofstream csv(particle_csv_path);
        csv << "i,x,y,z,r,region,Vol,phi_imag,phi_lhs,phi_rhs_prod,lin_res,edge_res,edge_src_phi0,"
               "div_j_re,div_j_im,|J|,recon_fallback,"
               "Jx_re,Jy_re,Jz_re,Jx_im,Jy_im,Jz_im\n";
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - center[0];
            const Real dy = pos[i][1] - center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real lin = lhs_prod[i] - rhs_prod[i];
            const Real jn = std::sqrt(j_re[i].squaredNorm() + j_im[i].squaredNorm());
            csv << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << ","
                << ophelieDoc199RegionTag(pos[i], center, radius, half_height) << "," << vol[i] << "," << phi[i] << ","
                << lhs_prod[i] << "," << rhs_prod[i] << "," << lin << "," << edge_res[i] << "," << edge_src_signed[i]
                << "," << div_j_re_host[i] << "," << div_j_im_host[i] << "," << jn << "," << fb_csv[i] << ","
                << j_re[i][0] << "," << j_re[i][1] << "," << j_re[i][2] << "," << j_im[i][0] << "," << j_im[i][1]
                << "," << j_im[i][2] << "\n";
        }
    }

    std::cout << "[ophelie][doc199][e1-cont] rhs=" << rep.phi_rhs_name << " lhs=" << rep.phi_lhs_name
              << " phi_eq_res_vol=" << rep.phi_eq_res_vol << " edge_res_abs_l2=" << rep.edge_res_l2
              << " edge_res_rel_l2=" << rep.edge_res_rel_l2 << " edge_global_int=" << rep.edge_global_integral
              << " antisym_I_rel_rms=" << rep.antisym_I_rel_rms
              << " rhsDiv_vs_edgeSrc_rel=" << rep.rhs_div_vs_edge_src_rel
              << " rhsLegacy_vs_edgeSrc_rel=" << rep.rhs_legacy_vs_edge_src_rel
              << " divJ_recon_abs_l2=" << rep.div_j_recon_abs_l2
              << " divJ_legacy_rel_l2=" << rep.div_j_recon_rel_l2 << " (unit-bearing; not %)"
              << " Lstar=" << rep.L_star << " Jstar=" << rep.J_star
              << " divJ_complex_abs_l2=" << rep.div_j_complex_abs_l2
              << " epsilon_div=" << rep.epsilon_div_complex
              << " epsilon_div_dual_gradphi=" << rep.epsilon_div_dual_gradphi
              << " j_dual_vs_recon_rel=" << rep.j_dual_vs_recon_rel
              << " n_edges=" << rep.n_edges_unique << std::endl;
    std::cout << "[ophelie][doc199][e1-cont] regions edge_rel_l2 interior=" << ophelieDoc199RelL2(rep.interior_edge)
              << " (n=" << rep.interior_edge.n << ") rim=" << ophelieDoc199RelL2(rep.rim_edge)
              << " (n=" << rep.rim_edge.n << ") all=" << ophelieDoc199RelL2(rep.all_edge) << std::endl;
    std::cout << "[ophelie][doc199][e1-cont][doc212] epsilon_div_interior=" << rep.epsilon_div_interior
              << " epsilon_div_rim=" << rep.epsilon_div_rim
              << " recon_fallback_frac=" << rep.recon_fallback_fraction
              << " recon_fallback_frac_rim=" << rep.recon_fallback_fraction_rim << std::endl;
    std::cout << "[ophelie][doc199][e1-cont][doc212] epsilon_div_edge_kirchhoff="
              << rep.epsilon_div_edge_kirchhoff << " epsilon_div_sph_recon=" << rep.epsilon_div_complex
              << " sph_vs_edge_eps_ratio=" << rep.sph_vs_edge_eps_ratio
              << " corr_edge_res_vs_sph_div_im=" << rep.corr_edge_res_vs_sph_div_im
              << " note: edge_kirchhoff uses Σq with same (L*,J*) as SPH epsilon_div; not a gate;"
              << " edge_rel_l2 keeps local EMF scale"
              << std::endl;
    // Soft gate note (contract): edge_rel <= 1% target; epsilon_div <= 1% (Doc200).
    // Interior edge_rel is reported separately — global edge can be rim-diluted.
    const bool edge_gate = rep.edge_res_rel_l2 <= Real(0.01);
    const bool edge_int_gate = ophelieDoc199RelL2(rep.interior_edge) <= Real(0.01);
    const bool lin_small = rep.phi_eq_res_vol <= Real(0.01);
    const bool eps_gate = rep.epsilon_div_complex <= Real(0.01);
    const bool eps_dual_gate = rep.epsilon_div_dual_gradphi <= Real(0.01);
    std::cout << "[ophelie][doc199][e1-cont] gates edge_rel_l2<=1%=" << (edge_gate ? 1 : 0)
              << " edge_rel_interior<=1%=" << (edge_int_gate ? 1 : 0)
              << " phi_eq_res_vol<=1%=" << (lin_small ? 1 : 0) << " epsilon_div<=1%=" << (eps_gate ? 1 : 0)
              << " epsilon_div_dual<=1%=" << (eps_dual_gate ? 1 : 0)
              << " note: EDGE_CONTINUITY_FAIL if lin_ok && edge_fail; epsilon_div is edge-recon SPH;"
              << " epsilon_div_dual is Gradφ+A; epsilon_div_edge_kirchhoff is Σq pairing (not a gate)"
              << std::endl;
    if (!particle_csv_path.empty())
        std::cout << "[ophelie][doc199][e1-cont] particle_csv=" << particle_csv_path << std::endl;
    if (!edge_csv_path.empty())
        std::cout << "[ophelie][doc199][e1-cont] edge_csv=" << edge_csv_path << std::endl;

    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC199_EDGE_CONTINUITY_AUDIT_H
