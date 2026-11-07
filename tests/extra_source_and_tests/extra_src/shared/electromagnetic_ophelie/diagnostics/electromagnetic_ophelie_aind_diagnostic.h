#ifndef ELECTROMAGNETIC_OPHELIE_AIND_DIAGNOSTIC_H
#define ELECTROMAGNETIC_OPHELIE_AIND_DIAGNOSTIC_H

#include "electromagnetic_ophelie_biot_savart.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "electromagnetic_ophelie_phi_gmres.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_postprocess.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_phi_rhs_diagnostics.h"

#include <algorithm>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real hostVecdVolWeightedNorm(BaseParticles &particles, const std::string &field_name, size_t n)
{
    return std::sqrt(hostVolWeightedVecdNormSquared(particles, field_name, n));
}

/** Doc 61: compare stored ACoil/BCoil to a fresh Biot evaluation at the same Position.
 *  Call immediately after applyMultiloopFilamentBiotToGlass and again after edge-flux. */
struct OphelieStoredCoilBiotAudit
{
    size_t n_checked = 0;
    Real mean_abs_B_ratio = 0.0;
    Real mean_abs_Bz_ratio = 0.0;
    Real max_rel_B_err = 0.0;
    Real midplane_mean_abs_B_ratio = 0.0;
    size_t n_midplane = 0;
    Real soft_m = 0.0;
    size_t n_moments = 0;
    Real sample_r = 0.0;
    Real sample_z = 0.0;
    Real sample_stored_Bz = 0.0;
    Real sample_fresh_Bz = 0.0;
    Real sample_stored_Br = 0.0;
    Real sample_fresh_Br = 0.0;
};

/** Doc 61: rebuild moments with the same L0/L1 excitation used to write ACoil/BCoil. */
inline OphelieStoredCoilBiotAudit auditStoredCoilBiotAgainstFresh(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieFrenchReducedCaseParams &french,
    Real mu0, Real softening_length, const char *tag)
{
    OphelieStoredCoilBiotAudit audit;
    audit.soft_m = softening_length;
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_coil_imag);

    StdVec<OphelieCurrentMomentSample> moments;
    buildFrenchCoilExcitationMoments(french.coil, french.coil_excitation, french.coil_section_radial_width_m,
                                     french.coil_section_axial_height_m, french.coil_section_n_radial,
                                     french.coil_section_n_axial, moments, /*quiet=*/true, french.coil_current_imag);
    audit.n_moments = moments.size();

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *b_stored = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    const Vecd *b_stored_imag = particles.getVariableDataByName<Vecd>(names.b_coil_imag);

    Real sum_B = 0.0;
    Real sum_Bz = 0.0;
    Real sum_B_mid = 0.0;
    Real best_target = 1.0e300;
    size_t best_i = 0;
    const Real target_r = Real(0.20);
    const Real target_z = Real(0.11);

    for (size_t i = 0; i < n; ++i)
    {
        Vecd a_r = Vecd::Zero();
        Vecd a_i = Vecd::Zero();
        Vecd b_r = Vecd::Zero();
        Vecd b_i = Vecd::Zero();
        accumulateBiotSavartFromMoments(moments, pos[i], mu0, softening_length, a_r, a_i, b_r, b_i);
        // Doc151: include imag channel so pure-imag I is auditable.
        const Real B_fresh = std::sqrt(b_r.squaredNorm() + b_i.squaredNorm());
        const Real B_store = std::sqrt(b_stored[i].squaredNorm() + b_stored_imag[i].squaredNorm());
        if (B_fresh <= TinyReal)
        {
            continue;
        }
        const Real ratio_B = B_store / B_fresh;
        const Real ratio_Bz =
            std::hypot(b_stored[i][2], b_stored_imag[i][2]) / (std::hypot(b_r[2], b_i[2]) + TinyReal);
        sum_B += ratio_B;
        sum_Bz += ratio_Bz;
        audit.max_rel_B_err = std::max(audit.max_rel_B_err, std::abs(ratio_B - Real(1)));
        ++audit.n_checked;

        const Real r = std::sqrt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]);
        const Real z = pos[i][2];
        if (z >= Real(0.09) && z <= Real(0.12) && r >= Real(0.05))
        {
            sum_B_mid += ratio_B;
            ++audit.n_midplane;
        }
        const Real d2 = (r - target_r) * (r - target_r) + (z - target_z) * (z - target_z);
        if (d2 < best_target)
        {
            best_target = d2;
            best_i = i;
            audit.sample_r = r;
            audit.sample_z = z;
            const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
            audit.sample_stored_Bz = b_stored[i][2];
            audit.sample_fresh_Bz = b_r[2];
            audit.sample_stored_Br = inv_r * (b_stored[i][0] * pos[i][0] + b_stored[i][1] * pos[i][1]);
            audit.sample_fresh_Br = inv_r * (b_r[0] * pos[i][0] + b_r[1] * pos[i][1]);
        }
        (void)best_i;
    }
    if (audit.n_checked > 0)
    {
        audit.mean_abs_B_ratio = sum_B / static_cast<Real>(audit.n_checked);
        audit.mean_abs_Bz_ratio = sum_Bz / static_cast<Real>(audit.n_checked);
    }
    if (audit.n_midplane > 0)
    {
        audit.midplane_mean_abs_B_ratio = sum_B_mid / static_cast<Real>(audit.n_midplane);
    }

    Real zmin = n > 0 ? pos[0][2] : Real(0);
    Real zmax = zmin;
    for (size_t i = 1; i < n; ++i)
    {
        zmin = std::min(zmin, pos[i][2]);
        zmax = std::max(zmax, pos[i][2]);
    }

    std::cout.setf(std::ios::scientific);
    std::cout.precision(6);
    std::cout << "[ophelie][doc61] stored_vs_fresh_biot tag=" << tag << " n_checked=" << audit.n_checked
              << " n_moments=" << audit.n_moments
              << " excitation=" << ophelieFrenchCoilExcitationLevelName(french.coil_excitation)
              << " soft=" << audit.soft_m
              << " mean_|B|_ratio=" << audit.mean_abs_B_ratio << " mean_|Bz|_ratio=" << audit.mean_abs_Bz_ratio
              << " midplane_mean_|B|_ratio=" << audit.midplane_mean_abs_B_ratio
              << " n_mid=" << audit.n_midplane << " max_rel_|B|_err=" << audit.max_rel_B_err
              << " sample_r=" << audit.sample_r << " sample_z=" << audit.sample_z
              << " sample_Bz_stored=" << audit.sample_stored_Bz << " sample_Bz_fresh=" << audit.sample_fresh_Bz
              << " sample_Br_stored=" << audit.sample_stored_Br << " sample_Br_fresh=" << audit.sample_fresh_Br
              << " I_re=" << french.coil.current_per_loop << " I_im=" << french.coil_current_imag
              << " R=" << french.coil.loop_radius << " n_loops=" << french.coil.num_loops
              << " pos_z_world=[" << zmin << "," << zmax << "]" << std::endl;
    return audit;
}

/** level0/post L2 edge residual ratio; returns -1 when level0 is negligible (ratio undefined). */
inline Real ophelieEdgeResidualReductionRatio(Real level0_l2, Real post_l2)
{
    constexpr Real level0_floor = TinyReal * Real(1e3);
    if (level0_l2 <= level0_floor)
    {
        return Real(-1);
    }
    return level0_l2 / (post_l2 + TinyReal);
}

/** sqrt( sum_i Vol_i * (|v_r_i|^2 + |v_i_i|^2) ) for complex vector potential stored as two Vecd fields. */
inline Real hostComplexVecdPairVolWeightedNorm(BaseParticles &particles, const std::string &field_real,
                                               const std::string &field_imag, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_real);
    syncVariableToHost<Vecd>(particles, field_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *v_r = particles.getVariableDataByName<Vecd>(field_real);
    const Vecd *v_i = particles.getVariableDataByName<Vecd>(field_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t k = 0; k < n; ++k)
    {
        sum += vol[k] * (v_r[k].squaredNorm() + v_i[k].squaredNorm());
    }
    return std::sqrt(sum);
}

inline Real hostComplexVecdPairMax(BaseParticles &particles, const std::string &field_real,
                                   const std::string &field_imag, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_real);
    syncVariableToHost<Vecd>(particles, field_imag);
    const Vecd *v_r = particles.getVariableDataByName<Vecd>(field_real);
    const Vecd *v_i = particles.getVariableDataByName<Vecd>(field_imag);
    Real max_norm = 0.0;
    for (size_t k = 0; k < n; ++k)
    {
        max_norm = std::max(max_norm, std::sqrt(v_r[k].squaredNorm() + v_i[k].squaredNorm()));
    }
    return max_norm;
}

inline void hostZeroVecdField(BaseParticles &particles, const std::string &field_name, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    for (size_t i = 0; i < n; ++i)
    {
        values[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, field_name);
}

inline void hostScaleVecdFieldInPlace(BaseParticles &particles, const std::string &field_name, Real scale, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    for (size_t i = 0; i < n; ++i)
    {
        values[i] *= scale;
    }
    syncVariableToDevice<Vecd>(particles, field_name);
}

/** Measure phi_rhs scale factor for edge-flux (does not modify fields). */
template <class ExecutionPolicy>
inline Real measureOphelieEdgeFluxRhsNormalizationScale(SolidBody &glass_body, Inner<> &glass_inner,
                                                        const OphelieGlassFieldNames &names, OphelieParameters &params,
                                                        size_t n, bool use_a_total_for_rhs,
                                                        Real safe_rhs_l2 = Real(1.0e4),
                                                        Real safe_rhs_max_abs = Real(1.2e3))
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const bool saved_use_a_total = params.use_a_total_for_edge_flux_;
    params.use_a_total_for_edge_flux_ = use_a_total_for_rhs;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, nullptr, 0.0, nullptr);
    const OpheliePhiRhsFingerprint rhs_fp = computeOpheliePhiRhsFingerprint(particles, names.phi_rhs_imag, n);
    params.use_a_total_for_edge_flux_ = saved_use_a_total;

    const Real rhs_max_abs = hostRobustPhiRhsMaxAbs(particles, names.phi_rhs_imag, n);
    const Real scale_l2 = computeOphelieEdgeFluxInputScaleFromRhsL2(rhs_fp.vol_weighted_l2, safe_rhs_l2);
    const Real scale_max = computeOphelieEdgeFluxInputScaleFromRhsMax(rhs_max_abs, safe_rhs_max_abs);
    return std::max(scale_l2, scale_max);
}

template <class ExecutionPolicy>
inline void applyOphelieEdgeFluxFieldScale(SolidBody &glass_body, const OphelieGlassFieldNames &names, Real scale)
{
    if (std::abs(scale - 1.0) <= TinyReal)
    {
        return;
    }
    StateDynamics<ExecutionPolicy, ScaleOphelieElectromagneticFieldsCK> scale_fields(glass_body, names, scale,
                                                                                    scale * scale);
    scale_fields.exec();
}

template <class ExecutionPolicy>
inline Real applyOphelieEdgeFluxInputNormalization(SolidBody &glass_body, Inner<> &glass_inner,
                                                 const OphelieGlassFieldNames &names, OphelieParameters &params,
                                                 size_t n, Real safe_rhs_l2 = Real(1.0e4),
                                                 Real safe_rhs_max_abs = Real(1.2e3),
                                                 Real *rhs_l2_pre_norm_out = nullptr, Real *scale_measured_out = nullptr)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const OpheliePhiRhsFingerprint rhs_fp = computeOpheliePhiRhsFingerprint(particles, names.phi_rhs_imag, n);
    if (rhs_l2_pre_norm_out != nullptr)
    {
        *rhs_l2_pre_norm_out = rhs_fp.vol_weighted_l2;
    }
    const Real scale_measured = measureOphelieEdgeFluxRhsNormalizationScale<ExecutionPolicy>(
        glass_body, glass_inner, names, params, n, false, safe_rhs_l2, safe_rhs_max_abs);
    if (scale_measured_out != nullptr)
    {
        *scale_measured_out = scale_measured;
    }
    const OphelieEdgeFluxNormalizationMode mode = params.edge_flux_normalization_mode_;
    if (!ophelieEdgeFluxUsesFieldScaleRestore(mode))
    {
        std::cout << "[ophelie] edge-flux normalization_mode="
                  << ophelieEdgeFluxNormalizationModeName(mode) << " rhs_l2=" << rhs_fp.vol_weighted_l2
                  << " scale_measured=" << scale_measured << " (no field scaling applied)" << std::endl;
        if (mode == OphelieEdgeFluxNormalizationMode::SolverLocal)
        {
            params.edge_flux_solver_local_rhs_scale_ = scale_measured;
            std::cout << "[ophelie] edge-flux solver-local: physical fields at true scale; "
                         "internal rhs_scale=" << scale_measured << " (phi PCG + edge recon)."
                      << std::endl;
        }
        return scale_measured;
    }
    if (std::abs(scale_measured - 1.0) <= TinyReal)
    {
        return 1.0;
    }
    const Real rhs_max_abs = hostRobustPhiRhsMaxAbs(particles, names.phi_rhs_imag, n);
    std::cout << "[ophelie] edge-flux input normalize mode=field-scale-restore rhs_l2=" << rhs_fp.vol_weighted_l2
              << " rhs_max_abs_p99.5=" << rhs_max_abs << " scale=" << scale_measured << std::endl;
    applyOphelieEdgeFluxFieldScale<ExecutionPolicy>(glass_body, names, scale_measured);
    return scale_measured;
}

struct OphelieEdgeFluxRestoreAudit
{
    Real input_scale = 1.0;
    Real input_scale_measured = 1.0;
    Real rhs_l2_raw = 0.0;
    Real rhs_l2_scaled = 0.0;
    Real j_imag_vol_before_restore = 0.0;
    Real j_imag_vol_after_restore = 0.0;
    Real j_max_before_restore = 0.0;
    Real j_max_after_restore = 0.0;
    Real j_restore_linearity = 0.0;
    Real restore_invariance_error_j = 0.0;
    Real restore_invariance_error_b = 0.0;
    Real restore_invariance_error_p = 0.0;
    Real j_edge_j_imag_max_ratio = 0.0;
    Real b_ind_over_b_coil_before_restore = 0.0;
    Real b_ind_over_b_coil_after_restore = 0.0;
    Real p_recon_before_restore = 0.0;
    Real p_recon_after_restore = 0.0;
    Real e_edge_em_vol_mismatch = 0.0;
    Real p_graph_over_recon = 0.0;
};

/** ||e_edge - (-grad_phi - omega*a_coil_real)||_vol / ||e_edge||_vol for imag chain. */
inline Real hostOphelieImagEdgeEmfMismatchVolRatio(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                   const OphelieParameters &params, size_t n)
{
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    const Vecd *e_edge = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real omega = params.omega();
    Real mismatch_sq = 0.0;
    Real edge_sq = 0.0;
    for (size_t k = 0; k < n; ++k)
    {
        const Vecd e_em = -grad_phi[k] - omega * a_coil_real[k];
        const Vecd delta = e_edge[k] - e_em;
        mismatch_sq += vol[k] * delta.squaredNorm();
        edge_sq += vol[k] * e_edge[k].squaredNorm();
    }
    return std::sqrt(mismatch_sq) / (std::sqrt(edge_sq) + TinyReal);
}

inline void hostImagEmfTermVolNorms(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                    const OphelieParameters &params, size_t n, Real &grad_phi_vol_norm,
                                    Real &omega_a_coil_vol_norm, Real &e_imag_vol_norm)
{
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.e_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(names.e_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real omega = params.omega();
    Real grad_sq = 0.0;
    Real omega_a_sq = 0.0;
    Real e_sq = 0.0;
    for (size_t k = 0; k < n; ++k)
    {
        const Vecd omega_a = omega * a_coil_real[k];
        grad_sq += vol[k] * grad_phi[k].squaredNorm();
        omega_a_sq += vol[k] * omega_a.squaredNorm();
        e_sq += vol[k] * e_imag[k].squaredNorm();
    }
    grad_phi_vol_norm = std::sqrt(grad_sq);
    omega_a_coil_vol_norm = std::sqrt(omega_a_sq);
    e_imag_vol_norm = std::sqrt(e_sq);
}

inline void printOphelieEdgeFluxRestoreAudit(const OphelieEdgeFluxRestoreAudit &audit)
{
    std::cout << "[ophelie] edge-flux restore audit: scale_applied=" << audit.input_scale
              << " scale_measured=" << audit.input_scale_measured << " rhs_l2_raw=" << audit.rhs_l2_raw
              << " rhs_l2_scaled=" << audit.rhs_l2_scaled << " J_vol_pre=" << audit.j_imag_vol_before_restore
              << " J_vol_post=" << audit.j_imag_vol_after_restore << " j_max_pre=" << audit.j_max_before_restore
              << " j_max_post=" << audit.j_max_after_restore << " j_linearity=" << audit.j_restore_linearity
              << " inv_err_J=" << audit.restore_invariance_error_j << " inv_err_B=" << audit.restore_invariance_error_b
              << " inv_err_P=" << audit.restore_invariance_error_p
              << " j_edge/j_imag_ratio=" << audit.j_edge_j_imag_max_ratio
              << " Bind/B_pre=" << audit.b_ind_over_b_coil_before_restore
              << " Bind/B_post=" << audit.b_ind_over_b_coil_after_restore << " P_pre=" << audit.p_recon_before_restore
              << " P_post=" << audit.p_recon_after_restore << " e_edge_em_mismatch=" << audit.e_edge_em_vol_mismatch
              << " p_graph/p_recon=" << audit.p_graph_over_recon << std::endl;
}

template <class ExecutionPolicy>
inline void restoreOphelieEdgeFluxInputNormalization(SolidBody &glass_body, const OphelieGlassFieldNames &names,
                                                   const OphelieParameters &params, Real input_scale, size_t n)
{
    if (std::abs(input_scale - 1.0) <= TinyReal)
    {
        return;
    }
    const Real inv = 1.0 / input_scale;
    const Real power_inv = inv * inv;
    BaseParticles &particles = glass_body.getBaseParticles();
    StateDynamics<ExecutionPolicy, ScaleOphelieElectromagneticFieldsCK> scale_fields(glass_body, names, inv, power_inv);
    scale_fields.exec();
    hostScaleVecdFieldInPlace(particles, names.e_edge_recon_imag, inv, n);
    hostScaleVecdFieldInPlace(particles, names.j_edge_recon_imag, inv, n);
    if (params.edge_flux_complex_)
    {
        hostScaleVecdFieldInPlace(particles, names.e_edge_recon_real, inv, n);
        hostScaleVecdFieldInPlace(particles, names.j_edge_recon_real, inv, n);
        hostScaleScalarFieldInPlace(particles, names.phi_real, inv, n);
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_imag, power_inv, n);
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_real, power_inv, n);
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_complex, power_inv, n);
    }
}

struct OphelieAIndOneWayDiagnostic
{
    Real a_coil_vol_norm = 0.0;
    Real a_ind_vol_norm = 0.0;
    Real b_coil_vol_norm = 0.0;
    Real b_ind_vol_norm = 0.0;
    Real a_coil_real_vol_norm = 0.0;
    Real a_coil_imag_vol_norm = 0.0;
    Real a_ind_real_vol_norm = 0.0;
    Real a_ind_imag_vol_norm = 0.0;
    Real a_src_real_vol_norm = 0.0;
    Real a_src_imag_vol_norm = 0.0;
    Real b_coil_real_vol_norm = 0.0;
    Real b_coil_imag_vol_norm = 0.0;
    Real b_ind_real_vol_norm = 0.0;
    Real b_ind_imag_vol_norm = 0.0;
    Real b_src_real_vol_norm = 0.0;
    Real b_src_imag_vol_norm = 0.0;
    Real max_a_coil = 0.0;
    Real max_a_ind = 0.0;
    Real max_a_ind_real = 0.0;
    Real max_a_ind_imag = 0.0;
    Real max_b_coil = 0.0;
    Real max_b_ind = 0.0;
    Real max_b_ind_real = 0.0;
    Real max_b_ind_imag = 0.0;
    Real a_ind_over_a_coil = 0.0;
    Real b_ind_over_b_coil = 0.0;
    Real a_ind_real_over_a_src_real = 0.0;
    Real a_ind_imag_over_a_src_real = 0.0;
    Real b_ind_real_over_b_src_real = 0.0;
    Real b_ind_imag_over_b_src_real = 0.0;
    Real phi_eq_res_vol = 0.0;
    Real joule_power_w = 0.0;
    Real max_j_real = 0.0;
    Real max_j_imag = 0.0;
    Real p_complex_coil_only = 0.0;
    Real p_complex_total_a = 0.0;
    Real phi_eq_res_vol_total = 0.0;
    Real phi_real_solver_rel_residual = 0.0;
    Real max_j_real_after_feedback = 0.0;
    Real max_j_imag_after_feedback = 0.0;
    Real edge_res_red_imag = 0.0;
    Real edge_res_red_real = 0.0;
    bool feedback_resolve_done = false;
    Real edge_flux_input_scale = 1.0;
    Real edge_flux_input_scale_measured = 1.0;
    Real edge_flux_rhs_l2_scaled = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real j_imag_vol_pre_restore = 0.0;
    Real joule_power_recon_w = 0.0;
    Real joule_power_pre_restore_w = 0.0;
    Real edge_flux_rhs_l2_pre_norm = 0.0;
    Real b_ind_over_b_coil_pre_restore = 0.0;
    Real e_edge_em_vol_mismatch = 0.0;
    Real j_restore_linearity = 0.0;
    Real restore_invariance_error_j = 0.0;
    Real restore_invariance_error_b = 0.0;
    Real restore_invariance_error_p = 0.0;
    Real p_graph_over_recon = 0.0;
    Real grad_phi_imag_vol_norm = 0.0;
    Real omega_a_coil_vol_norm = 0.0;
    Real e_imag_vol_norm = 0.0;
};

/**
 * One-way edge-flux + A_ind cycle after A_coil is already on the glass body.
 * Optionally rescales A_coil for float-safe edge-flux when |A_coil| is very large (native TEAM7 STL).
 */
template <class ExecutionPolicy>
inline OphelieAIndOneWayDiagnostic runOphelieEdgeFluxAIndOneWayAfterCoilBiot(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    Real dp, Real edge_flux_safe_max_a_coil = Real(1.0e4))
{
    OphelieAIndOneWayDiagnostic diag;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);

    Real cumulative_input_scale = applyOphelieEdgeFluxInputNormalization<ExecutionPolicy>(
        glass_body, glass_inner, names, params, n, params.edge_flux_safe_rhs_l2_, params.edge_flux_safe_rhs_max_abs_,
        &diag.edge_flux_rhs_l2_pre_norm, &diag.edge_flux_input_scale_measured);

    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a(glass_body, names);
    combine_a.exec();

    std::cout << "[ophelie] aind_one_way: current_form=" << ophelieCurrentFormKindName(params.ophelie_current_form_)
              << std::endl;

    if (params.enable_phi_correction_)
    {
        if (ophelieUseEdgeFluxElectromotiveRhs(params))
        {
            const bool saved_use_a_total = params.use_a_total_for_edge_flux_;
            params.use_a_total_for_edge_flux_ = false;
            OpheliePhiBoundaryGeometryContext boundary_geom;
            const OpheliePhiBoundaryGeometryContext *boundary_geom_ptr =
                prepareOpheliePhiBoundaryGeomForSolve(particles, names, params, boundary_geom, dp);
            OphelieComplexEdgeFluxSolveReport coil_only_report;
            diag.p_complex_coil_only = execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(
                glass_body, glass_inner, names, params, boundary_geom_ptr, dp, &coil_only_report);
            diag.phi_eq_res_vol = coil_only_report.phi_eq_res_vol_imag;
            diag.phi_real_solver_rel_residual = coil_only_report.phi_real_solver_rel_residual;
            std::cout << "[ophelie] aind_one_way: coil_only P_complex=" << diag.p_complex_coil_only
                      << " phi_eq_res_vol=" << diag.phi_eq_res_vol << std::endl;
            params.use_a_total_for_edge_flux_ = saved_use_a_total;
        }
        else
        {
            (void)solvePhiImag<ExecutionPolicy>(glass_body, glass_inner, names, params);
            applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, glass_inner, names, params);
            diag.phi_eq_res_vol = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);
            InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
                glass_inner, names);
            StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, names, params);
            UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
            UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
            update_cell_linked_list.exec();
            update_inner_relation.exec();
            compute_grad_phi.exec();
            compute_ejq_with_phi.exec();
        }
    }
    else
    {
        StateDynamics<ExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> compute_ejq(glass_body, names, params);
        compute_ejq.exec();
    }

    const std::string &j_imag_source = getOphelieAIndJImagFieldName(names, params);
    const std::string &j_real_source = getOphelieAIndJRealFieldName(names, params);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> compute_aind(
        glass_body, names, params, j_real_source, j_imag_source);
    compute_aind.exec();

    if (params.enable_phi_correction_ && ophelieUseEdgeFluxElectromotiveRhs(params) &&
        params.edge_flux_complex_ && params.aind_one_way_feedback_resolve_)
    {
        const OphelieEdgeFluxComponent imag_component = makeOphelieEdgeFluxImagComponent(names, params);
        const OphelieEdgeFluxComponent real_component = makeOphelieEdgeFluxRealComponent(names, params);
        const OphelieEdgeFluxResidualMetrics imag_level0 =
            evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, glass_inner, names, imag_component,
                                                                         params);
        const OphelieEdgeFluxResidualMetrics real_level0 =
            evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, glass_inner, names, real_component,
                                                                         params);

        StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a_total(glass_body, names);
        combine_a_total.exec();
        const bool saved_use_a_total = params.use_a_total_for_edge_flux_;
        params.use_a_total_for_edge_flux_ = true;
        const Real feedback_extra_scale = measureOphelieEdgeFluxRhsNormalizationScale<ExecutionPolicy>(
            glass_body, glass_inner, names, params, n, true, params.edge_flux_safe_rhs_l2_,
            params.edge_flux_safe_rhs_max_abs_);
        if (ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_) &&
            feedback_extra_scale < Real(1.0) - TinyReal)
        {
            std::cout << "[ophelie] edge-flux feedback normalize extra_scale=" << feedback_extra_scale << std::endl;
            applyOphelieEdgeFluxFieldScale<ExecutionPolicy>(glass_body, names, feedback_extra_scale);
            cumulative_input_scale *= feedback_extra_scale;
        }
        else if (!ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_) &&
                 feedback_extra_scale < Real(1.0) - TinyReal)
        {
            std::cout << "[ophelie] edge-flux feedback: skip extra field-scale (normalization_mode="
                      << ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_) << ")"
                      << std::endl;
        }
        const Real max_j_r_before = hostVecdFieldMax(particles, j_real_source, n);
        const Real max_j_i_before = hostVecdFieldMax(particles, j_imag_source, n);
        // Doc 23 P0-03: do not skip real phi/recon based on pre-feedback J_real.
        // With A_imag from one-way Biot, E_real can be nonzero even if J_real was ~0.
        const Real a_src_imag_norm = hostVecdVolWeightedNorm(particles, names.a_src_imag, n);
        if (max_j_r_before <= TinyReal * std::max(max_j_i_before, Real(1)))
        {
            std::cout << "[ophelie] edge-flux feedback: pre-feedback J_real~0 (max_Jr=" << max_j_r_before
                      << " max_Ji=" << max_j_i_before << ") but still solving real chain; A_src_imag_norm="
                      << a_src_imag_norm << std::endl;
        }
        OpheliePhiBoundaryGeometryContext boundary_geom_fb;
        const OpheliePhiBoundaryGeometryContext *boundary_geom_fb_ptr =
            prepareOpheliePhiBoundaryGeomForSolve(particles, names, params, boundary_geom_fb, dp);
        OphelieComplexEdgeFluxSolveReport total_report;
        diag.p_complex_total_a = execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(
            glass_body, glass_inner, names, params, boundary_geom_fb_ptr, dp, &total_report,
            /*skip_real_phi_solve=*/false);
        diag.phi_eq_res_vol_total = total_report.phi_eq_res_vol_imag;
        diag.phi_real_solver_rel_residual = total_report.phi_real_solver_rel_residual;
        params.use_a_total_for_edge_flux_ = saved_use_a_total;
        diag.feedback_resolve_done = true;

        const OphelieEdgeFluxResidualMetrics imag_post =
            evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, glass_inner, names, imag_component,
                                                                         params);
        const OphelieEdgeFluxResidualMetrics real_post =
            evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(glass_body, glass_inner, names, real_component,
                                                                         params);
        diag.edge_res_red_imag = ophelieEdgeResidualReductionRatio(imag_level0.edge_res_l2, imag_post.edge_res_l2);
        diag.edge_res_red_real = ophelieEdgeResidualReductionRatio(real_level0.edge_res_l2, real_post.edge_res_l2);

        const std::string &j_real_after = getOphelieAIndJRealFieldName(names, params);
        const std::string &j_imag_after = getOphelieAIndJImagFieldName(names, params);
        diag.max_j_real_after_feedback = hostVecdFieldMax(particles, j_real_after, n);
        diag.max_j_imag_after_feedback = hostVecdFieldMax(particles, j_imag_after, n);
        std::cout << "[ophelie] aind_one_way feedback: P_complex_total=" << diag.p_complex_total_a
                  << " P_coil_only=" << diag.p_complex_coil_only
                  << " edge_res_red_imag=" << diag.edge_res_red_imag << " edge_res_red_real=" << diag.edge_res_red_real
                  << " max_J_real=" << diag.max_j_real_after_feedback << " max_J_imag=" << diag.max_j_imag_after_feedback
                  << " A_src_imag_norm=" << hostVecdVolWeightedNorm(particles, names.a_src_imag, n) << std::endl;
    }

    OphelieEdgeFluxRestoreAudit restore_audit;
    restore_audit.input_scale = cumulative_input_scale;
    restore_audit.input_scale_measured = diag.edge_flux_input_scale_measured;
    restore_audit.rhs_l2_raw = diag.edge_flux_rhs_l2_pre_norm;
    restore_audit.rhs_l2_scaled =
        cumulative_input_scale > TinyReal ? diag.edge_flux_rhs_l2_pre_norm * cumulative_input_scale
                                          : diag.edge_flux_rhs_l2_pre_norm;
    restore_audit.j_imag_vol_before_restore = hostVecdVolWeightedNorm(particles, j_imag_source, n);
    restore_audit.j_max_before_restore = hostVecdFieldMax(particles, j_imag_source, n);
    restore_audit.b_ind_over_b_coil_before_restore =
        hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n) /
        (hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n) + TinyReal);
    restore_audit.p_recon_before_restore = hostEdgeFluxReconPower(particles, names, n, params);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list_audit(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation_audit(glass_inner);
    update_cell_linked_list_audit.exec();
    update_inner_relation_audit.exec();
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi_audit(
        glass_inner, names);
    compute_grad_phi_audit.exec();
    restore_audit.e_edge_em_vol_mismatch = hostOphelieImagEdgeEmfMismatchVolRatio(particles, names, params, n);

    if (ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_) &&
        std::abs(cumulative_input_scale - 1.0) > TinyReal)
    {
        restoreOphelieEdgeFluxInputNormalization<ExecutionPolicy>(glass_body, names, params, cumulative_input_scale,
                                                                    n);
    }
    else if (!ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_))
    {
        std::cout << "[ophelie] edge-flux restore skipped (normalization_mode="
                  << ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_) << ")" << std::endl;
    }
    diag.edge_flux_input_scale = cumulative_input_scale;
    diag.edge_flux_rhs_l2_scaled = restore_audit.rhs_l2_scaled;

    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> refresh_aind_after_restore(
        glass_body, names, params, j_real_source, j_imag_source);
    refresh_aind_after_restore.exec();
    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_after_restore(glass_body,
                                                                                                          names);
    combine_after_restore.exec();

    syncGlassElectromagneticFieldsToHost(particles, names);
    diag.a_coil_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_coil_real, names.a_coil_imag, n);
    diag.a_ind_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_ind_real, names.a_ind_imag, n);
    diag.b_coil_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n);
    diag.b_ind_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n);
    diag.a_coil_real_vol_norm = hostVecdVolWeightedNorm(particles, names.a_coil_real, n);
    diag.a_coil_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.a_coil_imag, n);
    diag.a_ind_real_vol_norm = hostVecdVolWeightedNorm(particles, names.a_ind_real, n);
    diag.a_ind_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.a_ind_imag, n);
    diag.a_src_real_vol_norm = hostVecdVolWeightedNorm(particles, names.a_src_real, n);
    diag.a_src_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.a_src_imag, n);
    diag.b_coil_real_vol_norm = hostVecdVolWeightedNorm(particles, names.b_coil_real, n);
    diag.b_coil_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.b_coil_imag, n);
    diag.b_ind_real_vol_norm = hostVecdVolWeightedNorm(particles, names.b_ind_real, n);
    diag.b_ind_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.b_ind_imag, n);
    diag.b_src_real_vol_norm = hostVecdVolWeightedNorm(particles, names.b_src_real, n);
    diag.b_src_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.b_src_imag, n);
    diag.max_a_coil = hostComplexVecdPairMax(particles, names.a_coil_real, names.a_coil_imag, n);
    diag.max_a_ind = hostComplexVecdPairMax(particles, names.a_ind_real, names.a_ind_imag, n);
    diag.max_a_ind_real = hostVecdFieldMax(particles, names.a_ind_real, n);
    diag.max_a_ind_imag = hostVecdFieldMax(particles, names.a_ind_imag, n);
    diag.max_b_coil = hostComplexVecdPairMax(particles, names.b_coil_real, names.b_coil_imag, n);
    diag.max_b_ind = hostComplexVecdPairMax(particles, names.b_ind_real, names.b_ind_imag, n);
    diag.max_b_ind_real = hostVecdFieldMax(particles, names.b_ind_real, n);
    diag.max_b_ind_imag = hostVecdFieldMax(particles, names.b_ind_imag, n);
    diag.max_j_imag = hostVecdFieldMax(particles, j_imag_source, n);
    diag.max_j_real = hostVecdFieldMax(particles, j_real_source, n);
    diag.j_imag_vol_norm = hostVecdVolWeightedNorm(particles, j_imag_source, n);
    diag.j_imag_vol_pre_restore = restore_audit.j_imag_vol_before_restore;
    diag.joule_power_recon_w = hostEdgeFluxReconPower(particles, names, n, params);
    diag.joule_power_pre_restore_w = restore_audit.p_recon_before_restore;
    diag.a_ind_over_a_coil = diag.a_ind_vol_norm / (diag.a_coil_vol_norm + TinyReal);
    diag.b_ind_over_b_coil = diag.b_ind_vol_norm / (diag.b_coil_vol_norm + TinyReal);
    diag.b_ind_over_b_coil_pre_restore = restore_audit.b_ind_over_b_coil_before_restore;
    restore_audit.j_imag_vol_after_restore = diag.j_imag_vol_norm;
    restore_audit.j_max_after_restore = hostVecdFieldMax(particles, j_imag_source, n);
    restore_audit.p_recon_after_restore = diag.joule_power_recon_w;
    const Real j_expected_after_restore =
        restore_audit.j_imag_vol_before_restore / (cumulative_input_scale + TinyReal);
    const Real p_expected_after_restore =
        restore_audit.p_recon_before_restore / (cumulative_input_scale * cumulative_input_scale + TinyReal);
    restore_audit.j_restore_linearity =
        restore_audit.j_max_after_restore /
        (restore_audit.j_max_before_restore / (cumulative_input_scale + TinyReal) + TinyReal);
    restore_audit.restore_invariance_error_j =
        std::abs(diag.j_imag_vol_norm - j_expected_after_restore) / (diag.j_imag_vol_norm + TinyReal);
    // Bind/Bcoil is a ratio: uniform field scaling leaves it unchanged pre/post restore.
    restore_audit.restore_invariance_error_b =
        std::abs(diag.b_ind_over_b_coil - restore_audit.b_ind_over_b_coil_before_restore) /
        (diag.b_ind_over_b_coil + TinyReal);
    restore_audit.restore_invariance_error_p =
        std::abs(diag.joule_power_recon_w - p_expected_after_restore) / (diag.joule_power_recon_w + TinyReal);
    const Real j_edge_max = hostVecdFieldMax(particles, names.j_edge_recon_imag, n);
    const Real j_imag_max = hostVecdFieldMax(particles, names.j_imag, n);
    restore_audit.j_edge_j_imag_max_ratio = j_edge_max / (j_imag_max + TinyReal);
    restore_audit.b_ind_over_b_coil_after_restore = diag.b_ind_over_b_coil;
    {
        const OphelieEdgeFluxPowerAuditDetail power_detail =
            hostOphelieEdgeFluxPowerAuditDetail(particles, names, n, params);
        restore_audit.p_graph_over_recon = power_detail.p_graph_over_recon;
        printOphelieEdgeFluxPowerAuditDetail(power_detail);
    }
    printOphelieEdgeFluxRestoreAudit(restore_audit);
    diag.e_edge_em_vol_mismatch = restore_audit.e_edge_em_vol_mismatch;
    diag.j_restore_linearity = restore_audit.j_restore_linearity;
    diag.restore_invariance_error_j = restore_audit.restore_invariance_error_j;
    diag.restore_invariance_error_b = restore_audit.restore_invariance_error_b;
    diag.restore_invariance_error_p = restore_audit.restore_invariance_error_p;
    diag.p_graph_over_recon = restore_audit.p_graph_over_recon;
    hostImagEmfTermVolNorms(particles, names, params, n, diag.grad_phi_imag_vol_norm, diag.omega_a_coil_vol_norm,
                            diag.e_imag_vol_norm);
    diag.a_ind_real_over_a_src_real = diag.a_ind_real_vol_norm / (diag.a_src_real_vol_norm + TinyReal);
    diag.a_ind_imag_over_a_src_real = diag.a_ind_imag_vol_norm / (diag.a_src_real_vol_norm + TinyReal);
    diag.b_ind_real_over_b_src_real = diag.b_ind_real_vol_norm / (diag.b_src_real_vol_norm + TinyReal);
    diag.b_ind_imag_over_b_src_real = diag.b_ind_imag_vol_norm / (diag.b_src_real_vol_norm + TinyReal);
    diag.joule_power_w = params.edge_flux_complex_ && diag.feedback_resolve_done
                             ? diag.p_complex_total_a
                             : (params.edge_flux_complex_ ? diag.p_complex_coil_only
                                                            : hostVolWeightedSum(particles, names.joule_heat, n));
    return diag;
}

/** Coil Biot -> phi solve (ACoil only) -> J -> A_ind/B_ind one-way; no feedback into A_src. */
template <class ExecutionPolicy>
inline OphelieAIndOneWayDiagnostic runFrenchReducedAIndOneWayDiagnostic(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    applyFrenchReducedCoilExcitationBiotToGlass(particles, names, french, params.mu0_, params.softening_length_);
    (void)auditStoredCoilBiotAgainstFresh(particles, names, french, params.mu0_, params.softening_length_,
                                          "post_biot_pre_edgeflux");
    OphelieAIndOneWayDiagnostic diag =
        runOphelieEdgeFluxAIndOneWayAfterCoilBiot<ExecutionPolicy>(glass_body, glass_inner, names, params, french.dp);
    (void)auditStoredCoilBiotAgainstFresh(particles, names, french, params.mu0_, params.softening_length_,
                                          "post_edgeflux_pre_dump");
    return diag;
}

struct OphelieFrenchSelfInductionPicardResult
{
    /** Unrelaxed fixed-point residual ||F(J)-J|| / ||F|| (doc 23 P0-04). */
    Real final_j_rel_change = 0.0;
    Real fixed_point_residual = 0.0;
    /** Under-relaxation step size; not the convergence metric. */
    Real relaxed_step_change = 0.0;
    size_t self_induction_iterations = 0;
    Real phi_solver_rel_residual = 0.0;
    Real phi_eq_res_vol = 0.0;
    bool picard_converged = false;
    Real joule_power_w = 0.0;
    Real a_ind_over_a_coil = 0.0;
    Real b_ind_over_b_coil = 0.0;
};

/** Multiloop A_coil + Picard A_ind feedback (experimental; not literature_passed). */
template <class ExecutionPolicy>
inline OphelieFrenchSelfInductionPicardResult runFrenchReducedSelfInductionPicard(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french)
{
    OphelieFrenchSelfInductionPicardResult result;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    applyFrenchReducedCoilExcitationBiotToGlass(particles, names, french, params.mu0_, params.softening_length_);
    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);

    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a(glass_body, names);
    combine_a.exec();

    result.final_j_rel_change = runOphelieSelfInductionWithPhiSolve<ExecutionPolicy>(
        glass_body, glass_inner, names, params, result.phi_solver_rel_residual, result.self_induction_iterations,
        result.phi_eq_res_vol, result.picard_converged, &result.relaxed_step_change);
    result.fixed_point_residual = result.final_j_rel_change;

    syncGlassElectromagneticFieldsToHost(particles, names);
    const Real a_coil_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_coil_real, names.a_coil_imag, n);
    const Real a_ind_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_ind_real, names.a_ind_imag, n);
    const Real b_coil_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n);
    const Real b_ind_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n);
    result.a_ind_over_a_coil = a_ind_norm / (a_coil_norm + TinyReal);
    result.b_ind_over_b_coil = b_ind_norm / (b_coil_norm + TinyReal);
    if (params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        result.joule_power_w = hostEdgeFluxReconPower(particles, names, n, params);
    }
    else
    {
        result.joule_power_w = hostVolWeightedSum(particles, names.joule_heat, n);
    }
    if (!result.picard_converged)
    {
        std::cout << "[ophelie] WARNING: Picard did not meet fixed-point/phi tolerances "
                     "(fixed_point_res="
                  << result.fixed_point_residual << " relaxed_step=" << result.relaxed_step_change
                  << " phi_eq_res_vol=" << result.phi_eq_res_vol << " iters=" << result.self_induction_iterations
                  << ")" << std::endl;
    }
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_AIND_DIAGNOSTIC_H
