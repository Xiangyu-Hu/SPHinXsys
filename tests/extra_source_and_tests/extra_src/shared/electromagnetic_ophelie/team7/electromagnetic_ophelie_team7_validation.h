#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_VALIDATION_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_VALIDATION_H

#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_biot_savart.h"
#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_racetrack_source.h"
#include "electromagnetic_ophelie_team7_coil_path_source.h"
#include "electromagnetic_ophelie_team7_probe.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class Team7ValidationLevel
{
    CoilOnly,
    OneWay,
    Picard
};

inline const char *team7ValidationLevelName(Team7ValidationLevel level)
{
    switch (level)
    {
    case Team7ValidationLevel::OneWay:
        return "one-way";
    case Team7ValidationLevel::Picard:
        return "picard";
    default:
        return "coil-only";
    }
}

inline Team7ValidationLevel parseTeam7ValidationLevel(const std::string &name)
{
    if (name == "one-way" || name == "oneway" || name == "l2")
    {
        return Team7ValidationLevel::OneWay;
    }
    if (name == "picard" || name == "l3")
    {
        return Team7ValidationLevel::Picard;
    }
    return Team7ValidationLevel::CoilOnly;
}

enum class Team7ValidationMode
{
    Smoke,
    Strict
};

enum class Team7Phase90Convention
{
    PlusImag,
    MinusImag
};

inline const char *team7ValidationModeName(Team7ValidationMode mode)
{
    return mode == Team7ValidationMode::Strict ? "strict" : "smoke";
}

inline Team7ValidationMode parseTeam7ValidationMode(const std::string &name)
{
    if (name == "strict")
    {
        return Team7ValidationMode::Strict;
    }
    return Team7ValidationMode::Smoke;
}

inline const char *team7Phase90ConventionName(Team7Phase90Convention convention)
{
    return convention == Team7Phase90Convention::MinusImag ? "minus-imag" : "plus-imag";
}

inline Team7Phase90Convention parseTeam7Phase90Convention(const std::string &name)
{
    if (name == "minus-imag" || name == "minus_imag" || name == "neg-ind")
    {
        return Team7Phase90Convention::MinusImag;
    }
    return Team7Phase90Convention::PlusImag;
}

inline bool team7Phase90ConventionUsesNegInd(Team7Phase90Convention convention)
{
    return convention == Team7Phase90Convention::MinusImag;
}

/** Layered TEAM7 pass semantics (smoke vs strict validation). */
struct Team7ValidationPassReport
{
    bool smoke_passed = false;
    bool team7_phase0_source_passed = false;
    bool team7_phase90_probe_passed = false;
    bool diagnostic_only = false;
    bool team7_validation_passed = false;
    Real phase90_rms_abs_mT = 0.0;
    Real phase90_rms_abs_threshold_mT = 0.5;
    StdVec<std::string> diagnostic_reasons;
};

inline std::string team7TaggedLevelBasename(Team7ValidationLevel level, const std::string &output_tag)
{
    const std::string level_name = team7ValidationLevelName(level);
    if (output_tag.empty())
    {
        return level_name;
    }
    return level_name + "_" + output_tag;
}

inline std::string team7OutputArtifactPath(const std::string &prefix, Team7ValidationLevel level,
                                           const std::string &output_tag, const std::string &suffix)
{
    return "./output/" + prefix + team7TaggedLevelBasename(level, output_tag) + suffix;
}

inline std::string makeTeam7AutoOutputTag(Team7ValidationLevel level, Real frequency_hz, Real imag_a_sign,
                                          Team7Phase90Convention phase90_convention, Real coil_source_scale,
                                          bool ind_j_post_scale_auto, Real ind_j_post_scale_applied,
                                          OphelieCoilSourceModel coil_source_model = OphelieCoilSourceModel::VolumeRacetrack)
{
    (void)level;
    std::ostringstream tag;
    tag << "f" << static_cast<int>(frequency_hz + Real(0.5)) << "_"
        << (coil_source_model == OphelieCoilSourceModel::FilamentRacetrack ? "fil" : "vol") << "_asign_"
        << (imag_a_sign < Real(0) ? "m1" : "p1") << "_p90_"
        << (team7Phase90ConventionUsesNegInd(phase90_convention) ? "minus" : "plus") << "_src";
    tag << std::fixed << std::setprecision(3) << coil_source_scale;
    if (ind_j_post_scale_auto)
    {
        tag << "_jauto";
    }
    else if (ind_j_post_scale_applied > TinyReal)
    {
        tag << "_jx" << std::setprecision(3) << ind_j_post_scale_applied;
    }
    else
    {
        tag << "_jraw";
    }
    return tag.str();
}

inline void printTeam7DiagnosticWarnings(Real coil_source_scale, Real imag_a_sign, bool ind_j_post_scale_auto,
                                         Real ind_j_post_scale_applied, Team7ValidationLevel level,
                                         OphelieEdgeFluxNormalizationMode normalization_mode, bool aind_feedback,
                                         bool normalization_sweep_pending)
{
    std::cout << "[team7] edge-flux normalization_mode="
              << ophelieEdgeFluxNormalizationModeName(normalization_mode) << std::endl;
    if (normalization_mode == OphelieEdgeFluxNormalizationMode::SolverLocal &&
        level != Team7ValidationLevel::Picard)
    {
        std::cout << "[team7][warning] solver-local normalization recommended for L3 Picard; L2 one-way remains "
                     "diagnostic."
                  << std::endl;
    }
    if (level == Team7ValidationLevel::Picard && ophelieEdgeFluxUsesFieldScaleRestore(normalization_mode))
    {
        std::cout << "[team7][warning] Picard with field-scale-restore is legacy; prefer solver-local for "
                     "physical-scale iteration."
                  << std::endl;
    }
    if (ophelieEdgeFluxUsesFieldScaleRestore(normalization_mode) &&
        (level == Team7ValidationLevel::Picard || aind_feedback))
    {
        std::cout << "[team7][warning] field-scale-restore with Picard/feedback is diagnostic-only until "
                     "normalization audit passes."
                  << std::endl;
    }
    if (normalization_sweep_pending)
    {
        std::cout << "[team7] normalization sweep enabled (diagnostic)." << std::endl;
    }
    if (std::abs(coil_source_scale - Real(1.0)) > TinyReal)
    {
        std::cout << "[team7][warning] coil_source_scale=" << coil_source_scale << " != 1.0.\n"
                  << "[team7][warning] This run uses post-calibrated source amplitude and cannot be used as strict "
                     "TEAM7 validation.\n"
                  << "[team7][warning] Use source_scale=1.0 with physical filament/validated coil source for TEAM7 "
                     "validation."
                  << std::endl;
    }
    if (std::abs(imag_a_sign - Real(1.0)) > TinyReal)
    {
        std::cout << "[team7][warning] imag_a_sign=" << imag_a_sign << " != +1.\n"
                  << "[team7][warning] This changes the solver equation (debug-only).\n"
                  << "[team7][warning] Do not use imag_a_sign != +1 for TEAM7 validation; use "
                     "--team7-phase90-convention= instead."
                  << std::endl;
    }
    if (ind_j_post_scale_auto || ind_j_post_scale_applied > TinyReal)
    {
        std::cout << "[team7][warning] team7-ind-j-post-scale is enabled (diagnostic J/B scaling).\n"
                  << "[team7][warning] Post-scaled probe metrics cannot be used for TEAM7 validation."
                  << std::endl;
    }
    if (level == Team7ValidationLevel::OneWay)
    {
        std::cout << "[team7] L2 one-way is smoke/diagnostic pipeline (not full TEAM7 eddy-current validation)."
                  << std::endl;
    }
}

inline Team7ValidationPassReport evaluateTeam7ValidationPassReport(
    Team7ValidationLevel level, Team7ValidationMode validation_mode, bool reference_ok, bool coil_path_audit_ok,
    bool em_ok, bool picard_ok, const Team7BzCompareMetrics &coil_metrics,
    const Team7BzCompareMetrics &total_phase0_metrics, const Team7BzCompareMetrics &total_phase90_metrics,
    Real coil_source_scale, Real imag_a_sign, bool ind_j_post_scale_auto, Real ind_j_post_scale_applied,
    OphelieEdgeFluxNormalizationMode normalization_mode, bool aind_feedback,
    Real phase90_rms_abs_threshold_mT = Real(0.5))
{
    Team7ValidationPassReport report;
    report.phase90_rms_abs_threshold_mT = phase90_rms_abs_threshold_mT;
    report.phase90_rms_abs_mT = total_phase90_metrics.rms_abs_error_mT;
    report.team7_phase90_probe_passed =
        reference_ok && total_phase90_metrics.n_probes > 0 &&
        total_phase90_metrics.rms_abs_error_mT <= phase90_rms_abs_threshold_mT;

    if (std::abs(coil_source_scale - Real(1.0)) > TinyReal)
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("coil_source_scale!=1.0");
    }
    if (ind_j_post_scale_auto || ind_j_post_scale_applied > TinyReal)
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("j_post_scale");
    }
    if (std::abs(imag_a_sign - Real(1.0)) > TinyReal)
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("imag_a_sign!=+1");
    }
    if (level == Team7ValidationLevel::OneWay)
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("L2_one-way");
    }
    if (normalization_mode == OphelieEdgeFluxNormalizationMode::SolverLocal &&
        level == Team7ValidationLevel::OneWay)
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("L2_one-way_solver_local");
    }
    if (ophelieEdgeFluxUsesFieldScaleRestore(normalization_mode) &&
        (level == Team7ValidationLevel::Picard || aind_feedback))
    {
        report.diagnostic_only = true;
        report.diagnostic_reasons.push_back("field_scale_restore_feedback_picard");
    }

    const bool bz_ok_for_smoke = !reference_ok || (level == Team7ValidationLevel::Picard ? total_phase0_metrics.passed
                                                                                          : coil_metrics.passed);
    report.smoke_passed = coil_path_audit_ok && em_ok && picard_ok && bz_ok_for_smoke;
    report.team7_phase0_source_passed = reference_ok && coil_metrics.passed;

    if (validation_mode == Team7ValidationMode::Strict && std::abs(coil_source_scale - Real(1.0)) > TinyReal)
    {
        report.smoke_passed = false;
    }

    const bool strict_candidate =
        !report.diagnostic_only && level == Team7ValidationLevel::Picard &&
        validation_mode == Team7ValidationMode::Strict;
    report.team7_validation_passed = strict_candidate && report.team7_phase0_source_passed &&
                                   report.team7_phase90_probe_passed && picard_ok && em_ok && coil_path_audit_ok;
    return report;
}

inline void printTeam7ValidationPassReport(const Team7ValidationPassReport &report)
{
    std::cout << "[team7] pass report: smoke_passed=" << (report.smoke_passed ? 1 : 0)
              << " team7_phase0_source_passed=" << (report.team7_phase0_source_passed ? 1 : 0)
              << " team7_phase90_probe_passed=" << (report.team7_phase90_probe_passed ? 1 : 0)
              << " (threshold_rms_abs_mT=" << report.phase90_rms_abs_threshold_mT
              << " actual=" << report.phase90_rms_abs_mT << ", diagnostic-only)"
              << " diagnostic_only=" << (report.diagnostic_only ? 1 : 0)
              << " team7_validation_passed=" << (report.team7_validation_passed ? 1 : 0);
    if (!report.diagnostic_reasons.empty())
    {
        std::cout << " diagnostic_reasons=";
        for (size_t i = 0; i < report.diagnostic_reasons.size(); ++i)
        {
            if (i > 0)
            {
                std::cout << ",";
            }
            std::cout << report.diagnostic_reasons[i];
        }
    }
    std::cout << std::endl;
}

struct Team7L1SourceValidationReport
{
    OphelieCoilSourceModel source_model = OphelieCoilSourceModel::VolumeRacetrack;
    Real coil_source_scale = 1.0;
    Team7BzCompareMetrics metrics;
    Team7BzCompareMetrics coil_x_span_metrics;
    Real peak_ref_x_mm = 0.0;
    Real peak_sim_x_mm = 0.0;
    Real peak_ref_mT = 0.0;
    Real peak_sim_mT = 0.0;
    Real best_fit_scale_peak = 1.0;
};

inline Team7L1SourceValidationReport makeTeam7L1SourceValidationReport(
    OphelieCoilSourceModel source_model, Real coil_source_scale, const StdVec<Team7BzProbePoint> &probes,
    const Team7BzCompareMetrics &metrics, const Team7BzCompareMetrics &coil_span_metrics)
{
    Team7L1SourceValidationReport report;
    report.source_model = source_model;
    report.coil_source_scale = coil_source_scale;
    report.metrics = metrics;
    report.coil_x_span_metrics = coil_span_metrics;
    if (metrics.n_probes > 0 && !probes.empty())
    {
        report.peak_ref_x_mm = probes[metrics.peak_ref_index].x_mm;
        report.peak_sim_x_mm = probes[metrics.peak_sim_index].x_mm;
        report.peak_ref_mT = metrics.peak_ref_mT;
        report.peak_sim_mT = metrics.peak_sim_mT;
        if (std::abs(metrics.peak_sim_mT) > TinyReal)
        {
            report.best_fit_scale_peak = metrics.peak_ref_mT / metrics.peak_sim_mT;
        }
    }
    return report;
}

inline void printTeam7L1SourceValidationReport(const Team7L1SourceValidationReport &report)
{
    std::cout << "[team7] L1 source validation model=" << ophelieCoilSourceModelName(report.source_model)
              << " source_scale=" << report.coil_source_scale << " rms_rel=" << report.metrics.rms_rel_error
              << " rms_abs_mT=" << report.metrics.rms_abs_error_mT
              << " coil_x_span_rms=" << report.coil_x_span_metrics.rms_rel_error << " peak_sim/ref="
              << report.peak_sim_mT << "/" << report.peak_ref_mT << " peak_x_sim/ref=" << report.peak_sim_x_mm << "/"
              << report.peak_ref_x_mm << " best_fit_scale_peak=" << report.best_fit_scale_peak
              << " (diagnostic; not used for validation_passed)" << std::endl;
}

inline bool parseTeam7ValidationLevelCli(int ac, char *av[], Team7ValidationLevel &level)
{
    static constexpr const char *k_prefix = "--team7-level=";
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], k_prefix, std::strlen(k_prefix)) == 0)
        {
            level = parseTeam7ValidationLevel(std::string(av[i] + std::strlen(k_prefix)));
            return true;
        }
    }
    return false;
}

struct Team7ProbeBzDecomposition
{
    Real bz_coil_mT = 0.0;
    Real bz_ind_real_mT = 0.0;
    Real bz_ind_imag_mT = 0.0;
    /** SPH-interpolated plate b_ind_imag at probe (sanity vs Biot(J)). */
    Real bz_ind_imag_particle_mT = 0.0;
    /** Biot(J) using plate-top skin layer only (diagnostic for depth leakage). */
    Real bz_ind_imag_skin_mT = 0.0;
    /** Biot(J_imag) with −1 sign (phasor convention probe). */
    Real bz_ind_imag_neg_mT = 0.0;
    /** Biot(J): same J·dV moment, source point relocated to plate top z (thin-sheet probe). */
    Real bz_ind_imag_sheet_top_mT = 0.0;
    /** Biot(J) from particles in z_mid ± band only. */
    Real bz_ind_imag_zmid_mT = 0.0;
    Real bz_total_real_mT = 0.0;
    Real bz_total_imag_mT = 0.0;
};

struct Team7ProbePhase90BiotVariantMetrics
{
    Team7BzCompareMetrics volume_full;
    Team7BzCompareMetrics volume_neg;
    Team7BzCompareMetrics skin_layer;
    Team7BzCompareMetrics sheet_top_relocate;
    Team7BzCompareMetrics sheet_top_neg;
    Team7BzCompareMetrics zmid_band;
};

struct Team7OneWayEdgeFluxSummary
{
    Real p_complex_coil_only = 0.0;
    Real p_complex_total = 0.0;
    Real phi_eq_res_vol_imag = 0.0;
    Real phi_eq_res_vol_real = 0.0;
    Real a_ind_over_a_coil = 0.0;
    Real b_ind_over_b_coil = 0.0;
    Real max_j_real = 0.0;
    Real max_j_imag = 0.0;
    Real edge_flux_input_scale = 1.0;
    Real edge_flux_input_scale_measured = 1.0;
    Real edge_flux_rhs_l2_scaled = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real j_imag_vol_pre_restore = 0.0;
    Real joule_power_recon_w = 0.0;
    Real joule_power_pre_restore_w = 0.0;
    Real b_ind_over_b_coil_pre_restore = 0.0;
    Real edge_flux_rhs_l2_pre_norm = 0.0;
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

struct Team7NormalizationSweepRecord
{
    std::string normalization_mode;
    Real safe_rhs_l2 = 0.0;
    Real safe_rhs_max_abs = 0.0;
    Real input_scale_applied = 1.0;
    Real input_scale_measured = 1.0;
    Real rhs_l2_raw = 0.0;
    Real rhs_l2_scaled = 0.0;
    Real j_imag_vol_pre_restore = 0.0;
    Real j_imag_vol_post_restore = 0.0;
    Real b_ind_over_bcoil_pre_restore = 0.0;
    Real b_ind_over_bcoil_post_restore = 0.0;
    Real p_recon_pre_restore = 0.0;
    Real p_recon_post_restore = 0.0;
    Real phase90_rms_post = 0.0;
    Real restore_invariance_error_j = 0.0;
    Real restore_invariance_error_b = 0.0;
    Real restore_invariance_error_p = 0.0;
    Real j_restore_linearity = 0.0;
    Real bind_over_bcoil_post = 0.0;
};

struct Team7EdgeFluxOmegaScalingRecord
{
    Real frequency_hz = 0.0;
    Real omega_rad_s = 0.0;
    Real skin_depth_mm = 0.0;
    Real edge_flux_input_scale = 1.0;
    Real rhs_l2_pre_norm = 0.0;
    Real grad_phi_imag_vol_norm = 0.0;
    Real omega_a_coil_vol_norm = 0.0;
    Real e_imag_vol_norm = 0.0;
    Real omega_a_over_grad_phi = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real j_imag_over_omega = 0.0;
    Real bind_over_bcoil = 0.0;
    Real p_recon_w = 0.0;
    Real e_edge_em_mismatch = 0.0;
    Real j_restore_linearity = 0.0;
    Real p_graph_over_recon = 0.0;
};

inline Real team7AnalyticSkinDepthMm(Real frequency_hz, Real sigma, Real mu0 = Real(4e-7) * Pi)
{
    if (frequency_hz <= TinyReal || sigma <= TinyReal)
    {
        return 0.0;
    }
    const Real omega = Real(2) * Pi * frequency_hz;
    const Real delta_m = std::sqrt(Real(2) / (omega * mu0 * sigma));
    return delta_m * 1000.0;
}

inline Team7EdgeFluxOmegaScalingRecord makeTeam7EdgeFluxOmegaScalingRecord(Real frequency_hz, Real sigma,
                                                                         const Team7OneWayEdgeFluxSummary &summary)
{
    Team7EdgeFluxOmegaScalingRecord record;
    record.frequency_hz = frequency_hz;
    record.omega_rad_s = Real(2) * Pi * frequency_hz;
    record.skin_depth_mm = team7AnalyticSkinDepthMm(frequency_hz, sigma);
    record.edge_flux_input_scale = summary.edge_flux_input_scale;
    record.rhs_l2_pre_norm = summary.edge_flux_rhs_l2_pre_norm;
    record.grad_phi_imag_vol_norm = summary.grad_phi_imag_vol_norm;
    record.omega_a_coil_vol_norm = summary.omega_a_coil_vol_norm;
    record.e_imag_vol_norm = summary.e_imag_vol_norm;
    record.omega_a_over_grad_phi =
        record.omega_a_coil_vol_norm / (record.grad_phi_imag_vol_norm + TinyReal);
    record.j_imag_vol_norm = summary.j_imag_vol_norm;
    record.j_imag_over_omega = record.j_imag_vol_norm / (record.omega_rad_s + TinyReal);
    record.bind_over_bcoil = summary.b_ind_over_b_coil;
    record.p_recon_w = summary.joule_power_recon_w;
    record.e_edge_em_mismatch = summary.e_edge_em_vol_mismatch;
    record.j_restore_linearity = summary.j_restore_linearity;
    record.p_graph_over_recon = summary.p_graph_over_recon;
    return record;
}

inline void writeTeam7EdgeFluxOmegaScalingCsv(const std::string &output_path,
                                              const StdVec<Team7EdgeFluxOmegaScalingRecord> &records)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 omega scaling CSV: " << output_path << std::endl;
        return;
    }
    out << "frequency_hz,omega_rad_s,skin_depth_mm,edge_flux_input_scale,rhs_l2_pre_norm,grad_phi_imag_vol,"
           "omega_a_coil_vol,e_imag_vol,omega_a_over_grad_phi,j_imag_vol,j_imag_over_omega,Bind_over_Bcoil,"
           "P_recon_W,e_edge_em_mismatch,j_restore_linearity,p_graph_over_p_recon\n";
    out << std::setprecision(10);
    for (const Team7EdgeFluxOmegaScalingRecord &record : records)
    {
        out << record.frequency_hz << "," << record.omega_rad_s << "," << record.skin_depth_mm << ","
            << record.edge_flux_input_scale << "," << record.rhs_l2_pre_norm << "," << record.grad_phi_imag_vol_norm
            << "," << record.omega_a_coil_vol_norm << "," << record.e_imag_vol_norm << "," << record.omega_a_over_grad_phi
            << "," << record.j_imag_vol_norm << "," << record.j_imag_over_omega << "," << record.bind_over_bcoil << ","
            << record.p_recon_w << "," << record.e_edge_em_mismatch << "," << record.j_restore_linearity << ","
            << record.p_graph_over_recon << "\n";
    }
    std::cout << "[ophelie] TEAM7 edge-flux omega scaling CSV: " << output_path << std::endl;
}

inline void printTeam7EdgeFluxOmegaScalingReport(const StdVec<Team7EdgeFluxOmegaScalingRecord> &records,
                                                 Real reference_frequency_hz = Real(50))
{
    const Team7EdgeFluxOmegaScalingRecord *ref = nullptr;
    for (const Team7EdgeFluxOmegaScalingRecord &record : records)
    {
        if (std::abs(record.frequency_hz - reference_frequency_hz) < Real(0.5))
        {
            ref = &record;
            break;
        }
    }
    std::cout << "[team7] edge-flux omega scaling audit (ref_f=" << reference_frequency_hz << " Hz):" << std::endl;
    for (const Team7EdgeFluxOmegaScalingRecord &record : records)
    {
        Real bind_ratio = 0.0;
        Real j_ratio = 0.0;
        Real f_ratio = 0.0;
        if (ref != nullptr && ref->bind_over_bcoil > TinyReal)
        {
            bind_ratio = record.bind_over_bcoil / ref->bind_over_bcoil;
            j_ratio = record.j_imag_vol_norm / (ref->j_imag_vol_norm + TinyReal);
            f_ratio = record.frequency_hz / (ref->frequency_hz + TinyReal);
        }
        std::cout << "[team7]   f=" << record.frequency_hz << " Hz delta_mm=" << record.skin_depth_mm
                  << " rhs_l2=" << record.rhs_l2_pre_norm << " input_scale=" << record.edge_flux_input_scale
                  << " omegaA/gradPhi=" << record.omega_a_over_grad_phi << " J_imag_vol=" << record.j_imag_vol_norm
                  << " Bind/B=" << record.bind_over_bcoil << " P_recon=" << record.p_recon_w
                  << " e_edge_em_mis=" << record.e_edge_em_mismatch;
        if (ref != nullptr && record.frequency_hz != ref->frequency_hz)
        {
            std::cout << " | Bind/B_ratio=" << bind_ratio << " J_ratio=" << j_ratio << " f_ratio=" << f_ratio
                      << " (Bind/B vs f_ratio=" << (bind_ratio / (f_ratio + TinyReal))
                      << ", vs f^2=" << (bind_ratio / (f_ratio * f_ratio + TinyReal)) << ")";
        }
        std::cout << std::endl;
    }
}

inline Real team7ProbeWendlandC4Weight(Real distance, Real kernel_h)
{
    const Real q = distance / (Real(2) * kernel_h + TinyReal);
    if (q >= Real(1))
    {
        return 0.0;
    }
    const Real t = Real(1) - q;
    return t * t * t * t * (Real(1) + Real(4) * q);
}

/** SPH kernel interpolation of plate Vecd field Bz at air probes (diagnostic). */
inline void evaluatePlateBzAtProbesSphInterp(BaseParticles &plate_particles, const std::string &field_name,
                                             Real kernel_h, const StdVec<Team7BzProbePoint> &probes,
                                             StdVec<Real> &bz_mT)
{
    bz_mT.assign(probes.size(), 0.0);
    syncVariableToHost<Vecd>(plate_particles, "Position");
    syncVariableToHost<Vecd>(plate_particles, field_name);
    const Vecd *pos = plate_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *field = plate_particles.getVariableDataByName<Vecd>(field_name);
    const size_t n_plate = plate_particles.TotalRealParticles();
    for (size_t ip = 0; ip < probes.size(); ++ip)
    {
        Real weight_sum = 0.0;
        Real bz_sum = 0.0;
        for (size_t i = 0; i < n_plate; ++i)
        {
            const Real weight = team7ProbeWendlandC4Weight((probes[ip].position_m - pos[i]).norm(), kernel_h);
            if (weight <= TinyReal)
            {
                continue;
            }
            weight_sum += weight;
            bz_sum += weight * field[i][2];
        }
        bz_mT[ip] = (weight_sum > TinyReal ? bz_sum / weight_sum : Real(0)) * 1000.0;
    }
}

/** Host Biot–Savart Bz from volume current density on a source body (probe line in air). */
inline void evaluateVolumeCurrentBiotSavartBzAtProbes(
    BaseParticles &source_particles, const std::string &j_real_field, const std::string &j_imag_field, Real mu0,
    Real softening_length, const StdVec<Team7BzProbePoint> &probes, StdVec<Real> &bz_real_mT, StdVec<Real> &bz_imag_mT,
    const std::string &skip_particle_mask_field = std::string(),
    Real source_z_min_m = -std::numeric_limits<Real>::infinity(),
    Real source_z_max_m = std::numeric_limits<Real>::infinity(),
    Real relocate_source_z_m = std::numeric_limits<Real>::quiet_NaN(), Real ind_imag_sign = Real(1))
{
    bz_real_mT.assign(probes.size(), 0.0);
    bz_imag_mT.assign(probes.size(), 0.0);
    syncVariableToHost<Vecd>(source_particles, "Position");
    syncVariableToHost<Vecd>(source_particles, j_real_field);
    syncVariableToHost<Vecd>(source_particles, j_imag_field);
    syncVariableToHost<Real>(source_particles, "VolumetricMeasure");
    const bool skip_masked = !skip_particle_mask_field.empty();
    const Real *skip_mask = nullptr;
    if (skip_masked)
    {
        syncVariableToHost<Real>(source_particles, skip_particle_mask_field);
        skip_mask = source_particles.getVariableDataByName<Real>(skip_particle_mask_field);
    }
    const Vecd *pos = source_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_real = source_particles.getVariableDataByName<Vecd>(j_real_field);
    const Vecd *j_imag = source_particles.getVariableDataByName<Vecd>(j_imag_field);
    const Real *vol = source_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n_src = source_particles.TotalRealParticles();
    const Real coeff = mu0 / (4.0 * Pi);
    const Real eps2 = softening_length * softening_length;
    const bool relocate_z = std::isfinite(relocate_source_z_m);
    for (size_t ip = 0; ip < probes.size(); ++ip)
    {
        Vecd b_real = Vecd::Zero();
        Vecd b_imag = Vecd::Zero();
        for (size_t j = 0; j < n_src; ++j)
        {
            if (skip_masked && skip_mask[j] > Real(0.5))
            {
                continue;
            }
            if (pos[j][2] < source_z_min_m || pos[j][2] > source_z_max_m)
            {
                continue;
            }
            const Vecd source_pos =
                relocate_z ? Vecd(pos[j][0], pos[j][1], relocate_source_z_m) : pos[j];
            const Vecd r = probes[ip].position_m - source_pos;
            const Real r2 = r.squaredNorm() + eps2;
            const Real inv_r3 = 1.0 / (std::sqrt(r2) * r2);
            const Vecd jvr = j_real[j] * vol[j];
            const Vecd jvi = j_imag[j] * vol[j];
            b_real += coeff * jvr.cross(r) * inv_r3;
            b_imag += coeff * jvi.cross(r) * inv_r3;
        }
        bz_real_mT[ip] = b_real[2] * 1000.0;
        bz_imag_mT[ip] = ind_imag_sign * b_imag[2] * 1000.0;
    }
}

inline StdVec<Team7ProbeBzDecomposition> evaluateTeam7ProbeBzDecomposition(
    BaseParticles &coil_particles, const OphelieCoilFieldNames &coil_names, BaseParticles &plate_particles,
    const OphelieGlassFieldNames &plate_names, const OphelieParameters &params,
    const std::string &plate_j_real_field, const std::string &plate_j_imag_field,
    const StdVec<Team7BzProbePoint> &probes, Real plate_z_skin_min_m = -std::numeric_limits<Real>::infinity(),
    Real plate_z_top_m = std::numeric_limits<Real>::quiet_NaN(),
    Real plate_z_mid_m = std::numeric_limits<Real>::quiet_NaN(), Real plate_z_band_half_m = Real(0),
    const StdVec<OphelieCurrentMomentSample> *filament_moments = nullptr)
{
    StdVec<Team7BzProbePoint> coil_probes;
    if (filament_moments != nullptr && !filament_moments->empty())
    {
        evaluateTeam7FilamentBiotSavartBzAtProbes(*filament_moments, params.mu0_, params.softening_length_, probes,
                                                  coil_probes);
    }
    else
    {
        evaluateCoilBiotSavartBzAtProbes(coil_particles, coil_names, params, probes, coil_probes);
    }
    StdVec<Real> bz_ind_real;
    StdVec<Real> bz_ind_imag;
    evaluateVolumeCurrentBiotSavartBzAtProbes(plate_particles, plate_j_real_field, plate_j_imag_field, params.mu0_,
                                              params.softening_length_, probes, bz_ind_real, bz_ind_imag,
                                              plate_names.edge_recon_fallback);
    StdVec<Real> bz_ind_imag_skin;
    if (std::isfinite(plate_z_skin_min_m))
    {
        StdVec<Real> bz_ind_real_skin_unused;
        evaluateVolumeCurrentBiotSavartBzAtProbes(
            plate_particles, plate_j_real_field, plate_j_imag_field, params.mu0_, params.softening_length_, probes,
            bz_ind_real_skin_unused, bz_ind_imag_skin, plate_names.edge_recon_fallback, plate_z_skin_min_m);
    }
    StdVec<Real> bz_ind_imag_particle;
    const Real interp_kernel_h = params.softening_length_ * Real(4);
    evaluatePlateBzAtProbesSphInterp(plate_particles, plate_names.b_ind_imag, interp_kernel_h, probes,
                                     bz_ind_imag_particle);
    StdVec<Real> bz_ind_imag_sheet_top;
    StdVec<Real> bz_ind_imag_zmid;
    StdVec<Real> bz_ind_real_unused;
    if (std::isfinite(plate_z_top_m))
    {
        evaluateVolumeCurrentBiotSavartBzAtProbes(plate_particles, plate_j_real_field, plate_j_imag_field, params.mu0_,
                                                 params.softening_length_, probes, bz_ind_real_unused,
                                                 bz_ind_imag_sheet_top, plate_names.edge_recon_fallback,
                                                 -std::numeric_limits<Real>::infinity(),
                                                 std::numeric_limits<Real>::infinity(), plate_z_top_m);
    }
    if (std::isfinite(plate_z_mid_m) && plate_z_band_half_m > TinyReal)
    {
        evaluateVolumeCurrentBiotSavartBzAtProbes(
            plate_particles, plate_j_real_field, plate_j_imag_field, params.mu0_, params.softening_length_, probes,
            bz_ind_real_unused, bz_ind_imag_zmid, plate_names.edge_recon_fallback, plate_z_mid_m - plate_z_band_half_m,
            plate_z_mid_m + plate_z_band_half_m);
    }
    StdVec<Team7ProbeBzDecomposition> out;
    out.resize(probes.size());
    for (size_t i = 0; i < probes.size(); ++i)
    {
        out[i].bz_coil_mT = coil_probes[i].bz_sim_mT;
        out[i].bz_ind_real_mT = bz_ind_real[i];
        out[i].bz_ind_imag_mT = bz_ind_imag[i];
        out[i].bz_ind_imag_skin_mT =
            std::isfinite(plate_z_skin_min_m) && i < bz_ind_imag_skin.size() ? bz_ind_imag_skin[i] : Real(0);
        out[i].bz_ind_imag_particle_mT = bz_ind_imag_particle[i];
        out[i].bz_ind_imag_neg_mT = -out[i].bz_ind_imag_mT;
        out[i].bz_ind_imag_sheet_top_mT =
            std::isfinite(plate_z_top_m) && i < bz_ind_imag_sheet_top.size() ? bz_ind_imag_sheet_top[i] : Real(0);
        out[i].bz_ind_imag_zmid_mT =
            std::isfinite(plate_z_mid_m) && i < bz_ind_imag_zmid.size() ? bz_ind_imag_zmid[i] : Real(0);
        out[i].bz_total_real_mT = out[i].bz_coil_mT + out[i].bz_ind_real_mT;
        out[i].bz_total_imag_mT = out[i].bz_ind_imag_mT;
    }
    return out;
}

inline Team7ProbePhase90BiotVariantMetrics evaluateTeam7ProbePhase90BiotVariants(
    const StdVec<Team7BzProbePoint> &reference_probes, const StdVec<Team7ProbeBzDecomposition> &decomp,
    Real smoke_rms_threshold = Real(1.0))
{
    Team7ProbePhase90BiotVariantMetrics metrics;
    StdVec<Real> phase90_ref_mT;
    phase90_ref_mT.resize(reference_probes.size());
    for (size_t i = 0; i < reference_probes.size(); ++i)
    {
        phase90_ref_mT[i] = reference_probes[i].bz_ref_phase90_mT;
    }
    auto compare_variant = [&](const auto &getter) {
        StdVec<Real> sim_mT;
        sim_mT.resize(decomp.size());
        for (size_t i = 0; i < decomp.size(); ++i)
        {
            sim_mT[i] = getter(decomp[i]);
        }
        return compareTeam7BzAgainstReference(reference_probes, sim_mT, phase90_ref_mT, smoke_rms_threshold);
    };
    metrics.volume_full = compare_variant([](const Team7ProbeBzDecomposition &d) { return d.bz_ind_imag_mT; });
    metrics.volume_neg = compare_variant([](const Team7ProbeBzDecomposition &d) { return d.bz_ind_imag_neg_mT; });
    metrics.skin_layer = compare_variant([](const Team7ProbeBzDecomposition &d) { return d.bz_ind_imag_skin_mT; });
    metrics.sheet_top_relocate =
        compare_variant([](const Team7ProbeBzDecomposition &d) { return d.bz_ind_imag_sheet_top_mT; });
    metrics.sheet_top_neg =
        compare_variant([](const Team7ProbeBzDecomposition &d) { return -d.bz_ind_imag_sheet_top_mT; });
    metrics.zmid_band = compare_variant([](const Team7ProbeBzDecomposition &d) { return d.bz_ind_imag_zmid_mT; });
    return metrics;
}

inline void printTeam7ProbePhase90BiotVariantReport(const Team7ProbePhase90BiotVariantMetrics &metrics)
{
    std::cout << "[ophelie] TEAM7 phase90 Biot variants RMS (vs ref phase90, lower=better):" << std::endl;
    std::cout << "[ophelie]   volume_full=" << metrics.volume_full.rms_rel_error
              << " volume_neg=" << metrics.volume_neg.rms_rel_error << " skin=" << metrics.skin_layer.rms_rel_error
              << " sheet_top=" << metrics.sheet_top_relocate.rms_rel_error
              << " sheet_top_neg=" << metrics.sheet_top_neg.rms_rel_error
              << " zmid_band=" << metrics.zmid_band.rms_rel_error << std::endl;
    std::cout << "[ophelie]   peak_sim full/neg/sheet_top_neg=" << metrics.volume_full.peak_sim_mT << "/"
              << metrics.volume_neg.peak_sim_mT << "/" << metrics.sheet_top_neg.peak_sim_mT
              << " ref_peak=" << metrics.volume_full.peak_ref_mT << std::endl;
}

inline Real hostOphelieEdgeFluxFallbackFraction(BaseParticles &particles, const OphelieGlassFieldNames &names)
{
    syncVariableToHost<Real>(particles, names.edge_recon_fallback);
    const Real *fallback = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        return 0.0;
    }
    size_t count = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (fallback[i] > Real(0.5))
        {
            ++count;
        }
    }
    return static_cast<Real>(count) / static_cast<Real>(n);
}

inline void team7ProbeComplexMagnitudeSimAndRef(const StdVec<Team7BzProbePoint> &reference_probes,
                                                const StdVec<Team7ProbeBzDecomposition> &decomp, StdVec<Real> &sim_mT,
                                                StdVec<Real> &ref_mT)
{
    sim_mT.resize(decomp.size());
    ref_mT.resize(reference_probes.size());
    for (size_t i = 0; i < decomp.size(); ++i)
    {
        sim_mT[i] = std::sqrt(decomp[i].bz_total_real_mT * decomp[i].bz_total_real_mT +
                               decomp[i].bz_total_imag_mT * decomp[i].bz_total_imag_mT);
        ref_mT[i] = std::sqrt(reference_probes[i].bz_ref_phase0_mT * reference_probes[i].bz_ref_phase0_mT +
                              reference_probes[i].bz_ref_phase90_mT * reference_probes[i].bz_ref_phase90_mT);
    }
}

inline StdVec<Team7BzProbePoint> team7ProbesFromDecompositionReal(const StdVec<Team7BzProbePoint> &reference_probes,
                                                                  const StdVec<Team7ProbeBzDecomposition> &decomp,
                                                                  bool use_total)
{
    StdVec<Team7BzProbePoint> probes = reference_probes;
    for (size_t i = 0; i < probes.size(); ++i)
    {
        probes[i].bz_sim_mT = use_total ? decomp[i].bz_total_real_mT : decomp[i].bz_coil_mT;
    }
    return probes;
}

inline void team7ProbePhase90SimAndRef(const StdVec<Team7BzProbePoint> &reference_probes,
                                       const StdVec<Team7ProbeBzDecomposition> &decomp, StdVec<Real> &sim_mT,
                                       StdVec<Real> &ref_mT, bool use_neg_ind_imag = false)
{
    sim_mT.resize(decomp.size());
    ref_mT.resize(reference_probes.size());
    for (size_t i = 0; i < decomp.size(); ++i)
    {
        sim_mT[i] = use_neg_ind_imag ? decomp[i].bz_ind_imag_neg_mT : decomp[i].bz_total_imag_mT;
        ref_mT[i] = reference_probes[i].bz_ref_phase90_mT;
    }
}

inline void team7ProbePhase90SkinSimAndRef(const StdVec<Team7BzProbePoint> &reference_probes,
                                           const StdVec<Team7ProbeBzDecomposition> &decomp, StdVec<Real> &sim_mT,
                                           StdVec<Real> &ref_mT)
{
    sim_mT.resize(decomp.size());
    ref_mT.resize(reference_probes.size());
    for (size_t i = 0; i < decomp.size(); ++i)
    {
        sim_mT[i] = decomp[i].bz_ind_imag_skin_mT;
        ref_mT[i] = reference_probes[i].bz_ref_phase90_mT;
    }
}

struct Team7ProbeXBandSpec
{
    std::string name;
    Real x_min_mm = 0.0;
    Real x_max_mm = 0.0;
};

struct Team7ProbeXBandMetrics
{
    Team7ProbeXBandSpec band;
    size_t n_probes = 0;
    Team7BzCompareMetrics phase0_coil;
    Team7BzCompareMetrics phase90_ind_full;
    Team7BzCompareMetrics phase90_ind_skin;
};

inline StdVec<Team7ProbeXBandSpec> defaultTeam7ProbePhase90XBands(Real coil_x_min_mm, Real coil_x_max_mm)
{
    return {{"x_left", 0.0, coil_x_min_mm - Real(0.5)},
            {"x_under_coil", coil_x_min_mm, coil_x_max_mm},
            {"x_peak_ref", 108.0, 144.0},
            {"x_high_bias", 144.0, 216.0},
            {"x_right", 216.0, 288.0}};
}

inline Team7BzCompareMetrics compareTeam7BzAgainstReferenceOverXRange(const StdVec<Team7BzProbePoint> &probes,
                                                                      const StdVec<Real> &sim_mT,
                                                                      const StdVec<Real> &ref_mT, Real x_min_mm,
                                                                      Real x_max_mm, Real smoke_rms_threshold = 0.5)
{
    StdVec<Real> sim_sub;
    StdVec<Real> ref_sub;
    StdVec<Team7BzProbePoint> probe_sub;
    for (size_t i = 0; i < probes.size(); ++i)
    {
        if (probes[i].x_mm >= x_min_mm && probes[i].x_mm <= x_max_mm)
        {
            probe_sub.push_back(probes[i]);
            sim_sub.push_back(sim_mT[i]);
            ref_sub.push_back(ref_mT[i]);
        }
    }
    return compareTeam7BzAgainstReference(probe_sub, sim_sub, ref_sub, smoke_rms_threshold);
}

inline StdVec<Team7ProbeXBandMetrics> evaluateTeam7ProbePhase90XBands(
    const StdVec<Team7BzProbePoint> &reference_probes, const StdVec<Team7ProbeBzDecomposition> &decomp,
    const StdVec<Team7ProbeXBandSpec> &bands, Real smoke_rms_threshold = 0.5)
{
    StdVec<Real> phase90_sim;
    StdVec<Real> phase90_ref;
    StdVec<Real> phase90_skin_sim;
    team7ProbePhase90SimAndRef(reference_probes, decomp, phase90_sim, phase90_ref);
    team7ProbePhase90SkinSimAndRef(reference_probes, decomp, phase90_skin_sim, phase90_ref);
    const StdVec<Team7BzProbePoint> coil_probes =
        team7ProbesFromDecompositionReal(reference_probes, decomp, false);

    StdVec<Team7ProbeXBandMetrics> out;
    out.reserve(bands.size());
    for (const Team7ProbeXBandSpec &band : bands)
    {
        Team7ProbeXBandMetrics metrics;
        metrics.band = band;
        for (size_t i = 0; i < reference_probes.size(); ++i)
        {
            if (reference_probes[i].x_mm >= band.x_min_mm && reference_probes[i].x_mm <= band.x_max_mm)
            {
                ++metrics.n_probes;
            }
        }
        metrics.phase0_coil =
            compareTeam7BzPhase0OverXRange(coil_probes, band.x_min_mm, band.x_max_mm, smoke_rms_threshold);
        metrics.phase90_ind_full = compareTeam7BzAgainstReferenceOverXRange(
            reference_probes, phase90_sim, phase90_ref, band.x_min_mm, band.x_max_mm, smoke_rms_threshold);
        metrics.phase90_ind_skin = compareTeam7BzAgainstReferenceOverXRange(
            reference_probes, phase90_skin_sim, phase90_ref, band.x_min_mm, band.x_max_mm, smoke_rms_threshold);
        out.push_back(metrics);
    }
    return out;
}

inline void printTeam7ProbePhase90XBandReport(const StdVec<Team7ProbeXBandMetrics> &bands)
{
    std::cout << "[ophelie] TEAM7 phase90 x-band RMS (phase0_coil / phase90_full / phase90_skin):" << std::endl;
    for (const Team7ProbeXBandMetrics &band_metrics : bands)
    {
        std::cout << "[ophelie]   " << band_metrics.band.name << " x=[" << band_metrics.band.x_min_mm << ","
                  << band_metrics.band.x_max_mm << "] n=" << band_metrics.n_probes
                  << " coil=" << band_metrics.phase0_coil.rms_rel_error
                  << " p90_full=" << band_metrics.phase90_ind_full.rms_rel_error
                  << " p90_skin=" << band_metrics.phase90_ind_skin.rms_rel_error << std::endl;
    }
}

inline void writeTeam7ProbePhase90XBandCsv(const std::string &output_path,
                                           const StdVec<Team7ProbeXBandMetrics> &bands)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 x-band CSV: " << output_path << std::endl;
        return;
    }
    out << "band,x_min_mm,x_max_mm,n_probes,rms_phase0_coil,rms_phase90_full,rms_phase90_skin\n";
    for (const Team7ProbeXBandMetrics &band_metrics : bands)
    {
        out << band_metrics.band.name << "," << band_metrics.band.x_min_mm << "," << band_metrics.band.x_max_mm << ","
            << band_metrics.n_probes << "," << band_metrics.phase0_coil.rms_rel_error << ","
            << band_metrics.phase90_ind_full.rms_rel_error << "," << band_metrics.phase90_ind_skin.rms_rel_error
            << "\n";
    }
    std::cout << "[ophelie] TEAM7 phase90 x-band CSV: " << output_path << std::endl;
}

struct Team7PlateDepthBinStats
{
    Real z_center_mm = 0.0;
    Real z_min_mm = 0.0;
    Real z_max_mm = 0.0;
    size_t n_particles = 0;
    Real j_imag_l2_vol = 0.0;
    Real b_ind_l2_vol = 0.0;
    Real b_coil_l2_vol = 0.0;
    Real bind_over_bcoil = 0.0;
    Real vol_sum_m3 = 0.0;
};

inline StdVec<Team7PlateDepthBinStats> computeTeam7PlateDepthProfile(
    BaseParticles &plate_particles, const OphelieGlassFieldNames &plate_names, const std::string &j_imag_field,
    Real z_lower_m, Real z_upper_m, Real bin_h_m, Real mm_to_m = Real(1.0e-3),
    StdVec<Real> *particle_depth_bin_center_mm = nullptr)
{
    if (bin_h_m <= TinyReal || z_upper_m <= z_lower_m + TinyReal)
    {
        return {};
    }
    const size_t n_bins =
        std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil((z_upper_m - z_lower_m) / bin_h_m)));
    StdVec<Team7PlateDepthBinStats> bins(n_bins);
    StdVec<Real> j_sq_sum(n_bins, 0.0);
    StdVec<Real> b_ind_sq_sum(n_bins, 0.0);
    StdVec<Real> b_coil_sq_sum(n_bins, 0.0);
    StdVec<Real> vol_sum(n_bins, 0.0);
    StdVec<size_t> count(n_bins, 0);

    syncVariableToHost<Vecd>(plate_particles, "Position");
    syncVariableToHost<Vecd>(plate_particles, j_imag_field);
    syncVariableToHost<Vecd>(plate_particles, plate_names.b_ind_imag);
    syncVariableToHost<Vecd>(plate_particles, plate_names.b_coil_real);
    syncVariableToHost<Real>(plate_particles, "VolumetricMeasure");
    const Vecd *pos = plate_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_imag = plate_particles.getVariableDataByName<Vecd>(j_imag_field);
    const Vecd *b_ind = plate_particles.getVariableDataByName<Vecd>(plate_names.b_ind_imag);
    const Vecd *b_coil = plate_particles.getVariableDataByName<Vecd>(plate_names.b_coil_real);
    const Real *vol = plate_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = plate_particles.TotalRealParticles();

    if (particle_depth_bin_center_mm != nullptr)
    {
        particle_depth_bin_center_mm->assign(n, 0.0);
    }

    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        if (z < z_lower_m || z > z_upper_m)
        {
            continue;
        }
        size_t bin = static_cast<size_t>(std::floor((z - z_lower_m) / bin_h_m));
        bin = std::min(bin, n_bins - 1);
        const Real v = vol[i];
        j_sq_sum[bin] += v * j_imag[i].squaredNorm();
        b_ind_sq_sum[bin] += v * b_ind[i].squaredNorm();
        b_coil_sq_sum[bin] += v * b_coil[i].squaredNorm();
        vol_sum[bin] += v;
        ++count[bin];
        if (particle_depth_bin_center_mm != nullptr)
        {
            const Real z_center_m = z_lower_m + (static_cast<Real>(bin) + Real(0.5)) * bin_h_m;
            (*particle_depth_bin_center_mm)[i] = z_center_m / mm_to_m;
        }
    }

    for (size_t b = 0; b < n_bins; ++b)
    {
        bins[b].z_min_mm = (z_lower_m + static_cast<Real>(b) * bin_h_m) / mm_to_m;
        bins[b].z_max_mm = (z_lower_m + static_cast<Real>(b + 1) * bin_h_m) / mm_to_m;
        bins[b].z_center_mm = (z_lower_m + (static_cast<Real>(b) + Real(0.5)) * bin_h_m) / mm_to_m;
        bins[b].n_particles = count[b];
        bins[b].vol_sum_m3 = vol_sum[b];
        bins[b].j_imag_l2_vol = std::sqrt(j_sq_sum[b]);
        bins[b].b_ind_l2_vol = std::sqrt(b_ind_sq_sum[b]);
        bins[b].b_coil_l2_vol = std::sqrt(b_coil_sq_sum[b]);
        bins[b].bind_over_bcoil = bins[b].b_ind_l2_vol / (bins[b].b_coil_l2_vol + TinyReal);
    }
    return bins;
}

inline void hostWriteScalarField(BaseParticles &particles, const std::string &field_name,
                                 const StdVec<Real> &host_values)
{
    syncVariableToHost<Real>(particles, field_name);
    Real *values = particles.getVariableDataByName<Real>(field_name);
    for (size_t i = 0; i < host_values.size(); ++i)
    {
        values[i] = host_values[i];
    }
    syncVariableToDevice<Real>(particles, field_name);
}

inline void evaluatePlateJyAtProbesSphInterp(BaseParticles &plate_particles, const std::string &j_real_field,
                                             const std::string &j_imag_field, Real kernel_h,
                                             const StdVec<Team7JeyProbePoint> &probes,
                                             StdVec<Team7JeyProbePoint> &results)
{
    results = probes;
    syncVariableToHost<Vecd>(plate_particles, "Position");
    syncVariableToHost<Vecd>(plate_particles, j_real_field);
    syncVariableToHost<Vecd>(plate_particles, j_imag_field);
    const Vecd *pos = plate_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_real = plate_particles.getVariableDataByName<Vecd>(j_real_field);
    const Vecd *j_imag = plate_particles.getVariableDataByName<Vecd>(j_imag_field);
    const size_t n_plate = plate_particles.TotalRealParticles();
    for (size_t ip = 0; ip < probes.size(); ++ip)
    {
        Real weight_sum = 0.0;
        Real jy_real_sum = 0.0;
        Real jy_imag_sum = 0.0;
        for (size_t i = 0; i < n_plate; ++i)
        {
            const Real weight = team7ProbeWendlandC4Weight((probes[ip].position_m - pos[i]).norm(), kernel_h);
            if (weight <= TinyReal)
            {
                continue;
            }
            weight_sum += weight;
            jy_real_sum += weight * j_real[i][1];
            jy_imag_sum += weight * j_imag[i][1];
        }
        results[ip].jy_sim_phase0_Am2 = weight_sum > TinyReal ? jy_real_sum / weight_sum : Real(0);
        results[ip].jy_sim_phase90_Am2 = weight_sum > TinyReal ? jy_imag_sum / weight_sum : Real(0);
    }
}

struct Team7JeySurfaceDepthDiagnostic
{
    Real jy_surface_phase0_l2_vol = 0.0;
    Real jy_surface_phase90_l2_vol = 0.0;
    Real jy_mid_depth_phase0_l2_vol = 0.0;
    Real jy_mid_depth_phase90_l2_vol = 0.0;
    Real jy_surface_over_mid_phase0 = 0.0;
    Real jy_surface_over_mid_phase90 = 0.0;
};

inline Team7JeySurfaceDepthDiagnostic computeTeam7JeySurfaceDepthDiagnostic(
    BaseParticles &plate_particles, const OphelieGlassFieldNames &plate_names, const std::string &j_real_field,
    const std::string &j_imag_field, Real z_lower_m, Real z_upper_m, Real surface_band_m, Real mid_band_half_m)
{
    Team7JeySurfaceDepthDiagnostic diag;
    syncVariableToHost<Vecd>(plate_particles, "Position");
    syncVariableToHost<Vecd>(plate_particles, j_real_field);
    syncVariableToHost<Vecd>(plate_particles, j_imag_field);
    syncVariableToHost<Real>(plate_particles, "VolumetricMeasure");
    const Vecd *pos = plate_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_real = plate_particles.getVariableDataByName<Vecd>(j_real_field);
    const Vecd *j_imag = plate_particles.getVariableDataByName<Vecd>(j_imag_field);
    const Real *vol = plate_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = plate_particles.TotalRealParticles();
    const Real z_mid_m = Real(0.5) * (z_lower_m + z_upper_m);
    Real surf_jy0_sq = 0.0;
    Real surf_jyi_sq = 0.0;
    Real surf_vol = 0.0;
    Real mid_jy0_sq = 0.0;
    Real mid_jyi_sq = 0.0;
    Real mid_vol = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        const Real v = vol[i];
        const Real jy0 = j_real[i][1];
        const Real jyi = j_imag[i][1];
        if (z >= z_upper_m - surface_band_m)
        {
            surf_jy0_sq += v * jy0 * jy0;
            surf_jyi_sq += v * jyi * jyi;
            surf_vol += v;
        }
        if (z >= z_mid_m - mid_band_half_m && z <= z_mid_m + mid_band_half_m)
        {
            mid_jy0_sq += v * jy0 * jy0;
            mid_jyi_sq += v * jyi * jyi;
            mid_vol += v;
        }
    }
    diag.jy_surface_phase0_l2_vol = std::sqrt(surf_jy0_sq);
    diag.jy_surface_phase90_l2_vol = std::sqrt(surf_jyi_sq);
    diag.jy_mid_depth_phase0_l2_vol = std::sqrt(mid_jy0_sq);
    diag.jy_mid_depth_phase90_l2_vol = std::sqrt(mid_jyi_sq);
    diag.jy_surface_over_mid_phase0 =
        diag.jy_surface_phase0_l2_vol / (diag.jy_mid_depth_phase0_l2_vol + TinyReal);
    diag.jy_surface_over_mid_phase90 =
        diag.jy_surface_phase90_l2_vol / (diag.jy_mid_depth_phase90_l2_vol + TinyReal);
    (void)surf_vol;
    (void)mid_vol;
    return diag;
}

inline void printTeam7JeyValidationReport(Real frequency_hz, bool jey_reference_ok, bool jey_frequency_ref_ok,
                                          const Team7BzCompareMetrics &jy_phase0_metrics,
                                          const Team7BzCompareMetrics &jy_phase90_metrics,
                                          const Team7JeySurfaceDepthDiagnostic &surface_depth)
{
    std::cout << "[team7] Jey validation f=" << frequency_hz << " Hz reference_loaded=" << (jey_reference_ok ? 1 : 0)
              << " frequency_ref_ok=" << (jey_frequency_ref_ok ? 1 : 0);
    if (jey_reference_ok && jey_frequency_ref_ok)
    {
        std::cout << " Jy_phase0_RMS=" << jy_phase0_metrics.rms_rel_error
                  << " Jy_phase90_RMS=" << jy_phase90_metrics.rms_rel_error
                  << " peak_sim/ref_phase90=" << jy_phase90_metrics.peak_sim_mT << "/"
                  << jy_phase90_metrics.peak_ref_mT;
    }
    else if (jey_reference_ok)
    {
        std::cout << " (skip RMS: no reference values for this frequency)";
    }
    std::cout << " Jy_surface/mid_phase0=" << surface_depth.jy_surface_over_mid_phase0
              << " Jy_surface/mid_phase90=" << surface_depth.jy_surface_over_mid_phase90
              << " (diagnostic)" << std::endl;
}

struct Team7JvsBProbeSplitRecord
{
    Real x_mm = 0.0;
    Real jy_sim_phase90_Am2 = 0.0;
    Real jy_ref_phase90_Am2 = 0.0;
    Real jy_sim_over_ref = 0.0;
    Real bz_ind_sim_phase90_mT = 0.0;
    Real bz_ref_phase90_mT = 0.0;
    Real bz_sim_over_ref = 0.0;
    Real bz_over_jy_ratio_sim = 0.0;
    Real bz_over_jy_ratio_ref = 0.0;
};

inline StdVec<Team7JvsBProbeSplitRecord> computeTeam7JvsBProbeSplit(
    const StdVec<Team7JeyProbePoint> &jey_probes, const StdVec<Team7BzProbePoint> &bz_reference_probes,
    const StdVec<Team7ProbeBzDecomposition> &decomp, bool use_neg_ind_phase90)
{
    StdVec<Team7JvsBProbeSplitRecord> records;
    const size_t n = std::min(jey_probes.size(), std::min(bz_reference_probes.size(), decomp.size()));
    records.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        Team7JvsBProbeSplitRecord record;
        record.x_mm = jey_probes[i].x_mm;
        record.jy_sim_phase90_Am2 = jey_probes[i].jy_sim_phase90_Am2;
        record.jy_ref_phase90_Am2 = jey_probes[i].jy_ref_phase90_Am2;
        record.jy_sim_over_ref = record.jy_sim_phase90_Am2 / (record.jy_ref_phase90_Am2 + TinyReal);
        record.bz_ind_sim_phase90_mT =
            use_neg_ind_phase90 ? decomp[i].bz_ind_imag_neg_mT : decomp[i].bz_total_imag_mT;
        record.bz_ref_phase90_mT = bz_reference_probes[i].bz_ref_phase90_mT;
        record.bz_sim_over_ref = record.bz_ind_sim_phase90_mT / (record.bz_ref_phase90_mT + TinyReal);
        record.bz_over_jy_ratio_sim =
            record.bz_ind_sim_phase90_mT / (std::abs(record.jy_sim_phase90_Am2) * 1.0e-6 + TinyReal);
        record.bz_over_jy_ratio_ref =
            record.bz_ref_phase90_mT / (std::abs(record.jy_ref_phase90_Am2) * 1.0e-6 + TinyReal);
        records.push_back(record);
    }
    return records;
}

inline void printTeam7JvsBProbeSplitReport(const StdVec<Team7JvsBProbeSplitRecord> &records,
                                           Real highlight_x_mm = Real(162))
{
    Real j_ratio_median = 0.0;
    Real b_ratio_median = 0.0;
    size_t n_valid = 0;
    for (const Team7JvsBProbeSplitRecord &record : records)
    {
        if (std::abs(record.jy_ref_phase90_Am2) > TinyReal && std::abs(record.bz_ref_phase90_mT) > TinyReal)
        {
            j_ratio_median += std::abs(record.jy_sim_over_ref);
            b_ratio_median += std::abs(record.bz_sim_over_ref);
            ++n_valid;
        }
    }
    if (n_valid > 0)
    {
        j_ratio_median /= Real(n_valid);
        b_ratio_median /= Real(n_valid);
    }
    std::cout << "[team7] J vs B probe split (imag chain): median_|J_sim/J_ref|=" << j_ratio_median
              << " median_|B_ind_sim/B_ref|=" << b_ratio_median
              << " (if J~ref and B>>ref => Biot/probe; if both >>ref => J-driven)" << std::endl;
    for (const Team7JvsBProbeSplitRecord &record : records)
    {
        if (std::abs(record.x_mm - highlight_x_mm) < Real(0.5))
        {
            std::cout << "[team7]   @x=" << record.x_mm << "mm Jy_sim/ref=" << record.jy_sim_phase90_Am2 << "/"
                      << record.jy_ref_phase90_Am2 << " ratio=" << record.jy_sim_over_ref
                      << " B_ind_sim/ref=" << record.bz_ind_sim_phase90_mT << "/" << record.bz_ref_phase90_mT
                      << " ratio=" << record.bz_sim_over_ref << " (B/J)_sim=" << record.bz_over_jy_ratio_sim
                      << " (B/J)_ref=" << record.bz_over_jy_ratio_ref << std::endl;
            break;
        }
    }
}

inline void writeTeam7JvsBProbeSplitCsv(const std::string &output_path,
                                        const StdVec<Team7JvsBProbeSplitRecord> &records)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        return;
    }
    out << "x_mm,Jy_sim_phase90_Am2,Jy_ref_phase90_Am2,Jy_sim_over_ref,B_ind_sim_phase90_mT,B_ref_phase90_mT,"
           "B_sim_over_ref,B_over_J_sim,B_over_J_ref\n";
    for (const Team7JvsBProbeSplitRecord &record : records)
    {
        out << record.x_mm << "," << record.jy_sim_phase90_Am2 << "," << record.jy_ref_phase90_Am2 << ","
            << record.jy_sim_over_ref << "," << record.bz_ind_sim_phase90_mT << "," << record.bz_ref_phase90_mT << ","
            << record.bz_sim_over_ref << "," << record.bz_over_jy_ratio_sim << "," << record.bz_over_jy_ratio_ref
            << "\n";
    }
    std::cout << "[ophelie] TEAM7 J-vs-B probe split CSV: " << output_path << std::endl;
}

inline void printTeam7PlateDepthProfileReport(const StdVec<Team7PlateDepthBinStats> &bins)
{
    std::cout << "[ophelie] TEAM7 plate depth profile (z_mm, n, J_imag_L2_vol, Bind/Bcoil):" << std::endl;
    for (const Team7PlateDepthBinStats &bin : bins)
    {
        if (bin.n_particles == 0)
        {
            continue;
        }
        std::cout << "[ophelie]   z=" << bin.z_center_mm << " n=" << bin.n_particles << " J_imag_L2=" << bin.j_imag_l2_vol
                  << " Bind/Bcoil=" << bin.bind_over_bcoil << std::endl;
    }
}

inline void writeTeam7PlateDepthProfileCsv(const std::string &output_path,
                                           const StdVec<Team7PlateDepthBinStats> &bins)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 depth profile CSV: " << output_path << std::endl;
        return;
    }
    out << "z_center_mm,z_min_mm,z_max_mm,n_particles,vol_sum_m3,J_imag_L2_vol,B_ind_L2_vol,B_coil_L2_vol,"
           "Bind_over_Bcoil\n";
    for (const Team7PlateDepthBinStats &bin : bins)
    {
        out << bin.z_center_mm << "," << bin.z_min_mm << "," << bin.z_max_mm << "," << bin.n_particles << ","
            << bin.vol_sum_m3 << "," << bin.j_imag_l2_vol << "," << bin.b_ind_l2_vol << "," << bin.b_coil_l2_vol << ","
            << bin.bind_over_bcoil << "\n";
    }
    std::cout << "[ophelie] TEAM7 plate depth profile CSV: " << output_path << std::endl;
}

inline void appendTeam7ValidationXBandRecord(const std::string &output_path,
                                             const StdVec<Team7ProbeXBandMetrics> &bands)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path, std::ios::app);
    if (!out)
    {
        return;
    }
    out << "phase90_x_bands:\n";
    for (const Team7ProbeXBandMetrics &band_metrics : bands)
    {
        out << "  " << band_metrics.band.name << " n=" << band_metrics.n_probes
            << " coil=" << band_metrics.phase0_coil.rms_rel_error
            << " p90_full=" << band_metrics.phase90_ind_full.rms_rel_error
            << " p90_skin=" << band_metrics.phase90_ind_skin.rms_rel_error << "\n";
    }
    out << "\n";
}

inline void writeTeam7LevelProbeCsv(const std::string &output_path, const StdVec<Team7BzProbePoint> &reference_probes,
                                    const StdVec<Team7ProbeBzDecomposition> &decomp, bool include_reference)
{
    if (decomp.empty() || reference_probes.empty())
    {
        std::cout << "[team7][warning] probe CSV skipped (empty decomposition or reference probes): " << output_path
                  << std::endl;
        return;
    }
    if (include_reference && decomp.size() != reference_probes.size())
    {
        std::cout << "[team7][warning] probe CSV skipped (reference/decomp size mismatch): " << output_path
                  << std::endl;
        return;
    }
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 probe CSV: " << output_path << std::endl;
        return;
    }
    out << "x_mm";
    if (include_reference)
    {
        out << ",Bz_ref_phase0_mT";
    }
    out << ",Bz_coil_mT,Bz_ind_real_mT,Bz_ind_imag_mT,Bz_ind_imag_skin_mT,Bz_ind_imag_particle_mT,Bz_total_real_mT,"
           "Bz_total_imag_mT";
    if (include_reference)
    {
        out << ",err_mT,rel_err";
    }
    out << "\n";
    for (size_t i = 0; i < decomp.size(); ++i)
    {
        out << reference_probes[i].x_mm;
        if (include_reference)
        {
            out << "," << reference_probes[i].bz_ref_phase0_mT;
        }
        out << "," << decomp[i].bz_coil_mT << "," << decomp[i].bz_ind_real_mT << "," << decomp[i].bz_ind_imag_mT
            << "," << decomp[i].bz_ind_imag_skin_mT << "," << decomp[i].bz_ind_imag_particle_mT << ","
            << decomp[i].bz_total_real_mT << "," << decomp[i].bz_total_imag_mT;
        if (include_reference)
        {
            const Real err = decomp[i].bz_total_real_mT - reference_probes[i].bz_ref_phase0_mT;
            const Real rel = err / (std::abs(reference_probes[i].bz_ref_phase0_mT) + TinyReal);
            out << "," << err << "," << rel;
        }
        out << "\n";
    }
    std::cout << "[ophelie] TEAM7 probe CSV: " << output_path << std::endl;
}

struct Team7ImagChainPowerAudit
{
    Real j_imag_vol_norm = 0.0;
    Real e_imag_vol_norm = 0.0;
    Real j_e_vol_dot = 0.0;
    Real joule_power_vol_w = 0.0;
    Real j_e_cos_angle = 0.0;
};

inline Real hostVecdVolWeightedDot(BaseParticles &particles, const std::string &lhs_field,
                                   const std::string &rhs_field, size_t n)
{
    syncVariableToHost<Vecd>(particles, lhs_field);
    syncVariableToHost<Vecd>(particles, rhs_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *lhs = particles.getVariableDataByName<Vecd>(lhs_field);
    const Vecd *rhs = particles.getVariableDataByName<Vecd>(rhs_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t k = 0; k < n; ++k)
    {
        sum += vol[k] * lhs[k].dot(rhs[k]);
    }
    return sum;
}

inline Team7ImagChainPowerAudit auditTeam7ImagChainPowerOnPlate(BaseParticles &particles,
                                                                const OphelieGlassFieldNames &names,
                                                                const std::string &j_imag_field, size_t n)
{
    Team7ImagChainPowerAudit audit;
    audit.j_imag_vol_norm = hostVecdVolWeightedNorm(particles, j_imag_field, n);
    audit.e_imag_vol_norm = hostVecdVolWeightedNorm(particles, names.e_imag, n);
    audit.j_e_vol_dot = hostVecdVolWeightedDot(particles, j_imag_field, names.e_imag, n);
    audit.joule_power_vol_w = std::max(audit.j_e_vol_dot, Real(0));
    const Real denom = audit.j_imag_vol_norm * audit.e_imag_vol_norm + TinyReal;
    audit.j_e_cos_angle = audit.j_e_vol_dot / denom;
    return audit;
}

inline void printTeam7ImagChainPowerAudit(const Team7ImagChainPowerAudit &audit)
{
    std::cout << "[ophelie] TEAM7 imag-chain audit: |J_imag|_vol=" << audit.j_imag_vol_norm
              << " |E_imag|_vol=" << audit.e_imag_vol_norm << " Re(J·E)_vol=" << audit.j_e_vol_dot
              << " P_vol_W=" << audit.joule_power_vol_w << " cos(J,E)=" << audit.j_e_cos_angle << std::endl;
}

struct Team7ImagJEmVsEdgeAudit
{
    Real j_edge_vol_norm = 0.0;
    Real j_em_vol_norm = 0.0;
    Real j_em_over_edge = 0.0;
    Real bz_ind_edge_mT = 0.0;
    Real bz_ind_em_mT = 0.0;
    Real bz_ind_em_neg_mT = 0.0;
    Real ref_phase90_mT = 0.0;
    Real probe_x_mm = 0.0;
};

inline Real hostBzImagBiotFromCustomJAtProbe(BaseParticles &particles, const StdVec<Vecd> &j_imag, const Vec3d &probe_pos,
                                             Real mu0, Real softening_length)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real coeff = mu0 / (4.0 * Pi);
    const Real eps2 = softening_length * softening_length;
    Vecd b_imag = Vecd::Zero();
    for (size_t j = 0; j < n; ++j)
    {
        const Vecd r = probe_pos - pos[j];
        const Real r2 = r.squaredNorm() + eps2;
        const Real inv_r3 = 1.0 / (std::sqrt(r2) * r2);
        b_imag += coeff * (j_imag[j] * vol[j]).cross(r) * inv_r3;
    }
    return b_imag[2] * 1000.0;
}

inline Team7ImagJEmVsEdgeAudit auditTeam7ImagJEmVsEdgeAtProbe(BaseParticles &particles,
                                                              const OphelieGlassFieldNames &names,
                                                              const OphelieParameters &params,
                                                              const std::string &j_edge_field,
                                                              const Team7BzProbePoint &probe)
{
    Team7ImagJEmVsEdgeAudit audit;
    audit.probe_x_mm = probe.x_mm;
    audit.ref_phase90_mT = probe.bz_ref_phase90_mT;
    const size_t n = particles.TotalRealParticles();
    audit.j_edge_vol_norm = hostVecdVolWeightedNorm(particles, j_edge_field, n);

    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Real>(particles, names.sigma);
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real omega = params.omega();
    StdVec<Vecd> j_em(n);
    Real j_em_sq = 0.0;
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t k = 0; k < n; ++k)
    {
        const Vecd e_em = -grad_phi[k] - omega * a_coil_real[k];
        j_em[k] = sigma[k] * e_em;
        j_em_sq += vol[k] * j_em[k].squaredNorm();
    }
    audit.j_em_vol_norm = std::sqrt(j_em_sq);
    audit.j_em_over_edge = audit.j_em_vol_norm / (audit.j_edge_vol_norm + TinyReal);

    StdVec<Real> bz_edge;
    StdVec<Real> bz_edge_unused;
    StdVec<Team7BzProbePoint> probes(1, probe);
    evaluateVolumeCurrentBiotSavartBzAtProbes(particles, names.j_edge_recon_real, j_edge_field, params.mu0_,
                                              params.softening_length_, probes, bz_edge_unused, bz_edge,
                                              names.edge_recon_fallback);
    audit.bz_ind_edge_mT = bz_edge.empty() ? Real(0) : bz_edge[0];
    audit.bz_ind_em_mT = hostBzImagBiotFromCustomJAtProbe(particles, j_em, probe.position_m, params.mu0_,
                                                          params.softening_length_);
    audit.bz_ind_em_neg_mT = -audit.bz_ind_em_mT;
    return audit;
}

inline void printTeam7ImagJEmVsEdgeAudit(const Team7ImagJEmVsEdgeAudit &audit)
{
    std::cout << "[ophelie] TEAM7 J_em vs J_edge @x=" << audit.probe_x_mm
              << "mm: |J_edge|_vol=" << audit.j_edge_vol_norm << " |J_em|_vol=" << audit.j_em_vol_norm
              << " J_em/edge=" << audit.j_em_over_edge << " Bz_ind_edge=" << audit.bz_ind_edge_mT
              << " Bz_ind_em=" << audit.bz_ind_em_mT << " Bz_ind_-em=" << audit.bz_ind_em_neg_mT
              << " ref_p90=" << audit.ref_phase90_mT << std::endl;
}

/** Diagnostic post-scale: scale induced J/E (not coil fields), refresh A_ind/B_ind. B_ind scales ~linearly with J. */
template <class ExecutionPolicy>
inline Real applyTeam7InducedJPostScaleAndRefreshAInd(SolidBody &plate_body, const OphelieGlassFieldNames &names,
                                                      const OphelieParameters &params, Real j_scale, size_t n)
{
    if (std::abs(j_scale - 1.0) <= TinyReal)
    {
        return 1.0;
    }
    BaseParticles &particles = plate_body.getBaseParticles();
    const Real power_scale = j_scale * j_scale;
    const std::string &j_real = getOphelieAIndJRealFieldName(names, params);
    const std::string &j_imag = getOphelieAIndJImagFieldName(names, params);

    hostScaleVecdFieldInPlace(particles, j_imag, j_scale, n);
    hostScaleVecdFieldInPlace(particles, j_real, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.e_imag, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.e_real, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.j_edge_recon_imag, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.j_edge_recon_real, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.e_edge_recon_imag, j_scale, n);
    hostScaleVecdFieldInPlace(particles, names.e_edge_recon_real, j_scale, n);
    hostScaleScalarFieldInPlace(particles, names.phi_imag, j_scale, n);
    hostScaleScalarFieldInPlace(particles, names.joule_heat, power_scale, n);
    if (params.edge_flux_complex_)
    {
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_imag, power_scale, n);
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_real, power_scale, n);
        hostScaleScalarFieldInPlace(particles, names.joule_heat_edge_recon_complex, power_scale, n);
    }

    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> refresh_aind(
        plate_body, names, params, j_real, j_imag);
    refresh_aind.exec();
    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine(plate_body, names);
    combine.exec();
    syncGlassElectromagneticFieldsToHost(particles, names);

    const Real b_ind_vol = hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n);
    const Real b_coil_vol = hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n);
    std::cout << "[ophelie] TEAM7 J_imag post-scale=" << j_scale << " Bind/Bcoil_after=" << b_ind_vol / (b_coil_vol + TinyReal)
              << std::endl;
    return j_scale;
}

struct Team7PicardIterationRecord
{
    size_t iter = 0;
    std::string normalization_mode;
    Real normalization_scale = 1.0;
    Real j_rel_raw = 0.0;
    Real j_rel_relaxed = 0.0;
    Real aind_rel_raw = 0.0;
    Real aind_rel_relaxed = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real p_recon_w = 0.0;
    Real a_ind_over_a_coil = 0.0;
    Real bind_over_bcoil = 0.0;
    Real bz_phase0_rms = 0.0;
    Real bz_phase90_rms = 0.0;
    Real max_j_real = 0.0;
    Real max_j_imag = 0.0;
    Real phi_eq_res_imag = 0.0;
    Real phi_eq_res_real = 0.0;
    Real edge_res_red_imag = 0.0;
    Real edge_res_red_real = 0.0;
};

struct Team7PicardRunResult
{
    Real final_j_rel_change = 0.0;
    size_t iterations_used = 0;
    Real phi_solver_rel_residual = 0.0;
    Real phi_eq_res_vol = 0.0;
    bool picard_converged = false;
    Team7OneWayEdgeFluxSummary final_summary;
    StdVec<Team7PicardIterationRecord> iteration_log;
};

inline void writeTeam7PicardIterationCsv(const std::string &output_path,
                                         const StdVec<Team7PicardIterationRecord> &records)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 Picard CSV: " << output_path << std::endl;
        return;
    }
    out << "iter,normalization_mode,normalization_scale,j_rel_raw,j_rel_relaxed,aind_rel_raw,aind_rel_relaxed,"
           "J_imag_vol_norm,P_recon,Aind_over_Acoil,Bind_over_Bcoil,Bz_phase0_RMS,Bz_phase90_RMS,max_J_real,max_J_imag,"
           "phi_eq_res_imag,phi_eq_res_real,edge_res_red_imag,edge_res_red_real\n";
    out << std::setprecision(10);
    for (const Team7PicardIterationRecord &record : records)
    {
        out << record.iter << "," << record.normalization_mode << "," << record.normalization_scale << ","
            << record.j_rel_raw << "," << record.j_rel_relaxed << "," << record.aind_rel_raw << ","
            << record.aind_rel_relaxed << "," << record.j_imag_vol_norm << "," << record.p_recon_w << ","
            << record.a_ind_over_a_coil << "," << record.bind_over_bcoil << "," << record.bz_phase0_rms << ","
            << record.bz_phase90_rms << "," << record.max_j_real << "," << record.max_j_imag << ","
            << record.phi_eq_res_imag << "," << record.phi_eq_res_real << "," << record.edge_res_red_imag << ","
            << record.edge_res_red_real << "\n";
    }
    std::cout << "[ophelie] TEAM7 Picard iteration CSV: " << output_path << std::endl;
}

inline void printTeam7PicardIterationReport(const StdVec<Team7PicardIterationRecord> &records)
{
    std::cout << "[team7] Picard iteration log (iter Aind_rel J_rel Bind/B Bz_p90 phi_eq_imag):"
              << std::endl;
    for (const Team7PicardIterationRecord &record : records)
    {
        std::cout << "[team7]   iter=" << record.iter << " Aind_rel_raw=" << record.aind_rel_raw
                  << " Aind_rel=" << record.aind_rel_relaxed << " J_rel_raw=" << record.j_rel_raw
                  << " J_rel=" << record.j_rel_relaxed << " Bind/Bcoil=" << record.bind_over_bcoil
                  << " Bz_phase90_RMS=" << record.bz_phase90_rms << " phi_eq_imag=" << record.phi_eq_res_imag
                  << std::endl;
    }
}

template <class ExecutionPolicy>
inline Team7PicardIterationRecord collectTeam7PicardIterationRecord(
    size_t iter, Real j_rel_raw, Real j_rel_relaxed, Real aind_rel_raw, Real aind_rel_relaxed,
    const OphelieComplexEdgeFluxSolveReport &solve_report,
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const StdVec<Team7BzProbePoint> &reference_probes, bool reference_ok,
    BaseParticles &coil_particles, const OphelieCoilFieldNames &coil_names,
    const StdVec<OphelieCurrentMomentSample> *filament_moments, Real plate_z_skin_min_m, Real plate_z_top_m,
    Real plate_z_mid_m, Real dp, Team7Phase90Convention phase90_convention)
{
    Team7PicardIterationRecord record;
    record.iter = iter;
    record.normalization_mode = ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_);
    record.normalization_scale = params.edge_flux_solver_local_rhs_scale_;
    record.j_rel_raw = j_rel_raw;
    record.j_rel_relaxed = j_rel_relaxed;
    record.aind_rel_raw = aind_rel_raw;
    record.aind_rel_relaxed = aind_rel_relaxed;
    record.phi_eq_res_imag = solve_report.phi_eq_res_vol_imag;
    record.phi_eq_res_real = solve_report.phi_eq_res_vol_real;

    const size_t n = particles.TotalRealParticles();
    syncGlassElectromagneticFieldsToHost(particles, names);
    const std::string &j_real = getOphelieAIndJRealFieldName(names, params);
    const std::string &j_imag = getOphelieAIndJImagFieldName(names, params);
    const Real a_coil_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_coil_real, names.a_coil_imag, n);
    const Real a_ind_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_ind_real, names.a_ind_imag, n);
    const Real b_coil_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n);
    const Real b_ind_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n);
    record.a_ind_over_a_coil = a_ind_norm / (a_coil_norm + TinyReal);
    record.bind_over_bcoil = b_ind_norm / (b_coil_norm + TinyReal);
    record.max_j_real = hostVecdFieldMax(particles, j_real, n);
    record.max_j_imag = hostVecdFieldMax(particles, j_imag, n);
    record.j_imag_vol_norm = hostVecdVolWeightedNorm(particles, j_imag, n);
    record.p_recon_w = hostEdgeFluxReconPower(particles, names, n, params);

    if (reference_ok && !reference_probes.empty())
    {
        const StdVec<Team7ProbeBzDecomposition> decomp = evaluateTeam7ProbeBzDecomposition(
            coil_particles, coil_names, particles, names, params, j_real, j_imag, reference_probes, plate_z_skin_min_m,
            plate_z_top_m, plate_z_mid_m, dp, filament_moments);
        const StdVec<Team7BzProbePoint> total_probes = team7ProbesFromDecompositionReal(reference_probes, decomp, true);
        const Team7BzCompareMetrics phase0_metrics = compareTeam7BzPhase0(total_probes, Real(1.0));
        record.bz_phase0_rms = phase0_metrics.rms_rel_error;
        StdVec<Real> phase90_sim_mT;
        StdVec<Real> phase90_ref_mT;
        const bool use_neg_ind = team7Phase90ConventionUsesNegInd(phase90_convention);
        team7ProbePhase90SimAndRef(reference_probes, decomp, phase90_sim_mT, phase90_ref_mT, use_neg_ind);
        const Team7BzCompareMetrics phase90_metrics =
            compareTeam7BzAgainstReference(reference_probes, phase90_sim_mT, phase90_ref_mT, Real(1.0));
        record.bz_phase90_rms = phase90_metrics.rms_rel_error;
    }
    return record;
}

template <class ExecutionPolicy>
inline Team7OneWayEdgeFluxSummary team7OneWaySummaryFromPicardRecord(const Team7PicardIterationRecord &record,
                                                                     Real p_complex_total)
{
    Team7OneWayEdgeFluxSummary summary;
    summary.p_complex_coil_only = 0.0;
    summary.p_complex_total = p_complex_total;
    summary.phi_eq_res_vol_imag = record.phi_eq_res_imag;
    summary.phi_eq_res_vol_real = record.phi_eq_res_real;
    summary.a_ind_over_a_coil = record.a_ind_over_a_coil;
    summary.b_ind_over_b_coil = record.bind_over_bcoil;
    summary.max_j_real = record.max_j_real;
    summary.max_j_imag = record.max_j_imag;
    summary.joule_power_recon_w = record.p_recon_w;
    return summary;
}

template <class ExecutionPolicy>
inline Team7PicardRunResult runTeam7ComplexEdgeFluxPicardWithLog(
    SolidBody &plate_body, Inner<> &plate_inner, const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
    Real dp, BaseParticles &coil_particles, const OphelieCoilFieldNames &coil_names,
    const StdVec<Team7BzProbePoint> &reference_probes, bool reference_ok, Team7Phase90Convention phase90_convention,
    const StdVec<OphelieCurrentMomentSample> *filament_moments, Real plate_z_skin_min_m, Real plate_z_top_m,
    Real plate_z_mid_m, const std::string &csv_path)
{
    Team7PicardRunResult result;
    BaseParticles &particles = plate_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const std::string &j_real_source = getOphelieAIndJRealFieldName(plate_names, params);
    const std::string &j_imag_source = getOphelieAIndJImagFieldName(plate_names, params);
    const bool relax_aind = params.self_induction_relax_aind_;
    const Real relax_alpha = params.self_induction_relaxation_factor_;

    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        plate_body, plate_names);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> compute_self_induced_biot(
        plate_body, plate_names, params, j_real_source, j_imag_source);

    StdVec<Real> previous_j_real_components;
    StdVec<Real> previous_j_imag_components;
    StdVec<Real> previous_a_ind_real_components;
    StdVec<Real> previous_a_ind_imag_components;
    hostZeroVecdField(particles, plate_names.a_ind_real, n);
    hostZeroVecdField(particles, plate_names.a_ind_imag, n);
    hostStoreVecdFieldComponents<ExecutionPolicy>(particles, plate_names.a_ind_real, previous_a_ind_real_components,
                                                  n);
    hostStoreVecdFieldComponents<ExecutionPolicy>(particles, plate_names.a_ind_imag, previous_a_ind_imag_components,
                                                  n);

    Real relative_j_change = 0.0;
    Real relative_aind_change = 0.0;
    Real j_rel_raw = 0.0;
    Real aind_rel_raw = 0.0;
    result.iterations_used = 0;
    OphelieProgressLogger self_ind_progress("team7_picard");
    self_ind_progress.log(std::string("physical-scale Picard relax=") + (relax_aind ? "A_ind" : "J") +
                          " alpha=" + std::to_string(relax_alpha) + " norm=" +
                          ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_));

    const bool saved_use_a_total = params.use_a_total_for_edge_flux_;
    params.use_a_total_for_edge_flux_ = true;

    OphelieComplexEdgeFluxSolveReport solve_report;
    Real last_p_complex_total = 0.0;

    for (size_t iteration = 0; iteration < params.self_induction_max_iterations_; ++iteration)
    {
        result.iterations_used = iteration + 1;
        combine_vector_potential.exec();
        last_p_complex_total = execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(
            plate_body, plate_inner, plate_names, params, nullptr, dp, &solve_report);
        result.phi_solver_rel_residual = solve_report.phi_imag_solver_rel_residual;
        result.phi_eq_res_vol = ophelieSelfInductionPicardPhiEqResVol(solve_report);
        compute_self_induced_biot.exec();

        j_rel_raw = 0.0;
        aind_rel_raw = 0.0;
        if (relax_aind)
        {
            if (iteration > 0)
            {
                j_rel_raw = hostRelativeComplexEdgeJChange<ExecutionPolicy>(
                    plate_names, params, particles, previous_j_real_components, previous_j_imag_components, n);
                aind_rel_raw = hostRelativeComplexAIndChange<ExecutionPolicy>(
                    particles, plate_names, previous_a_ind_real_components, previous_a_ind_imag_components, n);
            }
            hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, plate_names.a_ind_real,
                                                                previous_a_ind_real_components, relax_alpha, n);
            hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, plate_names.a_ind_imag,
                                                                previous_a_ind_imag_components, relax_alpha, n);
            relative_aind_change = iteration > 0
                                       ? hostRelativeComplexAIndChange<ExecutionPolicy>(
                                             particles, plate_names, previous_a_ind_real_components,
                                             previous_a_ind_imag_components, n)
                                       : Real(0);
            relative_j_change = j_rel_raw;
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, plate_names.a_ind_real,
                                                          previous_a_ind_real_components, n);
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, plate_names.a_ind_imag,
                                                          previous_a_ind_imag_components, n);
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_real_source, previous_j_real_components, n);
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components, n);
            combine_vector_potential.exec();
        }
        else
        {
            if (iteration > 0)
            {
                j_rel_raw = hostRelativeComplexEdgeJChange<ExecutionPolicy>(
                    plate_names, params, particles, previous_j_real_components, previous_j_imag_components, n);
            }
            hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, j_real_source, previous_j_real_components,
                                                              relax_alpha, n);
            hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components,
                                                              relax_alpha, n);
            relative_j_change = iteration > 0
                                    ? hostRelativeComplexEdgeJChange<ExecutionPolicy>(
                                          plate_names, params, particles, previous_j_real_components,
                                          previous_j_imag_components, n)
                                    : Real(0);
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_real_source, previous_j_real_components, n);
            hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components, n);
        }

        const Real convergence_metric = relax_aind ? relative_aind_change : relative_j_change;
        result.iteration_log.push_back(collectTeam7PicardIterationRecord<ExecutionPolicy>(
            iteration, j_rel_raw, relative_j_change, aind_rel_raw, relative_aind_change, solve_report, particles,
            plate_names, params, reference_probes, reference_ok, coil_particles, coil_names, filament_moments,
            plate_z_skin_min_m, plate_z_top_m, plate_z_mid_m, dp, phase90_convention));
        result.picard_converged =
            iteration > 0 && ophelieSelfInductionPicardConverged(convergence_metric, result.phi_eq_res_vol, params);
        self_ind_progress.log("outer iter " + std::to_string(result.iterations_used) +
                              (relax_aind ? " Aind_rel_raw=" + std::to_string(aind_rel_raw) +
                                                " Aind_rel=" + std::to_string(relative_aind_change)
                                          : " J_rel_raw=" + std::to_string(j_rel_raw) +
                                                " J_rel=" + std::to_string(relative_j_change)) +
                              " Bind/Bcoil=" + std::to_string(result.iteration_log.back().bind_over_bcoil) +
                              " Bz_p90_RMS=" + std::to_string(result.iteration_log.back().bz_phase90_rms) +
                              " picard_converged=" + std::to_string(result.picard_converged ? 1 : 0));
        if (result.picard_converged)
        {
            break;
        }
    }

    params.use_a_total_for_edge_flux_ = saved_use_a_total;
    result.final_j_rel_change = relax_aind ? relative_aind_change : relative_j_change;
    result.picard_converged =
        ophelieSelfInductionPicardConverged(relative_j_change, result.phi_eq_res_vol, params);
    if (!result.iteration_log.empty())
    {
        result.final_summary =
            team7OneWaySummaryFromPicardRecord<ExecutionPolicy>(result.iteration_log.back(), last_p_complex_total);
        result.final_summary.j_imag_vol_norm = hostVecdVolWeightedNorm(particles, j_imag_source, n);
    }
    self_ind_progress.finish("final_J_rel=" + std::to_string(result.final_j_rel_change) + " Bind/Bcoil=" +
                             std::to_string(result.final_summary.b_ind_over_b_coil) +
                             " picard_converged=" + std::to_string(result.picard_converged ? 1 : 0));
    if (!csv_path.empty())
    {
        writeTeam7PicardIterationCsv(csv_path, result.iteration_log);
    }
    printTeam7PicardIterationReport(result.iteration_log);
    std::cout << "[ophelie] TEAM7 L3 Picard final: iters=" << result.iterations_used
              << " J_rel=" << result.final_j_rel_change << " phi_eq_res_vol=" << result.phi_eq_res_vol
              << " converged=" << (result.picard_converged ? 1 : 0)
              << " Bind/Bcoil=" << result.final_summary.b_ind_over_b_coil
              << " Bz_phase90_RMS=" << result.iteration_log.back().bz_phase90_rms << std::endl;
    return result;
}

template <class ExecutionPolicy>
inline Team7OneWayEdgeFluxSummary runTeam7ComplexEdgeFluxOneWay(SolidBody &plate_body, Inner<> &plate_inner,
                                                                const OphelieGlassFieldNames &plate_names,
                                                                OphelieParameters &params, Real dp)
{
    const OphelieAIndOneWayDiagnostic diag = runOphelieEdgeFluxAIndOneWayAfterCoilBiot<ExecutionPolicy>(
        plate_body, plate_inner, plate_names, params, dp);

    Team7OneWayEdgeFluxSummary summary;
    summary.p_complex_coil_only = diag.p_complex_coil_only;
    summary.p_complex_total = diag.feedback_resolve_done ? diag.p_complex_total_a : 0.0;
    summary.phi_eq_res_vol_imag = diag.feedback_resolve_done ? diag.phi_eq_res_vol_total : diag.phi_eq_res_vol;
    summary.phi_eq_res_vol_real = diag.phi_real_solver_rel_residual;
    summary.a_ind_over_a_coil = diag.a_ind_over_a_coil;
    summary.b_ind_over_b_coil = diag.b_ind_over_b_coil;
    summary.max_j_real = diag.feedback_resolve_done ? diag.max_j_real_after_feedback : diag.max_j_real;
    summary.max_j_imag = diag.feedback_resolve_done ? diag.max_j_imag_after_feedback : diag.max_j_imag;
    summary.edge_flux_input_scale = diag.edge_flux_input_scale;
    summary.edge_flux_input_scale_measured = diag.edge_flux_input_scale_measured;
    summary.edge_flux_rhs_l2_scaled = diag.edge_flux_rhs_l2_scaled;
    summary.j_imag_vol_norm = diag.j_imag_vol_norm;
    summary.j_imag_vol_pre_restore = diag.j_imag_vol_pre_restore;
    summary.joule_power_recon_w = diag.joule_power_recon_w;
    summary.joule_power_pre_restore_w = diag.joule_power_pre_restore_w;
    summary.b_ind_over_b_coil_pre_restore = diag.b_ind_over_b_coil_pre_restore;
    summary.edge_flux_rhs_l2_pre_norm = diag.edge_flux_rhs_l2_pre_norm;
    summary.e_edge_em_vol_mismatch = diag.e_edge_em_vol_mismatch;
    summary.j_restore_linearity = diag.j_restore_linearity;
    summary.restore_invariance_error_j = diag.restore_invariance_error_j;
    summary.restore_invariance_error_b = diag.restore_invariance_error_b;
    summary.restore_invariance_error_p = diag.restore_invariance_error_p;
    summary.p_graph_over_recon = diag.p_graph_over_recon;
    summary.grad_phi_imag_vol_norm = diag.grad_phi_imag_vol_norm;
    summary.omega_a_coil_vol_norm = diag.omega_a_coil_vol_norm;
    summary.e_imag_vol_norm = diag.e_imag_vol_norm;

    std::cout << "[ophelie] TEAM7 one-way edge-flux P_coil_only=" << summary.p_complex_coil_only
              << " P_total=" << summary.p_complex_total << " phi_eq_res_imag=" << summary.phi_eq_res_vol_imag
              << " phi_eq_res_real=" << summary.phi_eq_res_vol_real << " Aind/Acoil=" << summary.a_ind_over_a_coil
              << " Bind/Bcoil=" << summary.b_ind_over_b_coil               << " max_J_real=" << summary.max_j_real << " max_J_imag=" << summary.max_j_imag
              << " edge_flux_scale=" << summary.edge_flux_input_scale << " J_imag_L2_vol=" << summary.j_imag_vol_norm
              << " P_recon_W=" << summary.joule_power_recon_w << " rhs_l2_pre=" << summary.edge_flux_rhs_l2_pre_norm
              << " omegaA/gradPhi="
              << (summary.omega_a_coil_vol_norm / (summary.grad_phi_imag_vol_norm + TinyReal))
              << " e_edge_em_mis=" << summary.e_edge_em_vol_mismatch << std::endl;
    return summary;
}

template <class ExecutionPolicy>
inline StdVec<Team7EdgeFluxOmegaScalingRecord> runTeam7EdgeFluxOmegaScalingSweep(
    SolidBody &plate_body, Inner<> &plate_inner, const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
    Real dp, const StdVec<Real> &frequencies_hz, Real sigma,
    const StdVec<OphelieCurrentMomentSample> &filament_moments, bool use_filament_source,
    SolidBody &coil_body, const OphelieCoilFieldNames &coil_names,
    const OphelieTeam7CoilPathPrepareSummary &coil_path_summary)
{
    StdVec<Team7EdgeFluxOmegaScalingRecord> records;
    records.reserve(frequencies_hz.size());
    const Real saved_frequency = params.frequency_;
    for (Real frequency_hz : frequencies_hz)
    {
        params.frequency_ = frequency_hz;
        if (use_filament_source)
        {
            applyOphelieTeam7FilamentRacetrackBiotToGlass(plate_body, plate_names, filament_moments, params.mu0_,
                                                          params.softening_length_);
        }
        else
        {
            StateDynamics<ExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> init_coil_source(
                coil_body, coil_names, coil_path_summary.j0);
            init_coil_source.exec();
            StateDynamics<ExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
                plate_body, coil_body, plate_names, coil_names, params);
            compute_biot.exec();
        }
        const Team7OneWayEdgeFluxSummary summary =
            runTeam7ComplexEdgeFluxOneWay<ExecutionPolicy>(plate_body, plate_inner, plate_names, params, dp);
        records.push_back(makeTeam7EdgeFluxOmegaScalingRecord(frequency_hz, sigma, summary));
    }
    params.frequency_ = saved_frequency;
    printTeam7EdgeFluxOmegaScalingReport(records, Real(50));
    return records;
}

inline Team7NormalizationSweepRecord makeTeam7NormalizationSweepRecord(
    OphelieEdgeFluxNormalizationMode normalization_mode, Real safe_rhs_l2, Real safe_rhs_max_abs,
    const Team7OneWayEdgeFluxSummary &summary, Real phase90_rms_post = 0.0)
{
    Team7NormalizationSweepRecord record;
    record.normalization_mode = ophelieEdgeFluxNormalizationModeName(normalization_mode);
    record.safe_rhs_l2 = safe_rhs_l2;
    record.safe_rhs_max_abs = safe_rhs_max_abs;
    record.input_scale_applied = summary.edge_flux_input_scale;
    record.input_scale_measured = summary.edge_flux_input_scale_measured;
    record.rhs_l2_raw = summary.edge_flux_rhs_l2_pre_norm;
    record.rhs_l2_scaled = summary.edge_flux_rhs_l2_scaled;
    record.j_imag_vol_pre_restore = summary.j_imag_vol_pre_restore;
    record.j_imag_vol_post_restore = summary.j_imag_vol_norm;
    record.b_ind_over_bcoil_pre_restore = summary.b_ind_over_b_coil_pre_restore;
    record.b_ind_over_bcoil_post_restore = summary.b_ind_over_b_coil;
    record.p_recon_pre_restore = summary.joule_power_pre_restore_w;
    record.p_recon_post_restore = summary.joule_power_recon_w;
    record.phase90_rms_post = phase90_rms_post;
    record.restore_invariance_error_j = summary.restore_invariance_error_j;
    record.restore_invariance_error_b = summary.restore_invariance_error_b;
    record.restore_invariance_error_p = summary.restore_invariance_error_p;
    record.j_restore_linearity = summary.j_restore_linearity;
    record.bind_over_bcoil_post = summary.b_ind_over_b_coil;
    return record;
}

inline void writeTeam7NormalizationSweepCsv(const std::string &output_path,
                                            const StdVec<Team7NormalizationSweepRecord> &records)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 normalization sweep CSV: " << output_path << std::endl;
        return;
    }
    out << "normalization_mode,safe_rhs_l2,safe_rhs_max_abs,input_scale_applied,input_scale_measured,"
           "rhs_l2_raw,rhs_l2_scaled,J_imag_vol_pre_restore,J_imag_vol_post_restore,"
           "Bind_over_Bcoil_pre_restore,Bind_over_Bcoil_post_restore,P_recon_pre_restore,P_recon_post_restore,"
           "phase90_RMS_post,restore_invariance_error_J,restore_invariance_error_B,restore_invariance_error_P,"
           "j_restore_linearity,Bind_over_Bcoil_post\n";
    out << std::setprecision(10);
    for (const Team7NormalizationSweepRecord &record : records)
    {
        out << record.normalization_mode << "," << record.safe_rhs_l2 << "," << record.safe_rhs_max_abs << ","
            << record.input_scale_applied << "," << record.input_scale_measured << "," << record.rhs_l2_raw << ","
            << record.rhs_l2_scaled << "," << record.j_imag_vol_pre_restore << "," << record.j_imag_vol_post_restore
            << "," << record.b_ind_over_bcoil_pre_restore << "," << record.b_ind_over_bcoil_post_restore << ","
            << record.p_recon_pre_restore << "," << record.p_recon_post_restore << "," << record.phase90_rms_post
            << "," << record.restore_invariance_error_j << "," << record.restore_invariance_error_b << ","
            << record.restore_invariance_error_p << "," << record.j_restore_linearity << ","
            << record.bind_over_bcoil_post << "\n";
    }
    std::cout << "[ophelie] TEAM7 normalization sweep CSV: " << output_path << std::endl;
}

inline void printTeam7NormalizationSweepReport(const StdVec<Team7NormalizationSweepRecord> &records)
{
    std::cout << "[team7] edge-flux normalization sweep audit (post-restore invariance vs safe_rhs_l2):" << std::endl;
    for (const Team7NormalizationSweepRecord &record : records)
    {
        std::cout << "[team7]   mode=" << record.normalization_mode << " safe_rhs_l2=" << record.safe_rhs_l2
                  << " scale_applied=" << record.input_scale_applied << " scale_measured=" << record.input_scale_measured
                  << " J_post=" << record.j_imag_vol_post_restore << " Bind/B_post=" << record.bind_over_bcoil_post
                  << " phase90_RMS_post=" << record.phase90_rms_post
                  << " inv_err_J=" << record.restore_invariance_error_j
                  << " inv_err_B=" << record.restore_invariance_error_b
                  << " inv_err_P=" << record.restore_invariance_error_p << std::endl;
    }
}

template <class ExecutionPolicy>
inline Real evaluateTeam7Phase90RmsAfterOneWay(
    SolidBody &plate_body, Inner<> &plate_inner, const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
    Real dp, const StdVec<Team7BzProbePoint> &reference_probes, Real plate_z_skin_min_m, Real plate_z_top_m,
    Real plate_z_mid_m, Team7Phase90Convention phase90_convention,
    const StdVec<OphelieCurrentMomentSample> *filament_moments, SolidBody &coil_body,
    const OphelieCoilFieldNames &coil_names)
{
    if (reference_probes.empty())
    {
        return 0.0;
    }
    BaseParticles &particles = plate_body.getBaseParticles();
    const std::string &j_real = getOphelieAIndJRealFieldName(plate_names, params);
    const std::string &j_imag = getOphelieAIndJImagFieldName(plate_names, params);
    const StdVec<Team7ProbeBzDecomposition> decomp = evaluateTeam7ProbeBzDecomposition(
        coil_body.getBaseParticles(), coil_names, particles, plate_names, params, j_real, j_imag, reference_probes,
        plate_z_skin_min_m, plate_z_top_m, plate_z_mid_m, dp, filament_moments);
    StdVec<Real> phase90_sim_mT;
    StdVec<Real> phase90_ref_mT;
    const bool use_neg_ind = team7Phase90ConventionUsesNegInd(phase90_convention);
    team7ProbePhase90SimAndRef(reference_probes, decomp, phase90_sim_mT, phase90_ref_mT, use_neg_ind);
    const Team7BzCompareMetrics phase90_metrics =
        compareTeam7BzAgainstReference(reference_probes, phase90_sim_mT, phase90_ref_mT, Real(1.0));
    (void)plate_inner;
    return phase90_metrics.rms_rel_error;
}

template <class ExecutionPolicy>
inline StdVec<Team7NormalizationSweepRecord> runTeam7EdgeFluxNormalizationSweep(
    SolidBody &plate_body, Inner<> &plate_inner, const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
    Real dp, const StdVec<Real> &safe_rhs_l2_values, const StdVec<Team7BzProbePoint> &reference_probes,
    bool reference_ok, Team7Phase90Convention phase90_convention, Real plate_z_skin_min_m, Real plate_z_top_m,
    Real plate_z_mid_m, const StdVec<OphelieCurrentMomentSample> &filament_moments, bool use_filament_source,
    SolidBody &coil_body, const OphelieCoilFieldNames &coil_names,
    const OphelieTeam7CoilPathPrepareSummary &coil_path_summary)
{
    StdVec<Team7NormalizationSweepRecord> records;
    records.reserve(safe_rhs_l2_values.size());
    const Real saved_safe_rhs_l2 = params.edge_flux_safe_rhs_l2_;
    const OphelieEdgeFluxNormalizationMode saved_mode = params.edge_flux_normalization_mode_;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::FieldScaleRestore;

    for (Real safe_rhs_l2 : safe_rhs_l2_values)
    {
        params.edge_flux_safe_rhs_l2_ = safe_rhs_l2;
        if (use_filament_source)
        {
            applyOphelieTeam7FilamentRacetrackBiotToGlass(plate_body, plate_names, filament_moments, params.mu0_,
                                                          params.softening_length_);
        }
        else
        {
            StateDynamics<ExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> init_coil_source(
                coil_body, coil_names, coil_path_summary.j0);
            init_coil_source.exec();
            StateDynamics<ExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
                plate_body, coil_body, plate_names, coil_names, params);
            compute_biot.exec();
        }
        const Team7OneWayEdgeFluxSummary summary =
            runTeam7ComplexEdgeFluxOneWay<ExecutionPolicy>(plate_body, plate_inner, plate_names, params, dp);
        Real phase90_rms_post = 0.0;
        if (reference_ok)
        {
            const StdVec<OphelieCurrentMomentSample> *filament_ptr =
                use_filament_source && !filament_moments.empty() ? &filament_moments : nullptr;
            phase90_rms_post = evaluateTeam7Phase90RmsAfterOneWay<ExecutionPolicy>(
                plate_body, plate_inner, plate_names, params, dp, reference_probes, plate_z_skin_min_m, plate_z_top_m,
                plate_z_mid_m, phase90_convention, filament_ptr, coil_body, coil_names);
        }
        records.push_back(makeTeam7NormalizationSweepRecord(params.edge_flux_normalization_mode_, safe_rhs_l2,
                                                            params.edge_flux_safe_rhs_max_abs_, summary,
                                                            phase90_rms_post));
    }
    params.edge_flux_safe_rhs_l2_ = saved_safe_rhs_l2;
    params.edge_flux_normalization_mode_ = saved_mode;
    printTeam7NormalizationSweepReport(records);
    return records;
}

inline void writeTeam7SummaryText(const std::string &output_path, const std::string &level_tag,
                                  const std::string &output_tag, const Team7ValidationPassReport &pass_report,
                                  const Team7BzCompareMetrics &coil_metrics,
                                  const Team7BzCompareMetrics &total_phase90_metrics,
                                  const Team7OneWayEdgeFluxSummary &one_way_summary, Real coil_source_scale)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        return;
    }
    out << "level=" << level_tag << " output_tag=" << output_tag << "\n";
    out << "smoke_passed=" << (pass_report.smoke_passed ? 1 : 0)
        << " team7_phase0_source_passed=" << (pass_report.team7_phase0_source_passed ? 1 : 0)
        << " team7_phase90_probe_passed=" << (pass_report.team7_phase90_probe_passed ? 1 : 0)
        << " diagnostic_only=" << (pass_report.diagnostic_only ? 1 : 0)
        << " team7_validation_passed=" << (pass_report.team7_validation_passed ? 1 : 0) << "\n";
    out << "phase90_rms_abs_mT=" << total_phase90_metrics.rms_abs_error_mT
        << " threshold=" << pass_report.phase90_rms_abs_threshold_mT << "\n";
    out << "bz_rms_coil=" << coil_metrics.rms_rel_error << " coil_source_scale=" << coil_source_scale
        << " Bind/Bcoil=" << one_way_summary.b_ind_over_b_coil << "\n";
    std::cout << "[ophelie] TEAM7 summary: " << output_path << std::endl;
}

/** Append one validation run record (human-readable, for regression tracking). */
inline void appendTeam7ValidationRunRecord(const std::string &output_path, const std::string &level_tag,
                                           const std::string &output_tag, const Team7ValidationPassReport &pass_report,
                                           const Team7BzCompareMetrics &coil_metrics,
                                           const Team7BzCompareMetrics &coil_span_metrics,
                                           const Team7BzCompareMetrics &total_metrics,
                                           const Team7BzCompareMetrics &total_phase90_metrics,
                                           const Team7BzCompareMetrics &total_magnitude_metrics,
                                           const Team7BzCompareMetrics &ind_imag_skin_metrics,
                                           const Team7OneWayEdgeFluxSummary &one_way_summary, Real coil_source_scale,
                                           Real edge_fallback_frac)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path, std::ios::app);
    if (!out)
    {
        std::cout << "[ophelie] could not append TEAM7 validation record: " << output_path << std::endl;
        return;
    }
    out << "=== TEAM7 validation " << level_tag;
    if (!output_tag.empty())
    {
        out << " tag=" << output_tag;
    }
    out << " ===\n";
    out << "smoke_passed=" << (pass_report.smoke_passed ? 1 : 0)
        << " team7_phase0_source_passed=" << (pass_report.team7_phase0_source_passed ? 1 : 0)
        << " team7_phase90_probe_passed=" << (pass_report.team7_phase90_probe_passed ? 1 : 0)
        << " diagnostic_only=" << (pass_report.diagnostic_only ? 1 : 0)
        << " team7_validation_passed=" << (pass_report.team7_validation_passed ? 1 : 0) << "\n";
    out << "coil_source_scale=" << coil_source_scale << "\n";
    out << "phase90_rms_abs_mT=" << total_phase90_metrics.rms_abs_error_mT
        << " phase90_threshold_mT=" << pass_report.phase90_rms_abs_threshold_mT << "\n";
    out << "bz_rms_coil=" << coil_metrics.rms_rel_error << " bz_rms_coil_x_span=" << coil_span_metrics.rms_rel_error
        << " bz_rms_total_phase0=" << total_metrics.rms_rel_error
        << " bz_rms_total_phase90=" << total_phase90_metrics.rms_rel_error
        << " bz_rms_total_mag=" << total_magnitude_metrics.rms_rel_error
        << " bz_rms_ind_imag_skin=" << ind_imag_skin_metrics.rms_rel_error << "\n";
    out << "Bind/Bcoil=" << one_way_summary.b_ind_over_b_coil << " Aind/Acoil=" << one_way_summary.a_ind_over_a_coil
        << " edge_flux_input_scale=" << one_way_summary.edge_flux_input_scale
        << " edge_fallback_frac=" << edge_fallback_frac << "\n";
    out << "max_J_imag=" << one_way_summary.max_j_imag << " J_imag_L2_vol=" << one_way_summary.j_imag_vol_norm
        << " P_recon_W=" << one_way_summary.joule_power_recon_w
        << " phi_eq_res_imag=" << one_way_summary.phi_eq_res_vol_imag << "\n\n";
    std::cout << "[ophelie] TEAM7 validation record appended: " << output_path << std::endl;
}

inline void printTeam7OneWayEdgeFluxSummary(const Team7OneWayEdgeFluxSummary &summary)
{
    std::cout << "[ophelie] TEAM7 one-way summary P_coil_only=" << summary.p_complex_coil_only
              << " P_total=" << summary.p_complex_total << " Aind/Acoil=" << summary.a_ind_over_a_coil
              << " Bind/Bcoil=" << summary.b_ind_over_b_coil << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_VALIDATION_H
