#ifndef ELECTROMAGNETIC_OPHELIE_CLI_H
#define ELECTROMAGNETIC_OPHELIE_CLI_H

#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_racetrack_source.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline bool ophelieCliStartsWith(const char *arg, const char *key)
{
    return std::strncmp(arg, key, std::strlen(key)) == 0;
}

inline const char *ophelieCliValueAfter(const char *arg, const char *key)
{
    return arg + std::strlen(key);
}

inline const char *phiLhsOperatorKindName(OpheliePhiLhsOperatorKind kind)
{
    return kind == OpheliePhiLhsOperatorKind::DivSigmaGrad ? "div-sigma-grad" : "legacy-pairwise";
}

inline const char *phiProjectionOperatorKindName(OpheliePhiProjectionOperatorKind kind)
{
    switch (kind)
    {
    case OpheliePhiProjectionOperatorKind::CompatibleDivGrad:
        return "compatible-div-grad";
    case OpheliePhiProjectionOperatorKind::EdgeFlux:
        return "edge-flux";
    default:
        return "div-grad";
    }
}

inline OpheliePhiProjectionOperatorKind parseOpheliePhiProjectionOperatorKind(const std::string &name)
{
    if (name == "compatible-div-grad" || name == "compatible_div_grad")
    {
        return OpheliePhiProjectionOperatorKind::CompatibleDivGrad;
    }
    if (name == "edge-flux" || name == "edge_flux" || name == "legacy-pairwise-flux")
    {
        return OpheliePhiProjectionOperatorKind::EdgeFlux;
    }
    return OpheliePhiProjectionOperatorKind::DivGrad;
}

inline void applyOpheliePhiProjectionOperatorKind(OphelieParameters &params, OpheliePhiProjectionOperatorKind kind)
{
    params.phi_projection_operator_kind_ = kind;
    switch (kind)
    {
    case OpheliePhiProjectionOperatorKind::CompatibleDivGrad:
        params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
        params.phi_compatible_correction_ = true;
        params.phi_gradient_correction_ = true;
        break;
    case OpheliePhiProjectionOperatorKind::EdgeFlux:
        params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
        params.phi_compatible_correction_ = false;
        params.phi_gradient_correction_ = false;
        break;
    default:
        params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
        params.phi_compatible_correction_ = false;
        params.phi_gradient_correction_ = false;
        break;
    }
}

inline OpheliePhiProjectionOperatorKind inferOpheliePhiProjectionOperatorKind(const OphelieParameters &params)
{
    if (params.phi_compatible_correction_ && params.phi_gradient_correction_)
    {
        return OpheliePhiProjectionOperatorKind::CompatibleDivGrad;
    }
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::LegacyPairwise &&
        params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux)
    {
        return OpheliePhiProjectionOperatorKind::EdgeFlux;
    }
    return OpheliePhiProjectionOperatorKind::DivGrad;
}

inline bool opheliePhiProjectionIsProduction(const OphelieParameters &params)
{
    if (params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux)
    {
        return true;
    }
    return inferOpheliePhiProjectionOperatorKind(params) == OpheliePhiProjectionOperatorKind::DivGrad &&
           !params.phi_compatible_correction_ && !params.phi_gradient_correction_;
}

inline const char *opheliePhiProjectionRouteStatusName(const OphelieParameters &params)
{
    return opheliePhiProjectionIsProduction(params) ? "production" : "diagnostic_only";
}

inline void logOpheliePhiProjectionRouteWarnings(const OphelieParameters &params)
{
    const OpheliePhiProjectionOperatorKind kind = inferOpheliePhiProjectionOperatorKind(params);
    if (params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux)
    {
        std::cout << "[ophelie] edge-flux production: --ophelie-current-form=edge-flux "
                     "(pairwise electromotive drop + JouleHeatEdge + edge diagnostics)." << std::endl;
        return;
    }
    if (kind == OpheliePhiProjectionOperatorKind::CompatibleDivGrad)
    {
        std::cout << "[ophelie][warning] compatible-div-grad is experimental diagnostic only; "
                     "not a production OPHELIE phi solver for real Biot A." << std::endl;
    }
    if (kind == OpheliePhiProjectionOperatorKind::EdgeFlux)
    {
        std::cout << "[ophelie][warning] edge-flux (phi-only) changes only the phi equation pairing "
                     "(legacy-pairwise + legacy-flux)." << std::endl;
        std::cout << "[ophelie][warning] E/J/JouleHeat postprocess still uses particle gradient/divergence; "
                     "use --ophelie-current-form=edge-flux for full edge-flux solver." << std::endl;
    }
}

inline const char *phiSolverKindName(OpheliePhiSolverKind kind)
{
    switch (kind)
    {
    case OpheliePhiSolverKind::GMRES:
        return "GMRES";
    case OpheliePhiSolverKind::PCG:
        return "PCG";
    default:
        return "Jacobi";
    }
}

struct OphelieTestCliOptions
{
    bool compare_level0 = false;
    /** TEAM7 integration smoke: do not fail on phi residual / divJ (fields must still be finite and positive). */
    bool ophelie_smoke = false;
    /** Use TEAM7 STL geometry (particle_generation_em / COMSOL mm mesh, scaled to m). */
    bool native_stl = false;
    bool native_small_air_box = true;
    bool native_relax_air = true;
    Real native_dp_mm = 6.0;
    bool sigma_user_set = false;
    bool enable_self_induction = false;
    bool skip_relaxation = false;
    size_t relaxation_steps = 0;
    size_t relaxation_log_every = 50;
    /** Default 100 (official relax tests); --relax-vtp-every=0 disables relax VTP dumps. */
    size_t relaxation_vtp_every = 100;
    /** Non-empty: Reload.xml folder under build bin (default ./reload when cwd is .../bin/). */
    std::string reload_dir;
    bool no_power_scaling = false;
    /** Jacoutot OPHELIE-aligned profile: phi on, no field scaling, optional coil-current calibration. */
    bool literature_mode = false;
    bool literature_calibrate_current = true;
    Real literature_div_j_l2_reduction_min = 1.25;
    bool compare_team7_bz = false;
    bool compare_team7_bz_loop = false;
    bool compare_team7_bz_rect_loop = false;
    std::string team7_reference_dir;
    Real team7_bz_rms_smoke_threshold = 0.5;
    Real team7_coil_turns = 0.0;
    OphelieCoilSourceModel coil_source_model = OphelieCoilSourceModel::VolumeETheta;
    bool coil_source_model_user_set = false;
    Real racetrack_inset_mm = 0.0;
    Real racetrack_z_mm = 99.0;
    Real racetrack_ds_mm = 2.0;
    bool racetrack_sweep = false;
    /** Append P0 phi compatibility / boundary Jn row to CSV (empty = default path). */
    std::string phi_p0_csv_path;
    /** Optional dp override for vector divergence MMS test. */
    Real vector_divergence_mms_dp = -1.0;
    bool vector_divergence_mms_scan = false;
    std::string vector_divergence_mms_csv_path = "output/ophelie_vector_divergence_mms.csv";
    bool phi_lhs_operator_user_set = false;
    bool phi_rhs_operator_user_set = false;
    bool phi_projection_operator_user_set = false;
    bool ophelie_current_form_user_set = false;
};

inline void finalizeOphelieCurrentFormConfiguration(OphelieParameters &params, OphelieTestCliOptions &cli_options)
{
    if (params.edge_flux_complex_ && !cli_options.ophelie_current_form_user_set)
    {
        params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    }
    if (params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux)
    {
        applyOpheliePhiProjectionOperatorKind(params, OpheliePhiProjectionOperatorKind::EdgeFlux);
        params.phi_edge_flux_diagnostics_ = true;
    }
}

inline bool isOphelieTestCommandLineOption(const char *arg)
{
    return std::strcmp(arg, "--no-phi") == 0 || std::strcmp(arg, "--ophelie-compare-level0") == 0 ||
           std::strcmp(arg, "--ophelie-smoke") == 0 || std::strcmp(arg, "--native-stl") == 0 ||
           std::strcmp(arg, "--native-standard-air") == 0 || std::strcmp(arg, "--native-no-air-relax") == 0 ||
           std::strncmp(arg, "--native-dp-mm=", 15) == 0 || std::strcmp(arg, "--no-power-scaling") == 0 ||
           std::strncmp(arg, "--target-power=", 15) == 0 || std::strcmp(arg, "--compare-team7-bz") == 0 ||
           std::strcmp(arg, "--compare-team7-bz-loop") == 0 || std::strcmp(arg, "--compare-team7-bz-rect-loop") == 0 || std::strncmp(arg, "--team7-coil-turns=", 19) == 0 ||
           std::strcmp(arg, "--coil-j-outer-shell") == 0 || std::strncmp(arg, "--coil-j-outer-shell=", 21) == 0 ||
           std::strncmp(arg, "--team7-level=", 14) == 0 ||
           std::strncmp(arg, "--team7-reference-dir=", 22) == 0 ||
           std::strncmp(arg, "--team7-bz-rms-threshold=", 25) == 0 ||
           std::strncmp(arg, "--team7-coil-source-scale=", 26) == 0 ||
           std::strncmp(arg, "--team7-ind-j-post-scale=", 25) == 0 ||
           std::strcmp(arg, "--team7-probe-phase90-neg-ind=1") == 0 ||
           std::strcmp(arg, "--team7-probe-phase90-neg-ind") == 0 ||
           std::strncmp(arg, "--team7-phase90-convention=", 27) == 0 ||
           std::strncmp(arg, "--team7-validation-mode=", 24) == 0 ||
           std::strncmp(arg, "--team7-output-tag=", 19) == 0 ||
           std::strncmp(arg, "--team7-phase90-rms-abs-threshold-mT=", 37) == 0 ||
           std::strncmp(arg, "--team7-frequency=", 18) == 0 ||
           std::strncmp(arg, "--team7-probe-line=", 20) == 0 ||
           std::strcmp(arg, "--team7-omega-sweep=1") == 0 ||
           std::strcmp(arg, "--team7-omega-sweep") == 0 ||
           std::strncmp(arg, "--ophelie-edge-flux-normalization-mode=", 39) == 0 ||
           std::strcmp(arg, "--team7-normalization-sweep=1") == 0 ||
           std::strcmp(arg, "--team7-normalization-sweep") == 0 ||
           std::strcmp(arg, "--team7-operator-audit=1") == 0 ||
           std::strcmp(arg, "--team7-operator-audit") == 0 ||
           std::strcmp(arg, "--team7-boundary-normal-audit=1") == 0 ||
           std::strcmp(arg, "--team7-boundary-normal-audit") == 0 ||
           std::strcmp(arg, "--team7-boundary-consistency-audit=1") == 0 ||
           std::strcmp(arg, "--team7-boundary-consistency-audit") == 0 ||
           std::strcmp(arg, "--team7-l1-source-audit=1") == 0 ||
           std::strcmp(arg, "--team7-l1-source-audit") == 0 ||
           std::strcmp(arg, "--p6b-boundary-sweep=1") == 0 ||
           std::strcmp(arg, "--p6b-boundary-sweep") == 0 ||
           std::strcmp(arg, "--team7-aind-lenz-audit=1") == 0 ||
           std::strcmp(arg, "--team7-aind-lenz-audit") == 0 ||
           std::strcmp(arg, "--self-induction") == 0 || std::strncmp(arg, "--phi-solver=", 13) == 0 ||
           std::strncmp(arg, "--self-induction-relax=", 23) == 0 ||
           std::strncmp(arg, "--self-induction-max-iter=", 26) == 0 ||
           std::strncmp(arg, "--self-induction-phi-tol=", 25) == 0 ||
           std::strncmp(arg, "--team7-picard-relax-aind=", 26) == 0 || std::strcmp(arg, "--skip-relax") == 0 ||
           std::strncmp(arg, "--relax-steps=", 14) == 0 || std::strncmp(arg, "--relax-log-every=", 18) == 0 ||
           std::strncmp(arg, "--relax-vtp-every=", 18) == 0 || std::strncmp(arg, "--reload-dir=", 13) == 0 ||
           std::strncmp(arg, "--sigma=", 8) == 0 || std::strncmp(arg, "--coil-source-model=", 20) == 0 ||
           std::strncmp(arg, "--racetrack-inset-mm=", 21) == 0 || std::strncmp(arg, "--racetrack-z-mm=", 17) == 0 ||
           std::strncmp(arg, "--racetrack-ds-mm=", 18) == 0 || std::strcmp(arg, "--racetrack-sweep") == 0 ||
           std::strcmp(arg, "--literature-mode") == 0 || std::strcmp(arg, "--no-literature-calibrate") == 0 ||
           ophelieCliStartsWith(arg, "--divj-l2-red-min=") || ophelieCliStartsWith(arg, "--phi-gauge-penalty=") ||
           ophelieCliStartsWith(arg, "--phi-lhs-operator=") ||
           ophelieCliStartsWith(arg, "--phi-projection-operator=") ||
           ophelieCliStartsWith(arg, "--phi-gmres-max-outer-iter=") ||
           ophelieCliStartsWith(arg, "--phi-gmres-restart=") || ophelieCliStartsWith(arg, "--phi-gmres-eq-res-tol=") ||
           ophelieCliStartsWith(arg, "--phi-gmres-device-ops=") ||
           ophelieCliStartsWith(arg, "--phi-gmres-device-krylov=") || ophelieCliStartsWith(arg, "--phi-eq-res-gate=") ||
           ophelieCliStartsWith(arg, "--phi-rhs-operator=") ||
           ophelieCliStartsWith(arg, "--phi-rhs-project-zero-mean=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-diagnostics=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-distance-factor=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-mode=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-normal-source=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-grad-neumann=") ||
           ophelieCliStartsWith(arg, "--phi-boundary-lhs-grad-neumann=") ||
           ophelieCliStartsWith(arg, "--phi-gradient-correction=") ||
           ophelieCliStartsWith(arg, "--phi-compatible-correction=") || ophelieCliStartsWith(arg, "--phi-p0-csv=") ||
           ophelieCliStartsWith(arg, "--phi-biot-divergence-diagnostics=") ||
           ophelieCliStartsWith(arg, "--phi-biot-divergence-csv=") ||
           ophelieCliStartsWith(arg, "--phi-edge-flux-diagnostics=") ||
           ophelieCliStartsWith(arg, "--phi-edge-flux-csv=") ||
           ophelieCliStartsWith(arg, "--ophelie-particle-gradient-diagnostics=") ||
           ophelieCliStartsWith(arg, "--ophelie-current-form=") ||
           ophelieCliStartsWith(arg, "--coil-current-scale=") ||
           ophelieCliStartsWith(arg, "--ophelie-use-a-total-for-edge-flux=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-flux-complex=") ||
           ophelieCliStartsWith(arg, "--ophelie-aind-one-way-feedback=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-flux-safe-rhs-l2=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-flux-safe-rhs-max=") ||
           ophelieCliStartsWith(arg, "--ophelie-pair-weight-regularization=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-flux-imag-a-sign=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-recon-boundary-mode=") ||
           ophelieCliStartsWith(arg, "--ophelie-edge-recon-boundary-width-factor=") ||
           ophelieCliStartsWith(arg, "--ophelie-tangent-ls-distance-norm=") ||
           ophelieCliStartsWith(arg, "--vector-divergence-mms-dp=") ||
           ophelieCliStartsWith(arg, "--vector-divergence-mms-csv=") ||
           std::strcmp(arg, "--vector-divergence-mms-scan") == 0;
}

inline void applyOphelieTestCommandLineOption(const char *arg, OphelieParameters &params, OphelieTestCliOptions &cli_options)
{
    if (std::strcmp(arg, "--no-phi") == 0)
    {
        params.enable_phi_correction_ = false;
    }
    else if (std::strcmp(arg, "--ophelie-compare-level0") == 0)
    {
        cli_options.compare_level0 = true;
    }
    else if (std::strcmp(arg, "--ophelie-smoke") == 0)
    {
        cli_options.ophelie_smoke = true;
    }
    else if (std::strcmp(arg, "--native-stl") == 0)
    {
        cli_options.native_stl = true;
        if (!cli_options.coil_source_model_user_set)
        {
            cli_options.coil_source_model = OphelieCoilSourceModel::VolumeRacetrack;
        }
    }
    else if (std::strcmp(arg, "--native-standard-air") == 0)
    {
        cli_options.native_small_air_box = false;
    }
    else if (std::strcmp(arg, "--native-no-air-relax") == 0)
    {
        cli_options.native_relax_air = false;
    }
    else if (std::strncmp(arg, "--native-dp-mm=", 15) == 0)
    {
        cli_options.native_dp_mm = static_cast<Real>(std::atof(arg + 15));
    }
    else if (std::strcmp(arg, "--no-power-scaling") == 0)
    {
        cli_options.no_power_scaling = true;
        params.enable_power_scaling_ = false;
    }
    else if (std::strncmp(arg, "--target-power=", 15) == 0)
    {
        params.target_joule_power_ = static_cast<Real>(std::atof(arg + 15));
        if (params.target_joule_power_ <= TinyReal)
        {
            params.enable_power_scaling_ = false;
        }
    }
    else if (std::strcmp(arg, "--compare-team7-bz") == 0)
    {
        cli_options.compare_team7_bz = true;
    }
    else if (std::strcmp(arg, "--compare-team7-bz-loop") == 0)
    {
        cli_options.compare_team7_bz = true;
        cli_options.compare_team7_bz_loop = true;
    }
    else if (std::strcmp(arg, "--compare-team7-bz-rect-loop") == 0)
    {
        cli_options.compare_team7_bz = true;
        cli_options.compare_team7_bz_rect_loop = true;
        cli_options.coil_source_model = OphelieCoilSourceModel::FilamentRacetrack;
    }
    else if (std::strncmp(arg, "--coil-source-model=", 20) == 0)
    {
        cli_options.coil_source_model_user_set = true;
        cli_options.coil_source_model = parseOphelieCoilSourceModel(std::string(arg + 20));
        if (cli_options.coil_source_model == OphelieCoilSourceModel::FilamentRacetrack)
        {
            cli_options.compare_team7_bz = true;
        }
    }
    else if (std::strncmp(arg, "--racetrack-inset-mm=", 21) == 0)
    {
        cli_options.racetrack_inset_mm = static_cast<Real>(std::atof(arg + 21));
    }
    else if (std::strncmp(arg, "--racetrack-z-mm=", 17) == 0)
    {
        cli_options.racetrack_z_mm = static_cast<Real>(std::atof(arg + 17));
    }
    else if (std::strncmp(arg, "--racetrack-ds-mm=", 18) == 0)
    {
        cli_options.racetrack_ds_mm = static_cast<Real>(std::atof(arg + 18));
    }
    else if (std::strcmp(arg, "--racetrack-sweep") == 0)
    {
        cli_options.racetrack_sweep = true;
        cli_options.compare_team7_bz = true;
        cli_options.coil_source_model = OphelieCoilSourceModel::FilamentRacetrack;
    }
    else if (std::strncmp(arg, "--team7-coil-turns=", 19) == 0)
    {
        cli_options.team7_coil_turns = static_cast<Real>(std::atof(arg + 19));
    }
    else if (std::strcmp(arg, "--coil-j-outer-shell") == 0)
    {
        params.coil_j_outer_shell_only_ = true;
        params.coil_j_outer_shell_radius_fraction_ = 0.85;
    }
    else if (std::strncmp(arg, "--coil-j-outer-shell=", 21) == 0)
    {
        params.coil_j_outer_shell_only_ = true;
        params.coil_j_outer_shell_radius_fraction_ = static_cast<Real>(std::atof(arg + 21));
    }
    else if (std::strncmp(arg, "--team7-reference-dir=", 22) == 0)
    {
        cli_options.team7_reference_dir = std::string(arg + 22);
    }
    else if (std::strncmp(arg, "--team7-bz-rms-threshold=", 25) == 0)
    {
        cli_options.team7_bz_rms_smoke_threshold = static_cast<Real>(std::atof(arg + 25));
    }
    else if (std::strcmp(arg, "--self-induction") == 0)
    {
        cli_options.enable_self_induction = true;
        params.enable_self_induction_ = true;
        params.enable_phi_correction_ = true;
        std::cout << "[ophelie] WARNING: --self-induction is EXPERIMENTAL (Picard; not in literature_passed). "
                     "Do not use for validation reports.\n";
    }
    else if (std::strncmp(arg, "--phi-solver=", 13) == 0)
    {
        const std::string solver_name(arg + 13);
        if (solver_name == "PCG")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::PCG;
        }
        else if (solver_name == "GMRES")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
        }
        else if (solver_name == "Jacobi")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::Jacobi;
        }
    }
    else if (std::strncmp(arg, "--self-induction-relax=", 23) == 0)
    {
        params.self_induction_relaxation_factor_ = static_cast<Real>(std::atof(arg + 23));
    }
    else if (std::strncmp(arg, "--self-induction-max-iter=", 26) == 0)
    {
        params.self_induction_max_iterations_ = static_cast<size_t>(std::atoi(arg + 26));
    }
    else if (std::strncmp(arg, "--self-induction-phi-tol=", 25) == 0)
    {
        params.self_induction_phi_eq_res_tolerance_ = static_cast<Real>(std::atof(arg + 25));
    }
    else if (std::strncmp(arg, "--team7-picard-relax-aind=", 26) == 0)
    {
        params.self_induction_relax_aind_ = std::atoi(arg + 26) != 0;
    }
    else if (std::strcmp(arg, "--skip-relax") == 0)
    {
        cli_options.skip_relaxation = true;
    }
    else if (std::strncmp(arg, "--relax-steps=", 14) == 0)
    {
        cli_options.relaxation_steps = static_cast<size_t>(std::atoi(arg + 14));
    }
    else if (std::strncmp(arg, "--relax-log-every=", 18) == 0)
    {
        cli_options.relaxation_log_every = static_cast<size_t>(std::atoi(arg + 18));
    }
    else if (std::strncmp(arg, "--relax-vtp-every=", 18) == 0)
    {
        cli_options.relaxation_vtp_every = static_cast<size_t>(std::atoi(arg + 18));
    }
    else if (std::strncmp(arg, "--reload-dir=", 13) == 0)
    {
        cli_options.reload_dir = std::string(arg + 13);
    }
    else if (std::strncmp(arg, "--sigma=", 8) == 0)
    {
        params.sigma_glass_ = static_cast<Real>(std::atof(arg + 8));
        cli_options.sigma_user_set = true;
    }
    else if (std::strcmp(arg, "--literature-mode") == 0)
    {
        cli_options.literature_mode = true;
    }
    else if (std::strcmp(arg, "--no-literature-calibrate") == 0)
    {
        cli_options.literature_calibrate_current = false;
    }
    else if (ophelieCliStartsWith(arg, "--divj-l2-red-min="))
    {
        cli_options.literature_div_j_l2_reduction_min =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--divj-l2-red-min=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-gauge-penalty="))
    {
        params.phi_gauge_penalty_ = static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--phi-gauge-penalty=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-projection-operator="))
    {
        const std::string operator_name(ophelieCliValueAfter(arg, "--phi-projection-operator="));
        applyOpheliePhiProjectionOperatorKind(params, parseOpheliePhiProjectionOperatorKind(operator_name));
        cli_options.phi_projection_operator_user_set = true;
        cli_options.phi_lhs_operator_user_set = true;
        cli_options.phi_rhs_operator_user_set = true;
    }
    else if (ophelieCliStartsWith(arg, "--phi-lhs-operator="))
    {
        const std::string operator_name(ophelieCliValueAfter(arg, "--phi-lhs-operator="));
        if (operator_name == "legacy-pairwise" || operator_name == "legacy_pairwise")
        {
            params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
        }
        else if (operator_name == "div-sigma-grad" || operator_name == "div_sigma_grad")
        {
            params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
        }
        cli_options.phi_lhs_operator_user_set = true;
    }
    else if (ophelieCliStartsWith(arg, "--phi-gmres-max-outer-iter="))
    {
        params.phi_gmres_max_outer_iterations_ =
            static_cast<UnsignedInt>(std::atoi(ophelieCliValueAfter(arg, "--phi-gmres-max-outer-iter=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-gmres-restart="))
    {
        params.phi_gmres_restart_dimension_ =
            static_cast<UnsignedInt>(std::atoi(ophelieCliValueAfter(arg, "--phi-gmres-restart=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-gmres-eq-res-tol="))
    {
        params.phi_gmres_eq_res_tolerance_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--phi-gmres-eq-res-tol=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-gmres-device-ops="))
    {
        params.phi_gmres_use_device_vector_ops_ = (std::atoi(ophelieCliValueAfter(arg, "--phi-gmres-device-ops=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-gmres-device-krylov="))
    {
        params.phi_gmres_use_device_krylov_storage_ =
            (std::atoi(ophelieCliValueAfter(arg, "--phi-gmres-device-krylov=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-eq-res-gate="))
    {
        params.phi_eq_res_vol_gate_ = static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--phi-eq-res-gate=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-rhs-operator="))
    {
        params.phi_rhs_operator_kind_ = parseOpheliePhiRhsOperatorKind(std::string(ophelieCliValueAfter(arg, "--phi-rhs-operator=")));
        cli_options.phi_rhs_operator_user_set = true;
    }
    else if (ophelieCliStartsWith(arg, "--phi-rhs-project-zero-mean="))
    {
        params.phi_rhs_project_zero_mean_ = (std::atoi(ophelieCliValueAfter(arg, "--phi-rhs-project-zero-mean=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-diagnostics="))
    {
        params.phi_boundary_diagnostics_ = (std::atoi(ophelieCliValueAfter(arg, "--phi-boundary-diagnostics=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-distance-factor="))
    {
        params.phi_boundary_distance_factor_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--phi-boundary-distance-factor=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-mode="))
    {
        params.phi_boundary_mode_ = parseOpheliePhiBoundaryMode(std::string(ophelieCliValueAfter(arg, "--phi-boundary-mode=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-normal-source="))
    {
        params.phi_boundary_normal_source_ =
            parseOpheliePhiBoundaryNormalSource(std::string(ophelieCliValueAfter(arg, "--phi-boundary-normal-source=")));
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-grad-neumann="))
    {
        params.phi_boundary_grad_neumann_projection_ =
            (std::atoi(ophelieCliValueAfter(arg, "--phi-boundary-grad-neumann=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-boundary-lhs-grad-neumann="))
    {
        params.phi_boundary_lhs_grad_neumann_ =
            (std::atoi(ophelieCliValueAfter(arg, "--phi-boundary-lhs-grad-neumann=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-gradient-correction="))
    {
        params.phi_gradient_correction_ = (std::atoi(ophelieCliValueAfter(arg, "--phi-gradient-correction=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-compatible-correction="))
    {
        params.phi_compatible_correction_ = (std::atoi(ophelieCliValueAfter(arg, "--phi-compatible-correction=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-p0-csv="))
    {
        cli_options.phi_p0_csv_path = std::string(ophelieCliValueAfter(arg, "--phi-p0-csv="));
    }
    else if (ophelieCliStartsWith(arg, "--phi-biot-divergence-diagnostics="))
    {
        params.phi_biot_divergence_diagnostics_ =
            (std::atoi(ophelieCliValueAfter(arg, "--phi-biot-divergence-diagnostics=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-biot-divergence-csv="))
    {
        params.phi_biot_divergence_csv_path_ = std::string(ophelieCliValueAfter(arg, "--phi-biot-divergence-csv="));
    }
    else if (ophelieCliStartsWith(arg, "--phi-edge-flux-diagnostics="))
    {
        params.phi_edge_flux_diagnostics_ =
            (std::atoi(ophelieCliValueAfter(arg, "--phi-edge-flux-diagnostics=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--phi-edge-flux-csv="))
    {
        params.phi_edge_flux_csv_path_ = std::string(ophelieCliValueAfter(arg, "--phi-edge-flux-csv="));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-particle-gradient-diagnostics="))
    {
        params.output_particle_gradient_diagnostics_ =
            (std::atoi(ophelieCliValueAfter(arg, "--ophelie-particle-gradient-diagnostics=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-current-form="))
    {
        params.ophelie_current_form_ =
            parseOphelieCurrentFormKind(std::string(ophelieCliValueAfter(arg, "--ophelie-current-form=")));
        cli_options.ophelie_current_form_user_set = true;
    }
    else if (ophelieCliStartsWith(arg, "--coil-current-scale="))
    {
        params.coil_current_scale_ = static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--coil-current-scale=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-use-a-total-for-edge-flux="))
    {
        params.use_a_total_for_edge_flux_ =
            (std::atoi(ophelieCliValueAfter(arg, "--ophelie-use-a-total-for-edge-flux=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-flux-complex="))
    {
        params.edge_flux_complex_ =
            (std::atoi(ophelieCliValueAfter(arg, "--ophelie-edge-flux-complex=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-aind-one-way-feedback="))
    {
        params.aind_one_way_feedback_resolve_ =
            (std::atoi(ophelieCliValueAfter(arg, "--ophelie-aind-one-way-feedback=")) != 0);
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-flux-normalization-mode="))
    {
        params.edge_flux_normalization_mode_ = parseOphelieEdgeFluxNormalizationMode(
            std::string(ophelieCliValueAfter(arg, "--ophelie-edge-flux-normalization-mode=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-flux-safe-rhs-l2="))
    {
        params.edge_flux_safe_rhs_l2_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--ophelie-edge-flux-safe-rhs-l2=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-flux-safe-rhs-max="))
    {
        params.edge_flux_safe_rhs_max_abs_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--ophelie-edge-flux-safe-rhs-max=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-pair-weight-regularization="))
    {
        params.pair_weight_regularization_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--ophelie-pair-weight-regularization=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-flux-imag-a-sign="))
    {
        params.edge_flux_imag_a_sign_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--ophelie-edge-flux-imag-a-sign=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-recon-boundary-mode="))
    {
        params.edge_recon_boundary_mode_ =
            parseOphelieEdgeReconBoundaryMode(std::string(ophelieCliValueAfter(arg, "--ophelie-edge-recon-boundary-mode=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-edge-recon-boundary-width-factor="))
    {
        params.edge_recon_boundary_width_factor_ =
            static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--ophelie-edge-recon-boundary-width-factor=")));
    }
    else if (ophelieCliStartsWith(arg, "--ophelie-tangent-ls-distance-norm="))
    {
        params.tangent_ls_distance_norm_ = parseOphelieTangentLsDistanceNorm(
            std::string(ophelieCliValueAfter(arg, "--ophelie-tangent-ls-distance-norm=")));
    }
    else if (ophelieCliStartsWith(arg, "--vector-divergence-mms-dp="))
    {
        cli_options.vector_divergence_mms_dp = static_cast<Real>(std::atof(ophelieCliValueAfter(arg, "--vector-divergence-mms-dp=")));
    }
    else if (ophelieCliStartsWith(arg, "--vector-divergence-mms-csv="))
    {
        cli_options.vector_divergence_mms_csv_path =
            std::string(ophelieCliValueAfter(arg, "--vector-divergence-mms-csv="));
    }
    else if (std::strcmp(arg, "--vector-divergence-mms-scan") == 0)
    {
        cli_options.vector_divergence_mms_scan = true;
    }
}

inline StdVec<std::string> filterOphelieTestCommandLine(int ac, char *av[], OphelieParameters &params,
                                                        OphelieTestCliOptions &cli_options)
{
    StdVec<std::string> filtered_arguments;
    filtered_arguments.emplace_back(av[0]);
    for (int arg_index = 1; arg_index < ac; ++arg_index)
    {
        const char *arg = av[arg_index];
        if (isOphelieTestCommandLineOption(arg))
        {
            applyOphelieTestCommandLineOption(arg, params, cli_options);
            continue;
        }
        filtered_arguments.emplace_back(arg);
    }
    finalizeOphelieCurrentFormConfiguration(params, cli_options);
    return filtered_arguments;
}

inline void logOphelieFinalParams(const OphelieParameters &params, const OphelieTestCliOptions &cli_options)
{
    std::cout << "[ophelie] final_params:"
              << " literature_mode=" << (cli_options.literature_mode ? 1 : 0)
              << " target_power=" << params.target_joule_power_
              << " enable_power_scaling=" << (params.enable_power_scaling_ ? 1 : 0)
              << " enable_phi_correction=" << (params.enable_phi_correction_ ? 1 : 0)
              << " phi_gauge_penalty=" << params.phi_gauge_penalty_
              << " phi_solver=" << phiSolverKindName(params.phi_solver_kind_)
              << " phi_gmres_restart_dimension=" << params.phi_gmres_restart_dimension_
              << " phi_gmres_max_outer_iterations=" << params.phi_gmres_max_outer_iterations_
              << " phi_gmres_eq_res_tolerance=" << params.phi_gmres_eq_res_tolerance_
              << " phi_eq_res_gate=" << params.phi_eq_res_vol_gate_
              << " phi_lhs_operator=" << phiLhsOperatorKindName(params.phi_lhs_operator_kind_)
              << " phi_rhs_operator="
              << (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux ? "legacy-flux"
                                                                                         : "div-sigma-a")
              << " phi_boundary_mode=" << phiBoundaryModeName(params.phi_boundary_mode_)
              << " phi_rhs_project_zero_mean=" << (params.phi_rhs_project_zero_mean_ ? 1 : 0)
              << " phi_boundary_grad_neumann=" << (params.phi_boundary_grad_neumann_projection_ ? 1 : 0)
              << " phi_boundary_lhs_grad_neumann=" << (params.phi_boundary_lhs_grad_neumann_ ? 1 : 0)
              << " phi_gradient_correction=" << (params.phi_gradient_correction_ ? 1 : 0)
              << " phi_compatible_correction=" << (params.phi_compatible_correction_ ? 1 : 0)
              << " phi_gmres_device_ops=" << (params.phi_gmres_use_device_vector_ops_ ? 1 : 0)
              << " phi_gmres_device_krylov=" << (params.phi_gmres_use_device_krylov_storage_ ? 1 : 0)
              << " phi_projection_operator=" << phiProjectionOperatorKindName(inferOpheliePhiProjectionOperatorKind(params))
              << " ophelie_current_form=" << ophelieCurrentFormKindName(params.ophelie_current_form_)
              << " edge_flux_complex=" << (params.edge_flux_complex_ ? 1 : 0)
              << " aind_one_way_feedback=" << (params.aind_one_way_feedback_resolve_ ? 1 : 0)
              << " use_a_total_for_edge_flux=" << (params.use_a_total_for_edge_flux_ ? 1 : 0)
              << " projection_operator_status=" << opheliePhiProjectionRouteStatusName(params)
              << std::endl;
    logOpheliePhiProjectionRouteWarnings(params);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLI_H
