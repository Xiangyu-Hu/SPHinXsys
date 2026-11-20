#ifndef ELECTROMAGNETIC_OPHELIE_PARAMETERS_H
#define ELECTROMAGNETIC_OPHELIE_PARAMETERS_H

#include "sphinxsys.h"

#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OpheliePhiSolverKind
{
    Jacobi,
    GMRES,
    PCG
};

/** Phi equation LHS/RHS discretization. DivSigmaGrad matches GradPhi + DivJ postprocess. */
enum class OpheliePhiLhsOperatorKind
{
    LegacyPairwise,
    DivSigmaGrad
};

/** Phi equation RHS from A_src: div(-omega*sigma*A) vs legacy pairwise flux. */
enum class OpheliePhiRhsOperatorKind
{
    DivSigmaA,
    LegacyFlux
};

/** Unified LHS/RHS projection pairing (see OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md). */
enum class OpheliePhiProjectionOperatorKind
{
    DivGrad,
    CompatibleDivGrad,
    EdgeFlux
};

enum class OpheliePhiBoundaryMode
{
    None,
    OneSidedNeumann,
    VirtualShellDiagnostic
};

/** Doc213: how OneSidedNeumann injects n·J=0 into φ RHS. */
enum class OpheliePhiNeumannRhsKind
{
    /** Legacy shell: amp*(V/d_b)*σ*g_n (under-scaled vs LegacyPairwise; G5b FAIL). */
    ShellVdB,
    /** Formal: amp * a_sign * ω * A_i · Σ C_ij (x_j−x_i) (matches EMF truncation moment). */
    PairwiseTruncationMoment,
    /**
     * Doc213 G12: DivΣA-family surface flux on mask shell:
     *   rhs += amp * σ * g_n * n·(−m),  m=Σ g_ij,  g_n=−a_sign ω n·A
     * Same g_ij as ComputeOphelieVecdDivergenceCK (not C_ij).
     */
    DivSigmaASurfaceFlux
};

/**
 * Doc214 B1b: restrict pairwise-moment shell by analytic |n_z|.
 * Side: |n_z| < thr (cylinder wall); End: |n_z| ≥ thr (caps). All = no extra gate.
 */
enum class OpheliePhiNeumannMomentShellRegion
{
    All = 0,
    Side = 1,
    End = 2
};

inline OpheliePhiNeumannMomentShellRegion parseOpheliePhiNeumannMomentShellRegion(const std::string &name)
{
    if (name == "side" || name == "wall" || name == "lateral")
    {
        return OpheliePhiNeumannMomentShellRegion::Side;
    }
    if (name == "end" || name == "cap" || name == "caps" || name == "ends")
    {
        return OpheliePhiNeumannMomentShellRegion::End;
    }
    return OpheliePhiNeumannMomentShellRegion::All;
}

inline const char *phiNeumannMomentShellRegionName(OpheliePhiNeumannMomentShellRegion region)
{
    switch (region)
    {
    case OpheliePhiNeumannMomentShellRegion::Side:
        return "side";
    case OpheliePhiNeumannMomentShellRegion::End:
        return "end";
    default:
        return "all";
    }
}

/**
 * Doc214 B1c/B1d: how A enters pairwise truncation moment.
 * Full: A·m, m=Σ C Δx (B1 / G9).
 * SideNormalEndFull (B1c): side (n·A)(n·m); end A·m.
 * SideNormalAbsMEndFull (B1d): side (n·A)|m|sign(n·m); end A·m (full-magnitude normal).
 */
enum class OpheliePhiNeumannMomentAForm
{
    Full = 0,
    SideNormalEndFull = 1,
    SideNormalAbsMEndFull = 2
};

inline OpheliePhiNeumannMomentAForm parseOpheliePhiNeumannMomentAForm(const std::string &name)
{
    if (name == "side-normal-end-full" || name == "side_normal_end_full" || name == "b1c" ||
        name == "normal-side")
    {
        return OpheliePhiNeumannMomentAForm::SideNormalEndFull;
    }
    if (name == "side-normal-absm-end-full" || name == "side_normal_absm_end_full" || name == "b1d" ||
        name == "normal-absm-side")
    {
        return OpheliePhiNeumannMomentAForm::SideNormalAbsMEndFull;
    }
    return OpheliePhiNeumannMomentAForm::Full;
}

inline const char *phiNeumannMomentAFormName(OpheliePhiNeumannMomentAForm form)
{
    switch (form)
    {
    case OpheliePhiNeumannMomentAForm::SideNormalEndFull:
        return "side-normal-end-full";
    case OpheliePhiNeumannMomentAForm::SideNormalAbsMEndFull:
        return "side-normal-absm-end-full";
    default:
        return "full";
    }
}

enum class OpheliePhiBoundaryNormalSource
{
    AnalyticBox,
    AnalyticCylinder,
    LevelSet,
    KernelDeficiency
};

/** Postprocess / phi-solve current formulation (see OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md). */
enum class OphelieCurrentFormKind
{
    ParticleGradient,
    EdgeFlux
};

/** Edge-flux RHS normalization strategy (TEAM7 P0 audit). */
enum class OphelieEdgeFluxNormalizationMode
{
    /** No field scaling or restore; RHS may overflow at high sigma / fine dp. */
    Off,
    /** Scale A/E/J fields before solve, restore after (legacy default). */
    FieldScaleRestore,
    /** Target: physical fields at true scale; solver-local RHS scaling only (P0: no field scale/restore). */
    SolverLocal
};

inline const char *ophelieEdgeFluxNormalizationModeName(OphelieEdgeFluxNormalizationMode mode)
{
    switch (mode)
    {
    case OphelieEdgeFluxNormalizationMode::Off:
        return "off";
    case OphelieEdgeFluxNormalizationMode::FieldScaleRestore:
        return "field-scale-restore";
    case OphelieEdgeFluxNormalizationMode::SolverLocal:
        return "solver-local";
    default:
        return "unknown";
    }
}

inline OphelieEdgeFluxNormalizationMode parseOphelieEdgeFluxNormalizationMode(const std::string &name)
{
    if (name == "off" || name == "0")
    {
        return OphelieEdgeFluxNormalizationMode::Off;
    }
    if (name == "solver-local" || name == "solver_local" || name == "solverlocal")
    {
        return OphelieEdgeFluxNormalizationMode::SolverLocal;
    }
    return OphelieEdgeFluxNormalizationMode::FieldScaleRestore;
}

/** Diagnostic-only edge-recon boundary treatment (P5 no-flux). Does not overwrite production fields. */
enum class OphelieEdgeReconBoundaryMode
{
    None,
    /** Post-hoc E <- E - n(n·E); diagnostic-only. */
    ProjectNormal,
    /** Tangent-plane 2×2 LS; diagnostic-only. */
    TangentLs,
    /** Operator-level: add e_ig=0 ghost edge in ReconstructOphelieEdgeFluxElectricCurrentCK (production). */
    NoFluxGhostEdge,
    /** P5.5a: missing tangential moment in boundary LS reconstruction (production). */
    NoFluxMissingMoment,
    /** P5.5b: ghost outward pair in phi-RHS assembly from A (production). */
    NoFluxPhiRhsGhost,
    /** P5.5: ghost edge + missing moment + phi-RHS ghost (production). */
    NoFluxFull
};

inline const char *ophelieEdgeReconBoundaryModeName(OphelieEdgeReconBoundaryMode mode)
{
    switch (mode)
    {
    case OphelieEdgeReconBoundaryMode::ProjectNormal:
        return "project-normal";
    case OphelieEdgeReconBoundaryMode::TangentLs:
        return "tangent-ls";
    case OphelieEdgeReconBoundaryMode::NoFluxGhostEdge:
        return "no-flux-ghost-edge";
    case OphelieEdgeReconBoundaryMode::NoFluxMissingMoment:
        return "no-flux-missing-moment";
    case OphelieEdgeReconBoundaryMode::NoFluxPhiRhsGhost:
        return "no-flux-phi-rhs-ghost";
    case OphelieEdgeReconBoundaryMode::NoFluxFull:
        return "no-flux-full";
    default:
        return "none";
    }
}

inline OphelieEdgeReconBoundaryMode parseOphelieEdgeReconBoundaryMode(const std::string &name)
{
    if (name == "project-normal" || name == "project_normal" || name == "projectnormal")
    {
        return OphelieEdgeReconBoundaryMode::ProjectNormal;
    }
    if (name == "tangent-ls" || name == "tangent_ls" || name == "tangentls")
    {
        return OphelieEdgeReconBoundaryMode::TangentLs;
    }
    if (name == "no-flux-ghost-edge" || name == "no_flux_ghost_edge" || name == "nofluxghostedge")
    {
        return OphelieEdgeReconBoundaryMode::NoFluxGhostEdge;
    }
    if (name == "no-flux-missing-moment" || name == "no_flux_missing_moment" || name == "missing-moment")
    {
        return OphelieEdgeReconBoundaryMode::NoFluxMissingMoment;
    }
    if (name == "no-flux-phi-rhs-ghost" || name == "no_flux_phi_rhs_ghost" || name == "phi-rhs-ghost")
    {
        return OphelieEdgeReconBoundaryMode::NoFluxPhiRhsGhost;
    }
    if (name == "no-flux-full" || name == "no_flux_full" || name == "nofluxfull")
    {
        return OphelieEdgeReconBoundaryMode::NoFluxFull;
    }
    return OphelieEdgeReconBoundaryMode::None;
}

inline bool ophelieEdgeReconBoundaryModeIsDiagnosticOnly(OphelieEdgeReconBoundaryMode mode)
{
    return mode == OphelieEdgeReconBoundaryMode::ProjectNormal || mode == OphelieEdgeReconBoundaryMode::TangentLs;
}

inline bool ophelieEdgeReconBoundaryModeUsesProductionGhostEdge(OphelieEdgeReconBoundaryMode mode)
{
    return mode == OphelieEdgeReconBoundaryMode::NoFluxGhostEdge || mode == OphelieEdgeReconBoundaryMode::NoFluxFull;
}

/** P5-fix: tangent-LS directional_e denominator (diagnostic-only). */
enum class OphelieTangentLsDistanceNorm
{
    /** Legacy: edge_drop / |r_ij| (3D distance). */
    ThreeD,
    /** edge_drop / |r_ij - (r_ij·n)n| (tangent-plane distance). */
    Tangent
};

inline const char *ophelieTangentLsDistanceNormName(OphelieTangentLsDistanceNorm mode)
{
    switch (mode)
    {
    case OphelieTangentLsDistanceNorm::Tangent:
        return "tangent";
    default:
        return "3d";
    }
}

inline OphelieTangentLsDistanceNorm parseOphelieTangentLsDistanceNorm(const std::string &name)
{
    if (name == "tangent" || name == "distance_t" || name == "t")
    {
        return OphelieTangentLsDistanceNorm::Tangent;
    }
    return OphelieTangentLsDistanceNorm::ThreeD;
}

inline bool ophelieEdgeFluxUsesFieldScaleRestore(OphelieEdgeFluxNormalizationMode mode)
{
    return mode == OphelieEdgeFluxNormalizationMode::FieldScaleRestore;
}

inline OphelieCurrentFormKind parseOphelieCurrentFormKind(const std::string &name)
{
    if (name == "edge-flux" || name == "edge_flux" || name == "edgeflux")
    {
        return OphelieCurrentFormKind::EdgeFlux;
    }
    return OphelieCurrentFormKind::ParticleGradient;
}

inline const char *ophelieCurrentFormKindName(OphelieCurrentFormKind kind)
{
    return kind == OphelieCurrentFormKind::EdgeFlux ? "edge-flux" : "particle-gradient";
}

inline OpheliePhiRhsOperatorKind parseOpheliePhiRhsOperatorKind(const std::string &name)
{
    if (name == "legacy-flux" || name == "legacy_flux" || name == "legacy")
    {
        return OpheliePhiRhsOperatorKind::LegacyFlux;
    }
    return OpheliePhiRhsOperatorKind::DivSigmaA;
}

inline OpheliePhiBoundaryMode parseOpheliePhiBoundaryMode(const std::string &name)
{
    if (name == "one-sided-neumann" || name == "one_sided_neumann" || name == "neumann")
    {
        return OpheliePhiBoundaryMode::OneSidedNeumann;
    }
    if (name == "virtual-shell-diagnostic" || name == "virtual_shell" || name == "virtual-shell")
    {
        return OpheliePhiBoundaryMode::VirtualShellDiagnostic;
    }
    return OpheliePhiBoundaryMode::None;
}

inline OpheliePhiNeumannRhsKind parseOpheliePhiNeumannRhsKind(const std::string &name)
{
    if (name == "pairwise-truncation-moment" || name == "pairwise_truncation_moment" ||
        name == "truncation-moment" || name == "truncation_moment" || name == "pairwise-moment" ||
        name == "pairwise_moment")
    {
        return OpheliePhiNeumannRhsKind::PairwiseTruncationMoment;
    }
    if (name == "div-sigma-a-surface-flux" || name == "div_sigma_a_surface_flux" ||
        name == "div-sigma-surface" || name == "div_sigma_surface" || name == "g-deficiency" ||
        name == "g_deficiency")
    {
        return OpheliePhiNeumannRhsKind::DivSigmaASurfaceFlux;
    }
    return OpheliePhiNeumannRhsKind::ShellVdB;
}

inline const char *phiNeumannRhsKindName(OpheliePhiNeumannRhsKind kind)
{
    switch (kind)
    {
    case OpheliePhiNeumannRhsKind::PairwiseTruncationMoment:
        return "pairwise-truncation-moment";
    case OpheliePhiNeumannRhsKind::DivSigmaASurfaceFlux:
        return "div-sigma-a-surface-flux";
    default:
        return "shell-vdb";
    }
}

inline OpheliePhiBoundaryNormalSource parseOpheliePhiBoundaryNormalSource(const std::string &name)
{
    if (name == "analytic-cylinder" || name == "analytic_cylinder" || name == "cylinder")
    {
        return OpheliePhiBoundaryNormalSource::AnalyticCylinder;
    }
    if (name == "level-set" || name == "level_set" || name == "levelset")
    {
        return OpheliePhiBoundaryNormalSource::LevelSet;
    }
    if (name == "kernel-deficiency" || name == "kernel_deficiency")
    {
        return OpheliePhiBoundaryNormalSource::KernelDeficiency;
    }
    return OpheliePhiBoundaryNormalSource::AnalyticBox;
}

inline const char *phiBoundaryModeName(OpheliePhiBoundaryMode mode)
{
    switch (mode)
    {
    case OpheliePhiBoundaryMode::OneSidedNeumann:
        return "one-sided-neumann";
    case OpheliePhiBoundaryMode::VirtualShellDiagnostic:
        return "virtual-shell-diagnostic";
    default:
        return "none";
    }
}

inline const char *phiBoundaryNormalSourceName(OpheliePhiBoundaryNormalSource source)
{
    switch (source)
    {
    case OpheliePhiBoundaryNormalSource::AnalyticCylinder:
        return "analytic-cylinder";
    case OpheliePhiBoundaryNormalSource::LevelSet:
        return "level-set";
    case OpheliePhiBoundaryNormalSource::KernelDeficiency:
        return "kernel-deficiency";
    default:
        return "analytic-box";
    }
}

struct OphelieParameters
{
    /** Glass melt box half-sizes (m). */
    Real glass_halfsize_x_ = 0.325;
    Real glass_halfsize_y_ = 0.25;
    Real glass_halfsize_z_ = 0.325;

    /** Coil shell box: outer half-sizes and inner cutout half-sizes (m). */
    Real coil_outer_halfsize_x_ = 0.475;
    Real coil_outer_halfsize_y_ = 0.275;
    Real coil_outer_halfsize_z_ = 0.475;
    Real coil_inner_halfsize_x_ = 0.38;
    Real coil_inner_halfsize_y_ = 0.26;
    Real coil_inner_halfsize_z_ = 0.38;

    Vecd glass_center_ = Vecd(0.6, 0.5, 0.5);
    Vecd coil_center_ = Vecd(0.6, 0.5, 0.5);

    Real current_amplitude_ = 1.0;
    Real number_of_turns_ = 8.0;
    /** If positive, overrides equivalentCurrentDensity() for coil J0. */
    Real coil_j0_override_ = 0.0;
    /** Native STL: assign J only on outer xy shell (fraction of max radius); rescales J0 via smaller A_cross. */
    bool coil_j_outer_shell_only_ = false;
    Real coil_j_outer_shell_radius_fraction_ = 0.85;
    Real coil_max_xy_radius_m_ = 0.0;
    Real frequency_ = 300.0e3;

    Real sigma_glass_ = 16.0;
    Real target_joule_power_ = 50.0e3;
    /** When false or target_joule_power_<=0, field_scale=power_scale=1 (required for TEAM7 reference compare). */
    bool enable_power_scaling_ = true;

    Real softening_length_ = 1.0e-6;
    Real mu0_ = 4.0 * Pi * 1.0e-7;

    bool enable_phi_correction_ = true;
    OpheliePhiSolverKind phi_solver_kind_ = OpheliePhiSolverKind::PCG;
    OpheliePhiLhsOperatorKind phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
    /** Doc198 production default: continuity-consistent −ω div(σA). LegacyFlux = diagnostic. */
    OpheliePhiRhsOperatorKind phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    Real phi_gauge_penalty_ = 1.0;
    Real pair_weight_regularization_ = 0.01;
    Real phi_jacobi_relaxation_ = 0.85;
    size_t phi_jacobi_max_iterations_ = 4000;
    Real phi_jacobi_tolerance_ = 1.0e-4;
    UnsignedInt phi_gmres_restart_dimension_ = 40;
    UnsignedInt phi_gmres_max_outer_iterations_ = 80;
    Real phi_gmres_tolerance_ = 1.0e-4;
    /** Extra GMRES stop on ||L phi - RHS||_vol / ||RHS||_vol (0 = L2 Krylov residual only). */
    Real phi_gmres_eq_res_tolerance_ = 0.0;
    /** Use SYCL reductions for GMRES volume-weighted dot/norm (experimental). */
    bool phi_gmres_use_device_vector_ops_ = false;
    /** Keep Arnoldi Krylov vectors in device buffers (experimental; requires device_ops). */
    /** Requires phi_gmres_use_device_vector_ops_; enable with --phi-gmres-device-krylov=1 */
    bool phi_gmres_use_device_krylov_storage_ = false;
    /** Literature / DivSigmaGrad acceptance gate on phi_eq_res_vol. */
    Real phi_eq_res_vol_gate_ = 0.65;
    /** Subtract volume-mean from phi RHS before solve (Neumann compatibility probe). */
    bool phi_rhs_project_zero_mean_ = false;
    /** Log boundary Jn / analytic-cylinder normal diagnostics. */
    bool phi_boundary_diagnostics_ = true;
    /** Boundary shell width for Jn diagnostics: width = factor * dp. */
    Real phi_boundary_distance_factor_ = 2.0;
    /** Doc213: multiply Neumann RHS load (shell needs ≫1; pairwise-moment ~1). */
    Real phi_boundary_neumann_rhs_amplification_ = 1.0;
    /** Doc213 G9: PairwiseTruncationMoment is production default; ShellVdB rollback only. */
    OpheliePhiNeumannRhsKind phi_boundary_neumann_rhs_kind_ =
        OpheliePhiNeumannRhsKind::PairwiseTruncationMoment;
    /**
     * Doc213 G11: if true, pairwise moment only where PhiBoundaryMask>0.5 (boundary shell).
     * Global (false) hurts cylinder rotational A (G10); shell-only keeps EC0/EC1 insulation.
     */
    bool phi_boundary_neumann_moment_shell_only_ = true;
    /**
     * Doc214 B1b: further gate shell moment by |n_z| (needs PhiBoundaryNormal).
     * Default All preserves Doc213 G9/G11. Threshold: |n_z| ≥ 0.5 → end cap.
     */
    OpheliePhiNeumannMomentShellRegion phi_boundary_neumann_moment_shell_region_ =
        OpheliePhiNeumannMomentShellRegion::All;
    /** Doc214 B1c: Full (default) or SideNormalEndFull hybrid. */
    OpheliePhiNeumannMomentAForm phi_boundary_neumann_moment_a_form_ = OpheliePhiNeumannMomentAForm::Full;
    OpheliePhiBoundaryMode phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
    OpheliePhiBoundaryNormalSource phi_boundary_normal_source_ = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
    /** Post-solve grad_phi Neumann projection (experimental; can harm divJ on real Biot). */
    bool phi_boundary_grad_neumann_projection_ = false;
    /** Apply Neumann grad projection inside DivSigmaGrad LHS (GMRES matvec) and E/J grad. */
    bool phi_boundary_lhs_grad_neumann_ = false;
    /** Linear-corrected SPH gradient (B-matrix) for grad_phi and DivSigmaGrad LHS. */
    bool phi_gradient_correction_ = false;
    /**
     * Doc215 E1e: after inverting B→C, if |det(C)| > this cap, replace C with I.
     * 0 = off. Evidence (E1d): shell |det(C)|→30 while interior≈1; cap≈2 kills over-amp.
     */
    Real phi_grad_correction_det_max_cap_ = Real(0);
    /** Paired G_c/D_c operators: grad, LHS div, RHS div, E/J grad, divJ diagnostic. */
    bool phi_compatible_correction_ = false;
    /** Log D(σA) vs D_c(σA) on real Biot A_src in French reduced. */
    bool phi_biot_divergence_diagnostics_ = false;
    std::string phi_biot_divergence_csv_path_;
    /** Legacy pairwise edge-current residual (diagnostic; not production). */
    bool phi_edge_flux_diagnostics_ = false;
    std::string phi_edge_flux_csv_path_;
    OpheliePhiProjectionOperatorKind phi_projection_operator_kind_ = OpheliePhiProjectionOperatorKind::DivGrad;
    OphelieCurrentFormKind ophelie_current_form_ = OphelieCurrentFormKind::ParticleGradient;
    /** Stage 1 edge-flux acceptance: edge_res_red_l2 must exceed this. */
    Real edge_res_red_min_ = 100.0;
    Real edge_power_over_recon_soft_min_ = 0.5;
    Real edge_power_over_recon_soft_max_ = 2.0;
    Real edge_recon_condition_threshold_ = 1.0e-6;
    OphelieEdgeFluxNormalizationMode edge_flux_normalization_mode_ =
        OphelieEdgeFluxNormalizationMode::FieldScaleRestore;
    /** Edge-flux input normalize: vol-weighted ||phi_rhs||_2 target (larger → weaker down-scale). */
    Real edge_flux_safe_rhs_l2_ = 1.0e4;
    /** Edge-flux input normalize: p99.5 |phi_rhs| target (larger → weaker down-scale). */
    Real edge_flux_safe_rhs_max_abs_ = 1.2e3;
    /** P2β solver-local: shared RHS scale for phi PCG and edge E reconstruction (physical fields unchanged). */
    Real edge_flux_solver_local_rhs_scale_ = 1.0;
    /** Imag edge-flux chain a_sign in edge_drop / phi RHS (TEAM7 diagnostic: try −1 for phasor). */
    Real edge_flux_imag_a_sign_ = 1.0;
    /** P5 diagnostic: project E_edge onto tangent plane (E <- E - n(n·E)) using SPHinXsys NormalDirection. */
    OphelieEdgeReconBoundaryMode edge_recon_boundary_mode_ = OphelieEdgeReconBoundaryMode::None;
    /** Boundary shell width for edge-recon boundary diagnostic: width = factor * dp. */
    Real edge_recon_boundary_width_factor_ = 2.0;
    /** P5-fix tangent-LS: directional_e uses 3d distance (default) or tangent-plane distance_t. */
    OphelieTangentLsDistanceNorm tangent_ls_distance_norm_ = OphelieTangentLsDistanceNorm::ThreeD;
    /** Edge-flux production: run particle-gradient E/J/Q into *ParticleDiag fields (not primary). */
    bool output_particle_gradient_diagnostics_ = false;
    /** Edge-flux: read A from a_total (A_coil+A_ind) instead of coil-only a_coil when set. */
    bool use_a_total_for_edge_flux_ = false;
    /** Pair flux antisymmetry gate: q_antisym_rel_l2 must stay below this in edge-flux production. */
    Real q_antisym_rel_l2_max_ = 1.0e-5;
    /** Fixed-current scaling multiplier applied to coil current before EM solve (scaling regression). */
    Real coil_current_scale_ = 1.0;
    /** Complex edge-flux: solve both phi_real and phi_imag chains (default off preserves H v4). */
    bool edge_flux_complex_ = false;
    /** Doc213: shared edge EMF imag+real — rollback via --phi-rhs-operator=legacy-flux.
     *  Production default is DivΣA (finalize sets edge=false, div=true). */
    bool doc212_shared_complex_edge_rhs_ = true;
    /** Doc213 production (finalize): unified DivΣA imag+real. Wins over edge when both set. */
    bool doc212_shared_complex_div_sigma_rhs_ = false;
    /** Diagnostic: unified pairwise (A_i−A_j)·g RHS imag+real. */
    bool doc212_shared_complex_pairwise_rhs_ = false;
    /**
     * Doc214 E2 opt-in: K=G^*(σG), b=−G^*(σ a_sign ω A) with exact adjoint of SPH G.
     * Boundary mode stays none; not an amp scan. Disables DivΣA/EMF shared RHS when set.
     */
    bool doc214_e2_variational_ = false;
    /**
     * Doc216 E5 opt-in: replace full-space GᵀσG PCG by P1 energy projection
     *   min_{φ∈span{1,x,y,z}} ||Gφ + a_sign ω A||_σ²
     * Requires doc214_e2_variational_ (same G / recon). Not a Neumann amp scan.
     */
    bool doc216_e5_p1_projection_ = false;
    /**
     * Doc217 D2 opt-in: full-space local quadratic Neumann collocation (host Eigen SparseLU).
     * Same discrete family as Python d2_cpu_local_poly_deepfix.py. Not production default.
     * Mutually exclusive with doc216_e5_p1_projection_. Geometry fields below must be set by the case.
     */
    bool doc217_d2_local_poly_ = false;
    Real doc217_d2_radius_ = Real(0.1);
    Real doc217_d2_half_height_ = Real(0.3);
    Vecd doc217_d2_center_ = Vecd::Zero();
    Real doc217_d2_dp_ = Real(0.02);
    Real doc217_d2_h_factor_ = Real(2.5);
    Real doc217_d2_h_max_factor_ = Real(4.0);
    /**
     * Doc218 wall F evaluation for D2 Neumann RHS n·F:
     *  0 = center_legacy (Doc217: F at particle center; LHS at wall) — diagnostic only
     *  1 = wall_reconstructed (formal: same local fit, evaluate F at true wall point)
     * Formal regression defaults to 1 after A repro; French entries must not expose legacy CLI.
     */
    int doc218_wall_F_mode_ = 1; // Doc218 formal default after A repro
    /**
     * Doc218 boundary compatibility handling for D2:
     *  0 = boundary_mean_legacy (Doc217: subtract arithmetic mean of scaled boundary b)
     *  1 = none (formal: do not alter physical boundary RHS)
     */
    int doc218_compat_mode_ = 1; // Doc218 formal default after A repro
    /** Doc219 B: test-only failure injection (never for production French runs). */
    bool doc219_inject_fit_fail_ = false;
    bool doc219_inject_solve_fail_ = false;
    /** A_ind one-way: after K[J], re-solve complex edge-flux with A_total (Stage 3.5). */
    bool aind_one_way_feedback_resolve_ = true;
    /** Q spatial soft gate: max(Q)/mean(Q) upper bound (edge-flux production sanity). */
    Real q_spatial_max_over_mean_max_ = 1.0e4;
    /** Q spatial soft gate: outer_mean/center_mean must exceed this. */
    Real q_spatial_outer_over_center_min_ = 1.0;
    /** Post-solve Jacobi smoothing after GMRES/PCG (DivSigmaGrad real cases). */
    size_t phi_post_jacobi_refinement_iterations_ = 0;
    Real phi_post_jacobi_relaxation_ = 0.6;
    size_t phi_pcg_max_iterations_ = 6000;
    Real phi_pcg_tolerance_ = 5.0e-4;

    bool enable_self_induction_ = false;
    size_t self_induction_max_iterations_ = 8;
    Real self_induction_j_tolerance_ = 0.05;
    /** Picard outer-loop: phi_eq_res_vol must be below this (joint with J_rel). */
    Real self_induction_phi_eq_res_tolerance_ = 0.01;
    Real self_induction_relaxation_factor_ = 0.15;
    /** P2 Picard: under-relax A_ind (physical scale); if false, legacy J under-relaxation. */
    bool self_induction_relax_aind_ = true;

    Real omega() const { return 2.0 * Pi * frequency_; }

    Real coilCrossSectionArea() const
    {
        return (coil_outer_halfsize_x_ - coil_inner_halfsize_x_) * 2.0 *
               (coil_outer_halfsize_y_ - coil_inner_halfsize_y_) * 2.0;
    }

    Real equivalentCurrentDensity() const
    {
        const Real coil_height = 2.0 * coil_outer_halfsize_y_;
        const Real radial_gap = 2.0 * (coil_outer_halfsize_x_ - coil_inner_halfsize_x_);
        return number_of_turns_ * current_amplitude_ / (radial_gap * coil_height + TinyReal);
    }
};

/** Effective internal RHS scale for solver-local phi PCG and edge reconstruction. */
inline Real ophelieEdgeFluxEffectiveSolverLocalRhsScale(const OphelieParameters &params)
{
    if (params.edge_flux_normalization_mode_ != OphelieEdgeFluxNormalizationMode::SolverLocal)
    {
        return Real(1);
    }
    const Real scale = params.edge_flux_solver_local_rhs_scale_;
    return (scale > TinyReal && scale < Real(1) - TinyReal) ? scale : Real(1);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PARAMETERS_H
