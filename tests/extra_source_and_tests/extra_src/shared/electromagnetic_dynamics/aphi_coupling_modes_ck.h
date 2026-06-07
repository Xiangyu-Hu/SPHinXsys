#ifndef APHI_COUPLING_MODES_CK_H
#define APHI_COUPLING_MODES_CK_H

namespace SPH
{
namespace electromagnetics
{

/** GradPhi coupling: sigma * grad(phi) on the A-equation LHS. */
enum class AphiGradPhiCouplingMode
{
    /** DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — not the OPERATOR default. */
    LinearGradientWithBCorrection,
    /** OPERATOR default for matrix-free K(X) and fused AphiApplyCK. Uncorrected pairwise g_ij. */
    PairwiseUncorrected
};

/** DivSigmaA coupling: omega * div(sigma * A) on the phi-equation LHS. */
enum class AphiDivCouplingMode
{
    /** DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — constant sigma only. */
    LinearGradientTraceWithBCorrection,
    /** OPERATOR default — sigma_ij harmonic mean, uncorrected pairwise g_ij. */
    PairwiseUncorrected
};

/** A-divergence penalty: LhsA -= lambda * grad(div A). */
enum class AphiADivergencePenaltyMode
{
    /** DIAGNOSTIC_B_CORRECTED — gradA trace + B-corrected grad(divA). */
    BCorrectedTrace,
    /** OPERATOR default — pairwise divA + pairwise grad(divA); matches JacobiGradDivABlock PC. */
    PairwiseUncorrected
};

/** Contact A-penalty divA/grad(divA) stencil (Stage 10.12-C4). Main K(X) apply stays Inner+Contact. */
enum class AphiContactADivergencePenaltyStencilMode
{
    /** Legacy / C1–C2 equivalence: penalty divA/grad(divA) uses Inner+Contact neighbors. */
    InnerContact,
    /** Contact research default (C4 validated): penalty uses Inner-only stencil per body. */
    InnerOnly
};

struct AphiLhsTermFlags
{
    bool laplace_a = true;
    bool laplace_phi = true;
    bool reaction = true;
    bool grad_phi_coupling = true;
    bool div_sigma_a_coupling = true;
};

/** Stage 10.8 Inner η_A sweep: stable GMRES interval (Contact baseline still keeps λ_A off). */
struct AphiADivergencePenaltyResearchDefaults
{
    static constexpr Real primary_eta_a = 0.1;
    static constexpr Real eta_a_min = 0.1;
    static constexpr Real eta_a_max = 0.3;
    static constexpr Real eta_a_mid = 0.2;
    static constexpr Real optional_eta_a = 0.2;
    static constexpr Real upper_eta_a = 0.3;
};

struct AphiLhsAssemblyOptions
{
    AphiLhsTermFlags terms{};
    AphiGradPhiCouplingMode grad_phi_mode = AphiGradPhiCouplingMode::PairwiseUncorrected;
    AphiDivCouplingMode div_mode = AphiDivCouplingMode::PairwiseUncorrected;
    Real omega = 0.0;
    /** Numerical phi gauge regularization (not a physical term). Adds penalty*phi to LhsPhi. */
    bool use_phi_gauge_penalty = false;
    Real phi_gauge_penalty = 0.0;
    /** Numerical Coulomb-gauge regularization (not a physical term). LhsA -= penalty*grad(div A). */
    bool use_a_divergence_penalty = false;
    Real a_divergence_penalty = 0.0;
    AphiADivergencePenaltyMode a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    /** Contact-only: penalty divA/grad(divA) stencil. Does not change Laplace/reaction/coupling apply. */
    AphiContactADivergencePenaltyStencilMode contact_a_divergence_penalty_stencil =
        AphiContactADivergencePenaltyStencilMode::InnerOnly;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_COUPLING_MODES_CK_H
