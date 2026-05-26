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

struct AphiLhsTermFlags
{
    bool laplace_a = true;
    bool laplace_phi = true;
    bool reaction = true;
    bool grad_phi_coupling = true;
    bool div_sigma_a_coupling = true;
};

struct AphiLhsAssemblyOptions
{
    AphiLhsTermFlags terms{};
    AphiGradPhiCouplingMode grad_phi_mode = AphiGradPhiCouplingMode::PairwiseUncorrected;
    AphiDivCouplingMode div_mode = AphiDivCouplingMode::PairwiseUncorrected;
    Real omega = 0.0;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_COUPLING_MODES_CK_H
