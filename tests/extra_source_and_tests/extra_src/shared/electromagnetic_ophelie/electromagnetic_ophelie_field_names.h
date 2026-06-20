#ifndef ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H
#define ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H

#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieCoilFieldNames
{
    std::string j_src_real = "JSrcReal";
    std::string j_src_imag = "JSrcImag";
};

struct OphelieGlassFieldNames
{
    std::string sigma = "Sigma";

    std::string a_coil_real = "ACoilReal";
    std::string a_coil_imag = "ACoilImag";
    std::string b_coil_real = "BCoilReal";
    std::string b_coil_imag = "BCoilImag";

    std::string a_ind_real = "AIndReal";
    std::string a_ind_imag = "AIndImag";
    std::string b_ind_real = "BIndReal";
    std::string b_ind_imag = "BIndImag";

    /** Working total vector potential: A_coil + A_ind. */
    std::string a_src_real = "ASrcReal";
    std::string a_src_imag = "ASrcImag";
    std::string b_src_real = "BSrcReal";
    std::string b_src_imag = "BSrcImag";

    std::string e_real = "EReal";
    std::string e_imag = "EImag";
    std::string j_real = "JReal";
    std::string j_imag = "JImag";
    std::string joule_heat = "JouleHeat";

    std::string phi_real = "PhiReal";
    std::string grad_phi_real = "GradPhiReal";
    std::string phi_rhs_real = "PhiRhsReal";
    std::string phi_lhs_real = "PhiLhsReal";
    std::string phi_imag = "PhiImag";
    std::string grad_phi_imag = "GradPhiImag";
    std::string phi_rhs_imag = "PhiRhsImag";
    std::string phi_lhs_imag = "PhiLhsImag";
    std::string phi_laplace_diag = "PhiLaplaceDiag";
    /** Diagnostic SPH divergence of JImag (not a governing unknown). */
    std::string div_j_imag = "DivJImag";

    /** SPH pairwise edge-flux current residual sum_j q_ij (imag chain). */
    std::string edge_flux_residual_imag = "EdgeFluxResidualImag";
    /** SPH pairwise edge-flux current residual sum_j q_ij (real chain). */
    std::string edge_flux_residual_real = "EdgeFluxResidualReal";
    /** Graph/Laplace discrete energy density (diagnostic only, not physical Joule heat). */
    std::string joule_heat_edge = "JouleHeatEdgeGraph";
    std::string power_edge_particle = "PowerEdgeGraphParticle";
    /** Per-component Joule heat from edge-reconstructed E_imag: 0.5*sigma*|E_i|². */
    std::string joule_heat_edge_recon_imag = "JouleHeatEdgeReconImag";
    /** Per-component Joule heat from edge-reconstructed E_real: 0.5*sigma*|E_r|². */
    std::string joule_heat_edge_recon_real = "JouleHeatEdgeReconReal";
    /** Complex Joule heat: 0.5*sigma*(|E_r|²+|E_i|²). */
    std::string joule_heat_edge_recon_complex = "JouleHeatEdgeReconComplex";
    std::string e_edge_recon_imag = "EEdgeReconImag";
    std::string j_edge_recon_imag = "JEdgeReconImag";
    std::string e_edge_recon_real = "EEdgeReconReal";
    std::string j_edge_recon_real = "JEdgeReconReal";
    /** Particle-gradient E/J/Q kept for diagnostic comparison in edge-flux mode. */
    std::string e_imag_particle_diag = "EImagParticleDiag";
    std::string j_imag_particle_diag = "JImagParticleDiag";
    std::string joule_heat_particle_diag = "JouleHeatParticleDiag";
    std::string edge_recon_condition = "EdgeReconCondition";
    std::string edge_recon_fallback = "EdgeReconFallback";
    /** Per-particle max |edge_drop_ij| over neighbors (pair-level sign diagnostic). */
    std::string edge_drop_abs_max = "EdgeDropAbsMax";
    /** Per-particle mean edge_drop_ij² over neighbors (pair-level sign diagnostic). */
    std::string edge_drop_sq_mean = "EdgeDropSqMean";
    /** Per-particle max |q_ij+q_ji| over neighbors (pair antisymmetry diagnostic). */
    std::string edge_q_antisym_max = "EdgeQAntisymMax";
    /** Per-particle Σ(q_ij+q_ji)² over neighbors. */
    std::string edge_q_antisym_sq_sum = "EdgeQAntisymSqSum";
    /** Per-particle Σ q_ij² over neighbors (antisymmetry relative scale). */
    std::string edge_q_scale_sq_sum = "EdgeQScaleSqSum";
    /** Per-particle neighbor count used in q antisymmetry reduction. */
    std::string edge_q_neighbor_count = "EdgeQNeighborCount";
    /** Per-particle count of non-finite q_ij/q_ji pairs. */
    std::string edge_q_nonfinite_count = "EdgeQNonfiniteCount";

    /** Boundary shell for one-sided Neumann (0/1 mask, outward normal, target g_n=-ω n·A). */
    std::string phi_boundary_mask = "PhiBoundaryMask";
    std::string phi_boundary_normal = "PhiBoundaryNormal";
    std::string phi_boundary_gn = "PhiBoundaryGn";
    std::string phi_grad_linear_correction = "PhiGradLinearCorrection";

    /** GMRES Krylov workspace / basis / dot scratch (particle fields, device-resident). */
    std::string krylov_workspace = "PhiKrylovWorkspace";
    std::string krylov_basis = "PhiKrylovBasis";
    std::string krylov_scratch_a = "PhiKrylovScratchA";
    std::string krylov_scratch_b = "PhiKrylovScratchB";
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H
