/**
 * @file electromagnetic_ophelie_doc201_unified_candidate.h
 * @brief Doc201 A2 frozen candidate: identical LegacyPairwise+LegacyFlux on real/imag;
 *        edge_flux_normalization=Off (no imag-only SolverLocal scale). Named path only —
 *        does not silently change production DivΣA default.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC201_UNIFIED_CANDIDATE_H
#define ELECTROMAGNETIC_OPHELIE_DOC201_UNIFIED_CANDIDATE_H

#include "electromagnetic_ophelie_parameters.h"

#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kOphelieDoc201UnifiedCandidateName = "doc201_unified_legacy_flux_norm_off";

inline void applyOphelieDoc201UnifiedCandidate(OphelieParameters &params)
{
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
    params.phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
    params.phi_rhs_project_zero_mean_ = true;
    params.phi_boundary_lhs_grad_neumann_ = false;
    params.phi_boundary_grad_neumann_projection_ = false;
    params.phi_gauge_penalty_ = Real(1.0e-8);
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
    params.edge_flux_solver_local_rhs_scale_ = Real(1);
    params.edge_flux_complex_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
}

inline void logOphelieDoc201UnifiedCandidate()
{
    std::cout << "[ophelie][doc201] path=" << kOphelieDoc201UnifiedCandidateName
              << " phi_lhs=legacy-pairwise phi_rhs=legacy-flux (imag+real) "
              << "edge_flux_normalization=off boundary=none" << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
