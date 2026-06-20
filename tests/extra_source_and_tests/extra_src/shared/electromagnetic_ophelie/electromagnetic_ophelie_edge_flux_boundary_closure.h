#ifndef ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_BOUNDARY_CLOSURE_H
#define ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_BOUNDARY_CLOSURE_H

#include "electromagnetic_ophelie_parameters.h"

#include <Eigen/Dense>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** P5 operator-level boundary closure flags (derived from edge_recon_boundary_mode_). */
struct OphelieEdgeFluxBoundaryClosureFlags
{
    bool ghost_edge_recon = false;
    bool missing_moment_recon = false;
    bool phi_rhs_ghost = false;
    Real boundary_width_m = 0.0;
};

inline OphelieEdgeFluxBoundaryClosureFlags ophelieEdgeFluxBoundaryClosureFlagsFromMode(
    OphelieEdgeReconBoundaryMode mode, Real boundary_width_m)
{
    OphelieEdgeFluxBoundaryClosureFlags flags;
    flags.boundary_width_m = boundary_width_m;
    switch (mode)
    {
    case OphelieEdgeReconBoundaryMode::NoFluxGhostEdge:
        flags.ghost_edge_recon = true;
        break;
    case OphelieEdgeReconBoundaryMode::NoFluxMissingMoment:
        flags.missing_moment_recon = true;
        break;
    case OphelieEdgeReconBoundaryMode::NoFluxPhiRhsGhost:
        flags.phi_rhs_ghost = true;
        break;
    case OphelieEdgeReconBoundaryMode::NoFluxFull:
        flags.ghost_edge_recon = true;
        flags.missing_moment_recon = true;
        flags.phi_rhs_ghost = true;
        break;
    default:
        break;
    }
    return flags;
}

inline bool ophelieEdgeFluxBoundaryClosureActive(const OphelieEdgeFluxBoundaryClosureFlags &flags)
{
    return flags.ghost_edge_recon || flags.missing_moment_recon || flags.phi_rhs_ghost;
}

inline bool ophelieEdgeReconBoundaryModeUsesProductionClosure(OphelieEdgeReconBoundaryMode mode)
{
    return mode == OphelieEdgeReconBoundaryMode::NoFluxGhostEdge ||
           mode == OphelieEdgeReconBoundaryMode::NoFluxMissingMoment ||
           mode == OphelieEdgeReconBoundaryMode::NoFluxPhiRhsGhost ||
           mode == OphelieEdgeReconBoundaryMode::NoFluxFull;
}

/** P5.4: e_ig=0 ghost edge adds w_ig * n n^T to LS matrix M (b unchanged). */
inline void accumulateOphelieNoFluxGhostEdgeLsContribution(
    Eigen::Matrix<double, Dimensions, Dimensions> &m_acc, Eigen::Matrix<double, Dimensions, 1> &b_acc,
    const Vecd &unit_normal, Real ghost_conductance)
{
    const Eigen::Matrix<double, Dimensions, 1> e_hat = unit_normal.template cast<double>();
    const double w_ig = static_cast<double>(std::abs(ghost_conductance));
    m_acc += w_ig * e_hat * e_hat.transpose();
    (void)b_acc;
}

/** P5.5: level-set missing-support moment — boost tangential rank in boundary LS (I - n n^T). */
inline void accumulateOphelieMissingMomentLsCorrection(Eigen::Matrix<double, Dimensions, Dimensions> &m_acc,
                                                       const Vecd &unit_normal, Real moment_weight)
{
    const Eigen::Matrix<double, Dimensions, 1> n = unit_normal.template cast<double>();
    const Eigen::Matrix<double, Dimensions, Dimensions> identity =
        Eigen::Matrix<double, Dimensions, Dimensions>::Identity();
    m_acc += static_cast<double>(moment_weight) * (identity - n * n.transpose());
}

inline Real ophelieBoundaryGhostConductanceFromMaxNeighbor(Real max_neighbor_conductance)
{
    return max_neighbor_conductance;
}

inline Real ophelieBoundaryMissingMomentWeight(Real max_neighbor_conductance, Real ghost_distance_m)
{
    return max_neighbor_conductance * ghost_distance_m * ghost_distance_m * Real(0.5);
}

inline bool ophelieBoundaryShellParticle(Real signed_distance, Real boundary_width_m)
{
    return boundary_width_m > TinyReal && std::abs(signed_distance) <= boundary_width_m;
}

inline bool ophelieTryUnitBoundaryNormal(const Vecd *normal_direction, size_t index_i, Vecd &unit_normal_out)
{
    if (normal_direction == nullptr)
    {
        return false;
    }
    const Real normal_norm = normal_direction[index_i].norm();
    if (normal_norm <= TinyReal)
    {
        return false;
    }
    unit_normal_out = normal_direction[index_i] / normal_norm;
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_BOUNDARY_CLOSURE_H
