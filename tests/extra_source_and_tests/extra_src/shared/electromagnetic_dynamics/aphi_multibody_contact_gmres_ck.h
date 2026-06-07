#ifndef APHI_MULTIBODY_CONTACT_GMRES_CK_H
#define APHI_MULTIBODY_CONTACT_GMRES_CK_H

#include "electromagnetic_dynamics/aphi_gmres_solver_ck.h"

#include <functional>

namespace SPH
{
namespace electromagnetics
{

/** One body participating in a coupled Contact GMRES solve (shared variable names / workspace layout). */
struct AphiMultiBodyContactEntry
{
    SPHBody &body;
    /** Null when contact dynamics are supplied via apply_z_override / compute_jacobi_override. */
    Contact<> *contact = nullptr;
    /** Uniform inner relation; null when compute_jacobi_override is set (e.g. adaptive air). */
    Inner<> *inner = nullptr;
    std::function<void()> compute_jacobi_override;
    std::function<void(const AphiBlockNames &, const AphiBlockNames &)> apply_z_override;
};

/**
 * Restarted GMRES over a coupled multi-body Contact system.
 * Global inner products / norms aggregate per-body vol-weighted contributions.
 * Default: equal-per-body average (see AphiGMRESSolverOptions::multibody_inner_product).
 * Arnoldi matvec uses full Contact apply (neighbor Krylov blocks are intentional).
 */
template <class ExecutionPolicy>
class AphiMultiBodyContactGMRESSolverCK
{
  public:
    AphiMultiBodyContactGMRESSolverCK(const StdVec<AphiMultiBodyContactEntry> &bodies, const AphiVariableNames &variable_names,
                                      const AphiGMRESWorkspaceNames &workspace_names,
                                      const AphiLhsAssemblyOptions &operator_options,
                                      const AphiGMRESSolverOptions &solver_options);

    AphiGMRESResult solve();

  protected:
    AphiGMRESResult finalizeBreakdown(AphiGMRESResult &result, AphiGMRESBreakdownCode code, Real residual_norm,
                                      Real initial_norm);

    void maybeRecomputeTrueResidual(UnsignedInt arnoldi_step, Real initial_norm, Real &true_relative_residual,
                                    Real &recursive_true_gap, Real &true_residual_norm);

    void recomputeFinalTrueResidual(AphiGMRESResult &result, Real initial_norm);

    Real globalNormSquared(const AphiBlockNames &block_names) const;
    Real globalDotProduct(const AphiBlockNames &block_x, const AphiBlockNames &block_y) const;

    void computeJacobiDiagonalAll();
    void applySolutionToLhsAll();
    void computeRecursiveResidualAll();
    void computeTrueResidualAll();
    void applyFullContactZToSearchAll(const AphiBlockNames &z_block);
    void preconditionVToZAll(const AphiBlockNames &v_block, const AphiBlockNames &z_block);
    void copyVToZAll(const AphiBlockNames &v_block, const AphiBlockNames &z_block);
    void orthogonalizeSearchAll(Real scale_search, Real scale_v, const AphiBlockNames &v_block);
    void scaleCopyAll(const AphiBlockNames &dst_block, Real alpha, const AphiBlockNames &src_block);
    void axpySolutionAll(Real alpha, const AphiBlockNames &z_block);
    void normalizeV0All(Real inv_beta);

    StdVec<AphiMultiBodyContactEntry> bodies_;
    AphiVariableNames names_;
    AphiGMRESWorkspaceNames workspace_;
    AphiLhsAssemblyOptions operator_options_;
    AphiGMRESSolverOptions solver_options_;
};

} // namespace electromagnetics
} // namespace SPH

#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.hpp"

#endif // APHI_MULTIBODY_CONTACT_GMRES_CK_H
