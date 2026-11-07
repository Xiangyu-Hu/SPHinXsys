/**
 * @file ophelie_central_riemann_ck.h
 * @brief Doc203 central-candidate hybrids for AcousticStep (OPHELIE-local, not framework).
 *
 * - OphelieCentralDensJumpCK: AverageP/AverageV (central momentum) + DissipativeUJump (dens)
 *   + DissipativePJump=0 (no momentum AV). Isolates whether C1 NoRiemann blow-up is dens-jump.
 * - DissipativeRiemannSolverCK remains available as full-AV NoLimiter control.
 */
#ifndef OPHELIE_CENTRAL_RIEMANN_CK_H
#define OPHELIE_CENTRAL_RIEMANN_CK_H

#include "riemann_solver_ck.h"
#include "riemann_solver_ck.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Central pressure average + density DissipativeUJump only. */
class OphelieCentralDensJumpCK
    : public RiemannSolver<Base, WeaklyCompressibleFluid, WeaklyCompressibleFluid>
{
    using BaseRS = RiemannSolver<Base, WeaklyCompressibleFluid, WeaklyCompressibleFluid>;

  public:
    typedef WeaklyCompressibleFluid SourceFluid;
    typedef WeaklyCompressibleFluid TargetFluid;

    OphelieCentralDensJumpCK(WeaklyCompressibleFluid &fluid_i, WeaklyCompressibleFluid &fluid_j)
        : BaseRS(fluid_i, fluid_j)
    {
    }

    Real DissipativePJump(UnsignedInt, UnsignedInt, const Real &) const { return Real(0); }

    Real DissipativeUJump(UnsignedInt i, UnsignedInt j, const Real &p_jump) const
    {
        return p_jump * this->InvImpedanceAve(i, j);
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH
#endif // OPHELIE_CENTRAL_RIEMANN_CK_H
