#ifndef APHI_KRYLOV_DIAGNOSTICS_CK_H
#define APHI_KRYLOV_DIAGNOSTICS_CK_H

#include "base_data_type.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{

inline bool aphiIsFinite(const Real value) { return std::isfinite(static_cast<double>(value)); }

inline bool aphiIsNearZeroAbs(const Real value, const Real tolerance) { return std::abs(value) <= tolerance; }

inline bool aphiIsResidualExplosion(const Real relative_residual, const Real explosion_factor)
{
    return !aphiIsFinite(relative_residual) || relative_residual >= explosion_factor;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_KRYLOV_DIAGNOSTICS_CK_H
