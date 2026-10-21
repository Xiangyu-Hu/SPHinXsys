/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPDX-License-Identifier: Apache-2.0
 * ------------------------------------------------------------------------- */
/**
 * @file  atomic_dissipation_utils.h
 * @brief Shared atomic accumulation helpers for damping / dissipation kernels.
 * @author OpenAI Codex
 */

#ifndef ATOMIC_DISSIPATION_UTILS_H
#define ATOMIC_DISSIPATION_UTILS_H

#include "base_data_type.h"

namespace SPH
{
inline void atomic_add(Real &variable, const Real &value)
{
    AtomicRef<Real> atomic_variable(variable);
    atomic_variable.fetch_add(value);
}

inline void atomic_add(Vecd &variable, const Vecd &value)
{
    for (size_t i = 0; i < Dimensions; ++i)
    {
        AtomicRef<Real> atomic_variable(variable[i]);
        atomic_variable.fetch_add(value[i]);
    }
}
} // namespace SPH

#endif // ATOMIC_DISSIPATION_UTILS_H