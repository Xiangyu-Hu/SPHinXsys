/**
 * @file spdb_slot_hold.h
 * @brief Blind-slot hold using ContactForce + Ren2014 physical increments.
 * @details Lives only under tests/extra_source_and_tests. Does not modify src/.
 *          Contact uses the validated penalty ContactForce path; RepulsionForceCK
 *          is not used on the flush join. Acoustic relaxation freezes inelastic
 *          history; each accepted physical step calls commit once.
 */
#ifndef SPDB_SLOT_HOLD_H
#define SPDB_SLOT_HOLD_H

#include "sphinxsys.h"

namespace SPH
{
int run_slot_hold(Real dp, Real end_time, Real physical_dt, Real stiffness_scale, int hold_steps);
}

#endif // SPDB_SLOT_HOLD_H
